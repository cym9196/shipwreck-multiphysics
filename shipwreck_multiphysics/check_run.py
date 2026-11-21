"""Check saved shipwreck fields, observed inlet history, loads and raw inventory."""
import csv
import json
import math
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

from check_impact import assess_impact
from check_timeline import clock_settings, matches_timeline, same_clock, time_directory

case = Path(sys.argv[1]).resolve()
scenario = json.loads((case/'scenario.json').read_text())
end = scenario['end_time_s']
log = (case/'log.solver').read_text()
mesh = (case/'log.checkMesh.final').read_text()

def rows(path):
    with path.open(newline='') as stream:
        records = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    if not records or any(not math.isfinite(v) for row in records for v in row.values()):
        raise ValueError(f'Empty or nonfinite CSV: {path}')
    return records

flow = rows(case/'inflow.csv')
flow_times = [row['time_s'] for row in flow]
if (len(flow_times) < 2 or flow_times[0] != 0 or flow_times[-1] < end
        or any(a >= b for a, b in zip(flow_times, flow_times[1:]))):
    raise ValueError('Inflow times must increase strictly and cover 0 through end time')
load_root = case/'postProcessing/shipwreckLoads'
load_directory = time_directory(load_root, 0)
if load_directory is None:
    raise ValueError(f'Expected one initial-time load ledger directory in {load_root}')
loads = rows(load_directory/'loads.csv')
ledger = rows(case/'postProcessing/sedimentBudget/budget.csv')
timeline = [float(t) for t in re.findall(r'^Time = (\S+)\s*$', log, re.M)]

def inlet_errors():
    bracket = 0
    errors = {'speed_m_s': 0., 'Cs_volume_fraction': 0.}
    for row in loads:
        t = row['time_s']
        while bracket+2 < len(flow) and t > flow[bracket+1]['time_s']:
            bracket += 1
        lo, hi = flow[bracket:bracket+2]
        weight = max(0., min(1., (t-lo['time_s'])/(hi['time_s']-lo['time_s'])))
        for key, observed in [('speed_m_s','inlet_Ux_m_s'), ('Cs_volume_fraction','inlet_Cs_volume_fraction')]:
            expected = lo[key]+weight*(hi[key]-lo[key])
            errors[key] = max(errors[key],abs(row[observed]-expected))
    return errors

inlet = inlet_errors()
clock = clock_settings(case)
at_end = bool(timeline) and same_clock(timeline[-1], end, *clock)
# Check the actual endpoint in the high-precision CSVs as well as its printed name.
end_tolerance = max(2*math.ulp(end), min(1e-10, abs(end)*1e-10))
exact_steps = at_end and all(
    abs(dataset[-1]['time_s']-end) <= end_tolerance and matches_timeline(timeline, dataset, *clock)
    for dataset in (loads, ledger))
ratios = []
for row in ledger:
    scale = max(abs(row['suspended_solid_m3']),
                abs(row['actual_bed_change_solid_m3']+row['pending_point_target_solid_m3']),
                row['cumulative_gross_boundary_solid_m3'])
    if scale > 1024*sys.float_info.epsilon:
        ratios.append(abs(row['target_raw_residual_solid_m3'])/scale)
max_ratio = max(ratios,default=None)
angle = re.search(r'Mesh non-orthogonality Max:\s*(\S+)',mesh)
angle = float(angle.group(1)) if angle else math.nan
if not math.isfinite(angle):
    angle = None
co = [float(v) for v in re.findall(r'^Courant Number mean: \S+ max: (\S+)', log, re.M)]
valid_courant = bool(co) and all(math.isfinite(v) and v >= 0 for v in co)
max_courant = max(co) if valid_courant else None
final = time_directory(case, timeline[-1]) if at_end else None
fields = ['U','p','Cs','k','omega','nut','shipwreckMask','polyMesh/points']
saved = final is not None and all((final/name).is_file() and (final/name).stat().st_size > 64 for name in fields)
native_passed = 'Mesh OK.' in mesh
determinant_count = re.search(r'Cells with small determinant.*?number of cells:\s*(\d+)',mesh)
anisotropy = None
anisotropy_passed = False
if not native_passed and saved and re.search(r'^Failed 1 mesh checks\.',mesh,re.M) and determinant_count:
    attribution = case/'quality_attribution'/final.name/'summary.json'
    if attribution.is_file():
        anisotropy = json.loads(attribution.read_text())
        anisotropy_passed = (int(determinant_count.group(1))==anisotropy['flagged_cells']
                             and anisotropy['unexplained_flags']==0 and anisotropy['threshold_unchanged']==.001)
mesh_gate = ('Cell volumes OK.' in mesh and 'Non-orthogonality check OK.' in mesh
             and 'Face tets OK.' in mesh and angle is not None and angle<65
             and (native_passed or anisotropy_passed))
result = {
    'solver_ended': bool(re.search(r'^End\s*$',log,re.M)),
    'time_steps': len(timeline), 'end_time_s': end, 'load_and_budget_times_match_solver': bool(exact_steps),
    'saved_final_fields': saved, 'max_inlet_history_error': inlet,
    'max_water_force_N': max(math.hypot(*(r[f'water_F{i}_N'] for i in 'xyz')) for r in loads),
    'max_incident_sediment_kg_s': max(r['incident_Cs_kg_s'] for r in loads),
    'max_incident_sediment_momentum_N': max(math.hypot(*(r[f'incident_Cs_momentum_{i}_N'] for i in 'xyz')) for r in loads),
    'max_solid_core_suspended_m3': max(r['solid_core_Cs_m3'] for r in loads),
    'max_Courant': max_courant, 'max_nonorthogonality_deg': angle,
    'native_extended_mesh_passed': native_passed,
    'extended_determinant_attribution': anisotropy,
    'mesh_gate_passed': mesh_gate,
    'max_active_target_raw_inventory_ratio': max_ratio,
    'impact': assess_impact(case, scenario.get('impact', {'model': 'off'}), timeline, log, final, clock),
    'experimental_validation_passed': False,
}
result['numerical_checks_passed'] = bool(
    result['solver_ended'] and exact_steps and saved and max(inlet.values())<1e-8
    and max_courant is not None and max_courant < 1 and mesh_gate and result['impact']['passed']
    and (max_ratio is None or max_ratio<.005))
encoded = json.dumps(result, indent=2, allow_nan=False)
with TemporaryDirectory(dir=case, prefix='.run_summary.') as staging:
    temporary = Path(staging)/'run_summary.json'
    temporary.write_text(encoded + '\n', encoding='utf-8')
    temporary.replace(case/'run_summary.json')
print(encoded)
if not result['numerical_checks_passed']:
    sys.exit('Saved results did not meet numerical checks')
