"""Verify one complete refactor test suite, including the saved OpenFOAM field."""
import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    with path.open(newline='') as stream:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def compare_fields(case, reference):
    files = ['0.006/'+name for name in ['U','p','Cs','k','omega','nut','phi','pointDisplacement','polyMesh/points']]
    files += ['postProcessing/sedimentBudget/budget.csv','postProcessing/shipwreckLoads/0/loads.csv']
    return {name:digest(case/name)==digest(reference/name) for name in files}


def impact_rows(case):
    return read_rows(case/'postProcessing/shipwreckImpact/0/impact.csv')


def compare_serial(serial_case, parallel_case):
    serial, parallel = impact_rows(serial_case), impact_rows(parallel_case)
    if [row['time_s'] for row in serial]!=[row['time_s'] for row in parallel]:
        raise ValueError('Serial and parallel output timelines differ')
    measures = {
        'impact_force': [f'impact_F{axis}_N' for axis in 'xyz'],
        'impact_moment': [f'impact_M{axis}_Nm' for axis in 'xyz'],
        'incident_mass_rate': ['incident_kg_s'],
        'dissipated_power': ['dissipated_power_W'],
    }
    differences = {}
    for metric, keys in measures.items():
        peak = max(math.hypot(*(row[key] for key in keys)) for row in parallel)
        error = max(math.hypot(*(a[key]-b[key] for key in keys)) for a,b in zip(serial,parallel))
        differences[metric] = {'max_absolute_difference':error,'parallel_peak':peak,
                               'difference_over_peak':error/max(peak,1e-12)}
    return differences


def field_integral(baseline, case, evidence):
    dict_contents = '''FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object impactIntegral;
}
functions
{
    shellPower
    {
        type volFieldValue;
        libs (fieldFunctionObjects);
        regionType all;
        operation volIntegrate;
        fields (shipwreckImpactPowerDensity);
        writeFields false;
    }
}
'''
    with tempfile.TemporaryDirectory(dir=evidence) as tmp:
        dictionary = Path(tmp)/'impactIntegral'
        dictionary.write_text(dict_contents)
        process = subprocess.run(['bash',str(baseline/'scripts/run.sh'),'postProcess','-case',str(case),
                                  '-dict',str(dictionary),'-fields','(shipwreckImpactPowerDensity)','-latestTime'],
                                 capture_output=True,text=True)
    (evidence/'impact_field_integral.log').write_text(process.stdout+process.stderr)
    process.check_returncode()
    values = re.findall(r'of shipwreckImpactPowerDensity = (\S+)',process.stdout)
    if len(values)!=1:
        raise ValueError('OpenFOAM did not integrate the saved power field')
    return float(values[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,required=True)
    parser.add_argument('--legacy-case',type=Path,required=True)
    parser.add_argument('--integration',type=Path,required=True)
    parser.add_argument('--serial',type=Path,required=True)
    parser.add_argument('--extended',type=Path,required=True)
    parser.add_argument('--evidence',type=Path,required=True)
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    evidence = args.evidence.resolve()
    evidence.mkdir(exist_ok=False,parents=True)
    batches = {name:json.loads((folder/'batch_summary.json').read_text())
               for name,folder in [('integration',args.integration),('serial',args.serial),('extended',args.extended)]}
    if not all(b['status']=='completed' for b in batches.values()):
        raise RuntimeError('Every planned solve must complete before verification')
    regular = args.integration/'hull10_yaw0'
    elastic = args.integration/'hull10_yaw0_elastic'
    field_checks = {
        'impact_disabled':compare_fields(regular,args.integration/'hull10_yaw0_off'),
        'elastic_material':compare_fields(regular,elastic),
        'previous_cpp_solver':compare_fields(regular,args.legacy_case),
    }
    reference_rows, elastic_rows = impact_rows(regular), impact_rows(elastic)
    lossless = all(row['dissipated_power_W']==0 for row in elastic_rows)
    expected_scale = all(math.isclose(4*a['normal_force_sum_N']/3,b['normal_force_sum_N'],rel_tol=1e-12)
                         and math.isclose(a['incident_kg_s'],b['incident_kg_s'],rel_tol=1e-12)
                         for a,b in zip(reference_rows,elastic_rows))
    old_loads = read_rows(regular/'postProcessing/shipwreckLoads/0/loads.csv')
    consistent_sampling = all(math.isclose(a['incident_Cs_kg_s'],b['incident_kg_s'],rel_tol=1e-12)
                              for a,b in zip(old_loads,reference_rows))
    observed_power = field_integral(args.base.resolve(),regular.resolve(),evidence)
    expected_power = reference_rows[-1]['dissipated_power_W']
    field_matches = math.isclose(observed_power,expected_power,rel_tol=1e-10,abs_tol=1e-12)
    mpi = compare_serial(args.serial/'hull10_yaw0',regular)
    dt = {name:batches['extended']['cases'][name]['impact'] for name in ['hull10_ramp','hull10_ramp_half_dt']}
    time_step_differences = {}
    for metric in ['peak_total_impact_force_N','cumulative_impact_impulse_N_s','cumulative_dissipated_energy_J','cumulative_incident_mass_kg']:
        coarse,fine = [data[metric] for data in dt.values()]
        if isinstance(coarse,list):
            difference = math.hypot(*(a-b for a,b in zip(coarse,fine)))
            scale = math.hypot(*fine)
        else:
            difference = abs(coarse-fine)
            scale = abs(fine)
        time_step_differences[metric] = {'absolute_difference':difference,'difference_over_half_step_value':difference/max(scale,1e-12)}
    provenance_files = ['rust/impact.rs','src/ImpactABI.H','src/ImpactDiagnostics.H','src/HullSurface.H','src/Shipwreck.H',
                        'src/prepare.cpp','run.sh','build.sh','env.sh','verify_refactor.py','install_solver.py',
                        'check_run.py','check_impact.py','run_batch.py','build/bin/shipwreckSedExnerFoam',
                        'build/bin/shipwreck-case','build/lib/libshipwreck_impact.a','build/baseline_sources.json']
    report = {'numerical_integration_passed':all(c.get('numerical_checks_passed') for b in batches.values() for c in b['cases'].values()),
              'published_field_equality':field_checks,'elastic_no_dissipation':lossless,
              'elastic_normal_impulse_scaling_passed':expected_scale,'original_incident_rate_equal':consistent_sampling,
              'field_integral_W':observed_power,'ledger_final_power_W':expected_power,'saved_field_integral_matches':field_matches,
              'serial_parallel':mpi,'half_dt_observations':time_step_differences,
              'source_and_build_sha256':{n:digest(project/n) for n in provenance_files},'experimental_validation_passed':False}
    passed = (report['numerical_integration_passed'] and all(all(checks.values()) for checks in field_checks.values())
              and lossless and expected_scale and consistent_sampling and field_matches
              and all(x['difference_over_peak']<1e-6 for x in mpi.values()))
    report['refactor_checks_passed'] = passed
    (evidence/'refactor_checks.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if not passed:
        raise RuntimeError('Refactor checks failed, see refactor_checks.json')
    print(evidence/'refactor_checks.json')


if __name__=='__main__':
    main()
