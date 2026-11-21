"""Check the C++/Rust impact ledger independently of the solver."""
import csv
import math
import re
from check_timeline import matches_timeline, time_directory


def assess_impact(case, model, timeline, log, final_directory, clock):
    root = case/'postProcessing/shipwreckImpact'
    if model['model']=='off':
        return {'enabled': False, 'passed': not any(root.glob('*/impact.csv'))}
    if model['model']!='impulse':
        raise ValueError('Unknown case impact model')
    folder = time_directory(root, 0)
    if folder is None:
        raise ValueError(f'Expected one initial-time impact ledger directory in {root}')
    with (folder/'impact.csv').open(newline='') as stream:
        records = [{name:float(value) for name,value in row.items()} for row in csv.DictReader(stream)]
    if not records or any(not math.isfinite(v) for row in records for v in row.values()):
        raise ValueError('Empty or nonfinite impact ledger')
    deltas = [float(v) for v in re.findall(r'^deltaT = (\S+)\s*$',log,re.M)]
    timing_ok = (matches_timeline(timeline, records, *clock) and len(deltas)==len(records)
                 and all(r['delta_t_s']>0 and math.isclose(r['delta_t_s'],step,rel_tol=1e-9,abs_tol=1e-12)
                         for r,step in zip(records,deltas))
                 and all(math.isclose(b['time_s']-a['time_s'],b['delta_t_s'],rel_tol=1e-9,abs_tol=1e-12)
                         for a,b in zip(records,records[1:])))
    max_momentum_error = max_energy_error = max_coulomb_excess = 0.
    params_ok = nonnegative = True
    for r in records:
        params_ok &= (math.isclose(r['normal_restitution'],model['normal_restitution'],abs_tol=1e-12)
                      and math.isclose(r['friction_coefficient'],model['friction_coefficient'],abs_tol=1e-12))
        nonnegative &= all(r[k]>=0 for k in ['incident_kg_s','equivalent_grains_s','incident_power_W',
                          'rebound_power_W','dissipated_power_W','normal_force_sum_N','tangential_force_sum_N','weighted_area_m2'])
        momentum = [r[f'incident_p{axis}_N']-r[f'rebound_p{axis}_N']-r[f'impact_F{axis}_N'] for axis in 'xyz']
        force_scale = max(1.,r['normal_force_sum_N']+r['tangential_force_sum_N'])
        max_momentum_error = max(max_momentum_error, math.hypot(*momentum)/force_scale)
        incoming = r['incident_power_W']
        energy = incoming-r['rebound_power_W']-r['dissipated_power_W']
        max_energy_error = max(max_energy_error,abs(energy)/max(1.,incoming))
        bound = model['friction_coefficient']*r['normal_force_sum_N']
        max_coulomb_excess = max(max_coulomb_excess,(r['tangential_force_sum_N']-bound)/max(1.,bound))
    power_field = final_directory/'shipwreckImpactPowerDensity' if final_directory is not None else None
    field_exists = power_field is not None and power_field.is_file() and power_field.stat().st_size > 64
    result = {
        'enabled':True,'records':len(records),'timeline_matches_solver':timing_ok,
        'parameters_match':params_ok,'nonnegative_scalar_outputs':nonnegative,
        'max_normalized_momentum_residual':max_momentum_error,
        'max_normalized_energy_residual':max_energy_error,
        'max_normalized_coulomb_excess':max_coulomb_excess,
        'saved_power_density_field':field_exists,
        'peak_total_impact_force_N':max(math.hypot(*(r[f'impact_F{axis}_N'] for axis in 'xyz')) for r in records),
        'peak_dissipated_power_W':max(r['dissipated_power_W'] for r in records),
        'cumulative_dissipated_energy_J':math.fsum(r['dissipated_power_W']*r['delta_t_s'] for r in records),
        'cumulative_incident_mass_kg':math.fsum(r['incident_kg_s']*r['delta_t_s'] for r in records),
        'cumulative_impact_impulse_N_s':[
            math.fsum(r[f'impact_F{axis}_N']*r['delta_t_s'] for r in records) for axis in 'xyz'],
    }
    result['passed'] = (timing_ok and params_ok and nonnegative and field_exists
                        and max(max_momentum_error,max_energy_error,max_coulomb_excess)<1e-10)
    return result
