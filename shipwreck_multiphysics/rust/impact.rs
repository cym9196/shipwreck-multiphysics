//! One-way fixed-wall impact model in SI units. No OpenFOAM dependencies.
#![deny(unsafe_op_in_unsafe_fn)]

type Vector = [f64; 3];

#[repr(C)]
pub struct ImpactSample {
    velocity_m_s: Vector,
    outward_normal: Vector,
    area_m2: f64,
    sediment_density_kg_m3: f64,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct ImpactResult {
    incident_kg_s: f64,
    rebound_velocity_m_s: Vector,
    wall_force_n: Vector,
    incident_power_w: f64,
    rebound_power_w: f64,
    dissipated_power_w: f64,
    normal_force_n: f64,
    tangential_force_n: f64,
}

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: Vector) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}

fn scale(a: Vector, multiplier: f64) -> Vector {
    a.map(|x| x * multiplier)
}

fn add(a: Vector, b: Vector) -> Vector {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn evaluate(sample: &ImpactSample, restitution: f64, friction: f64) -> Option<ImpactResult> {
    let v = sample.velocity_m_s;
    let normal_length = norm(sample.outward_normal);
    if !v.iter().all(|x| x.is_finite()) || !normal_length.is_finite() || normal_length <= 1e-14
        || !sample.area_m2.is_finite() || sample.area_m2 < 0.0
        || !sample.sediment_density_kg_m3.is_finite() || sample.sediment_density_kg_m3 < 0.0
    {
        return None;
    }
    let normal = scale(sample.outward_normal, 1.0 / normal_length);
    let normal_velocity = dot(v, normal);
    if !normal_velocity.is_finite() {
        return None;
    }
    if normal_velocity >= 0.0 || sample.area_m2 == 0.0 || sample.sediment_density_kg_m3 == 0.0 {
        return Some(ImpactResult::default());
    }

    let approach_speed = -normal_velocity;
    let tangent = add(v, scale(normal, approach_speed));
    let tangent_speed = norm(tangent);
    let normal_impulse_per_mass = (1.0 + restitution) * approach_speed;
    // Limit the tangential impulse by both Coulomb friction and sticking.
    let tangent_decrement = (friction * normal_impulse_per_mass).min(tangent_speed);
    let tangent_direction = if tangent_speed > 0.0 {
        scale(tangent, 1.0 / tangent_speed)
    } else {
        [0.0; 3]
    };
    let rebound_tangent = tangent_speed - tangent_decrement;
    let rebound_normal = restitution * approach_speed;
    let rebound = add(scale(normal, rebound_normal), scale(tangent_direction, rebound_tangent));
    let mass = sample.sediment_density_kg_m3 * sample.area_m2 * approach_speed;
    let normal_force = mass * normal_impulse_per_mass;
    let tangential_force = mass * tangent_decrement;
    let force = add(scale(normal, -normal_force), scale(tangent_direction, tangential_force));
    let incident = 0.5 * mass * dot(v, v);
    let rebound_power = 0.5 * mass * dot(rebound, rebound);
    // Difference of squares in each component avoids cancellation near e = 1.
    let dissipated = 0.5 * mass
        * ((1.0 - restitution) * approach_speed * normal_impulse_per_mass
            + tangent_decrement * (tangent_speed + rebound_tangent));
    if !force.iter().chain(rebound.iter()).all(|x| x.is_finite())
        || ![mass, incident, rebound_power, dissipated, normal_force, tangential_force]
            .iter().all(|x| x.is_finite())
    {
        return None;
    }
    Some(ImpactResult {
        incident_kg_s: mass,
        rebound_velocity_m_s: rebound,
        wall_force_n: force,
        incident_power_w: incident,
        rebound_power_w: rebound_power,
        dissipated_power_w: dissipated,
        normal_force_n: normal_force,
        tangential_force_n: tangential_force,
    })
}

/// Returns 0 for success, 1 for invalid parameters/pointers, 2 for invalid sample data.
/// An error invalidates the entire output. No allocation occurs across the ABI.
///
/// # Safety
/// When n > 0 both buffers must be valid, aligned, non-overlapping arrays of n elements.
#[no_mangle]
pub unsafe extern "C" fn shipwreck_impact_batch_v1(
    samples: *const ImpactSample, output: *mut ImpactResult, n: usize,
    restitution: f64, friction: f64,
) -> i32 {
    if !restitution.is_finite() || !(0.0..=1.0).contains(&restitution)
        || !friction.is_finite() || friction < 0.0
    {
        return 1;
    }
    if n == 0 {
        return 0;
    }
    if samples.is_null() || output.is_null() || n > isize::MAX as usize / std::mem::size_of::<ImpactResult>() {
        return 1;
    }
    // Use pointer writes so C callers may provide uninitialized output memory.
    for i in 0..n {
        let sample = unsafe { &*samples.add(i) };
        match evaluate(sample, restitution, friction) {
            Some(result) => unsafe { output.add(i).write(result) },
            None => return 2,
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(velocity: Vector) -> ImpactSample {
        ImpactSample { velocity_m_s: velocity, outward_normal: [2.0, 0.0, 0.0],
                       area_m2: 0.25, sediment_density_kg_m3: 4.0 }
    }

    fn close(value: f64, expected: f64) {
        assert!((value - expected).abs() <= 1e-12 * (1.0 + expected.abs()), "{value} versus {expected}");
    }

    #[test]
    fn head_on_rebounds_and_transfers_momentum() {
        let r = evaluate(&sample([-2.0, 0.0, 0.0]), 0.5, 0.2).unwrap();
        close(r.incident_kg_s, 2.0);
        close(r.rebound_velocity_m_s[0], 1.0);
        close(r.wall_force_n[0], -6.0);
        close(r.incident_power_w, 4.0);
        close(r.rebound_power_w, 1.0);
        close(r.dissipated_power_w, 3.0);
    }

    #[test]
    fn elastic_frictionless_retains_energy() {
        let r = evaluate(&sample([-2.0, 3.0, 4.0]), 1.0, 0.0).unwrap();
        for (actual, expected) in r.rebound_velocity_m_s.iter().zip([2.0, 3.0, 4.0]) { close(*actual, expected); }
        close(r.wall_force_n[0], -8.0);
        close(r.dissipated_power_w, 0.0);
    }

    #[test]
    fn sliding_and_sticking_follow_coulomb_limit() {
        let s = sample([-2.0, 3.0, 4.0]);
        let sliding = evaluate(&s, 0.5, 0.25).unwrap();
        close(sliding.tangential_force_n, 1.5);
        close(sliding.rebound_velocity_m_s[1], 2.55);
        let sticking = evaluate(&s, 0.5, 2.0).unwrap();
        close(sticking.tangential_force_n, 10.0);
        assert_eq!(sticking.rebound_velocity_m_s, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn no_incoming_flux_for_outgoing_tangent_or_clear_water() {
        for v in [[2.0, 1.0, 1.0], [0.0, 1.0, 1.0]] {
            assert_eq!(evaluate(&sample(v), 0.5, 0.5).unwrap().wall_force_n, [0.0; 3]);
        }
        let mut s = sample([-1.0, 0.0, 0.0]);
        s.sediment_density_kg_m3 = 0.0;
        close(evaluate(&s, 0.5, 0.5).unwrap().incident_kg_s, 0.0);
    }

    #[test]
    fn balances_momentum_and_energy_for_rotated_contacts() {
        for index in 1..100 {
            let angle = (index as f64) * 0.231;
            let n = [angle.cos(), angle.sin(), 0.0];
            let mut s = sample(add(scale(n, -3.0), [0.0, 0.0, 2.0]));
            s.outward_normal = n;
            let r = evaluate(&s, index as f64 / 100.0, 0.7).unwrap();
            let momentum_change = scale(add(s.velocity_m_s, scale(r.rebound_velocity_m_s, -1.0)), r.incident_kg_s);
            for (actual, expected) in r.wall_force_n.iter().zip(momentum_change) { close(*actual, expected); }
            close(r.incident_power_w, r.rebound_power_w + r.dissipated_power_w);
            assert!(dot(r.rebound_velocity_m_s, n) >= 0.0 && r.dissipated_power_w >= 0.0);
        }
    }

    #[test]
    fn ffi_rejects_invalid_parameters_and_writes_batch() {
        let mut result = [ImpactResult::default(); 2];
        let mut samples = [sample([-2.0, 0.0, 0.0]), sample([2.0, 0.0, 0.0])];
        unsafe {
            assert_eq!(shipwreck_impact_batch_v1(samples.as_ptr(), result.as_mut_ptr(), 2, 0.5, 0.2), 0);
            close(result[0].wall_force_n[0], -6.0);
            close(result[1].incident_kg_s, 0.0);
            assert_eq!(shipwreck_impact_batch_v1(samples.as_ptr(), result.as_mut_ptr(), 2, 1.1, 0.0), 1);
            samples[0].outward_normal = [0.0; 3];
            assert_eq!(shipwreck_impact_batch_v1(samples.as_ptr(), result.as_mut_ptr(), 2, 0.5, 0.2), 2);
            assert_eq!(shipwreck_impact_batch_v1(std::ptr::null(), std::ptr::null_mut(), 0, 0.5, 0.2), 0);
        }
    }
}
