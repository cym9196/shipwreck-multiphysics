#ifndef ShipwreckImpactABI_H
#define ShipwreckImpactABI_H

#include <cstddef>
#include <cstdint>

// Keep field order and SI units synchronized with rust/impact.rs.
struct ImpactSample {
    double velocity[3];
    double outwardNormal[3];
    double area;
    double sedimentMassDensity;
};

struct ImpactResult {
    double massRate;
    double reboundVelocity[3];
    double force[3];
    double incidentPower;
    double reboundPower;
    double dissipatedPower;
    double normalForce;
    double tangentialForce;
};

static_assert(sizeof(ImpactSample) == 8 * sizeof(double));
static_assert(sizeof(ImpactResult) == 12 * sizeof(double));

extern "C" std::int32_t shipwreck_impact_batch_v1(
    const ImpactSample*, ImpactResult*, std::size_t, double restitution, double friction);

#endif
