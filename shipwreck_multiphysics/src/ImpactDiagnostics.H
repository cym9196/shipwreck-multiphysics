#ifndef ShipwreckImpactDiagnostics_H
#define ShipwreckImpactDiagnostics_H

#include "ImpactABI.H"
#include "OFstream.H"
#include <cmath>
#include <vector>

namespace Foam {

class ImpactDiagnostics {
    const fvMesh& mesh_;
    const vector centre_;
    const scalar restitution_;
    const scalar friction_;
    const scalar density_;
    const scalar grainMass_;
    volScalarField powerDensity_;
    autoPtr<OFstream> log_;
    std::vector<ImpactSample> samples_;
    std::vector<ImpactResult> results_;
    DynamicList<label> addressing_;
    DynamicList<point> contactPoints_;

public:
    ImpactDiagnostics(const fvMesh& mesh, const dictionary& properties, const vector& centre,
                      scalar density, scalar diameter)
    : mesh_(mesh), centre_(centre),
      restitution_(properties.getOrDefault<scalar>("normalRestitution",.5)),
      friction_(properties.getOrDefault<scalar>("frictionCoefficient",.3)), density_(density),
      grainMass_(density * std::acos(-1.0) * diameter * diameter * diameter / 6),
      powerDensity_(IOobject("shipwreckImpactPowerDensity",mesh.time().timeName(),mesh,
                             IOobject::NO_READ,IOobject::AUTO_WRITE),
                    mesh,dimensionedScalar(dimMass/dimLength/pow3(dimTime),0),"zeroGradient")
    {
        if (!(restitution_ >= 0 && restitution_ <= 1 && std::isfinite(friction_) && friction_ >= 0
              && grainMass_ > 0 && std::isfinite(grainMass_))) {
            FatalErrorInFunction << "Invalid impact coefficients or grain mass" << exit(FatalError);
        }
        if (Pstream::master()) {
            const auto folder = mesh.time().globalPath()/"postProcessing"/"shipwreckImpact"/mesh.time().timeName();
            mkDir(folder);
            log_.reset(new OFstream(folder/"impact.csv"));
            log_().precision(16);
            log_() << "time_s,delta_t_s,incident_kg_s,equivalent_grains_s,"
                   << "incident_px_N,incident_py_N,incident_pz_N,rebound_px_N,rebound_py_N,rebound_pz_N,"
                   << "impact_Fx_N,impact_Fy_N,impact_Fz_N,impact_Mx_Nm,impact_My_Nm,impact_Mz_Nm,"
                   << "incident_power_W,rebound_power_W,dissipated_power_W,"
                   << "normal_force_sum_N,tangential_force_sum_N,weighted_area_m2,normal_restitution,friction_coefficient" << endl;
        }
    }

    void sample(const volScalarField& concentration, const volVectorField& water, const volVectorField& settling,
                const labelUList& impactCells, const vectorField& gradient, const scalarField& distances) {
        samples_.clear();
        addressing_.clear();
        contactPoints_.clear();
        powerDensity_.primitiveFieldRef() = Zero;
        scalar areaSum = 0;
        for (const label cell : impactCells) {
            const scalar gradientMagnitude = mag(gradient[cell]);
            if (gradientMagnitude <= VSMALL) continue;
            const scalar area = 2 * gradientMagnitude * mesh_.V()[cell];
            const vector velocity = water[cell] + settling[cell];
            // Use one consistent resolved interface for area and outward normal.
            // At a facet junction an STL nearest-triangle normal is discontinuous.
            const vector normal = -gradient[cell]/gradientMagnitude;
            samples_.push_back({{velocity.x(),velocity.y(),velocity.z()}, {normal.x(),normal.y(),normal.z()},
                                area, density_ * concentration[cell]});
            addressing_.append(cell);
            contactPoints_.append(mesh_.C()[cell] - distances[cell]*normal);
            areaSum += area;
        }
        results_.resize(samples_.size());
        const auto code = shipwreck_impact_batch_v1(samples_.data(),results_.data(),samples_.size(),restitution_,friction_);
        if (returnReduceOr(code != 0)) {
            FatalErrorInFunction << "Rust impact evaluation rejected a sample or coefficient. Local code "
                                 << code << exit(FatalError);
        }

        vector incident(Zero), rebound(Zero), force(Zero), moment(Zero);
        scalar massRate=0, incomingEnergy=0, outgoingEnergy=0, dissipated=0, normalForce=0, tangentialForce=0;
        forAll(addressing_, i) {
            const label cell = addressing_[i];
            const auto& result = results_[i];
            const vector reaction(result.force[0],result.force[1],result.force[2]);
            incident += result.massRate * vector(samples_[i].velocity[0],samples_[i].velocity[1],samples_[i].velocity[2]);
            rebound += result.massRate * vector(result.reboundVelocity[0],result.reboundVelocity[1],result.reboundVelocity[2]);
            force += reaction;
            moment += (contactPoints_[i] - centre_) ^ reaction;
            massRate += result.massRate;
            incomingEnergy += result.incidentPower;
            outgoingEnergy += result.reboundPower;
            dissipated += result.dissipatedPower;
            normalForce += result.normalForce;
            tangentialForce += result.tangentialForce;
            powerDensity_[cell] = result.dissipatedPower / mesh_.V()[cell];
        }
        powerDensity_.correctBoundaryConditions();
        reduce(incident,sumOp<vector>());
        reduce(rebound,sumOp<vector>());
        reduce(force,sumOp<vector>());
        reduce(moment,sumOp<vector>());
        for (auto* metric : {&massRate,&incomingEnergy,&outgoingEnergy,&dissipated,&normalForce,&tangentialForce,&areaSum}) {
            reduce(*metric,sumOp<scalar>());
        }
        if (Pstream::master()) {
            auto& log = log_();
            log << mesh_.time().value() << ',' << mesh_.time().deltaTValue() << ',' << massRate << ',' << massRate/grainMass_;
            for (const vector& v : {incident,rebound,force,moment}) log << ',' << v.x() << ',' << v.y() << ',' << v.z();
            log << ',' << incomingEnergy << ',' << outgoingEnergy << ',' << dissipated << ','
                << normalForce << ',' << tangentialForce << ',' << areaSum << ',' << restitution_ << ',' << friction_ << endl;
        }
    }
};

} // namespace Foam
#endif
