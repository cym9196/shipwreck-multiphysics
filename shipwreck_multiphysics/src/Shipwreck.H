// Couples the fixed hull and load diagnostics to the OpenFOAM time loop.
#ifndef Shipwreck_H
#define Shipwreck_H

#include "HullSurface.H"
#include "ImpactDiagnostics.H"
#include "OFstream.H"

namespace Foam {

class Shipwreck {
    const fvMesh& mesh_;
    const IOdictionary properties_;
    const vector centre_;
    const scalar thickness_;
    const dimensionedScalar penaltyTime_;
    const HullSurface surface_;
    volScalarField mask_;
    autoPtr<ImpactDiagnostics> impact_;
    scalarField distances_;
    DynamicList<label> impactCells_;
    struct ShieldedPatch {
        bool valid = false;
        DynamicList<label> faces;
    };
    mutable List<ShieldedPatch> shieldedPatches_;
    autoPtr<OFstream> loads_;

public:
    Shipwreck(const fvMesh& mesh, scalar sedimentDensity, scalar diameter)
    : mesh_(mesh),
      properties_(IOobject("shipwreckProperties",mesh.time().constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE)),
      centre_(properties_.get<vector>("centre")),
      thickness_(properties_.get<scalar>("interfaceThickness")),
      penaltyTime_("penaltyTime",dimTime,properties_.get<scalar>("penaltyTime")),
      surface_(mesh,properties_.get<word>("surfaceFile"),thickness_),
      mask_(IOobject("shipwreckMask",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
            mesh,dimensionedScalar(dimless,0),"zeroGradient")
    {
        if (!(penaltyTime_.value() > 0)) {
            FatalErrorInFunction << "A positive penalty time is required" << exit(FatalError);
        }
        if (properties_.getOrDefault<bool>("impactEnabled",false)) {
            impact_.reset(new ImpactDiagnostics(mesh,properties_,centre_,sedimentDensity,diameter));
        }
        update();
        label interior = 0;
        for (const scalar m : mask_.primitiveField()) if (m>.5) ++interior;
        reduce(interior,sumOp<label>());
        if (interior == 0) FatalErrorInFunction << "Hull has no resolved interior cells; refine the mesh" << exit(FatalError);
        Info << "Shipwreck initial interior cell count: " << interior << nl;
        if (Pstream::master()) {
            const auto folder = mesh.time().globalPath()/"postProcessing"/"shipwreckLoads"/mesh.time().timeName();
            mkDir(folder);
            loads_.reset(new OFstream(folder/"loads.csv"));
            loads_().precision(16);
            loads_() << "time_s,inlet_Ux_m_s,inlet_Cs_volume_fraction,"
                     << "water_Fx_N,water_Fy_N,water_Fz_N,water_Mx_Nm,water_My_Nm,water_Mz_Nm,"
                     << "incident_Cs_kg_s,incident_Cs_momentum_x_N,incident_Cs_momentum_y_N,incident_Cs_momentum_z_N,"
                     << "solid_core_Cs_m3" << endl;
        }
    }

    void update() {
        auto distances = surface_.signedDistances(mesh_.C().primitiveField());
        auto& mask = mask_.primitiveFieldRef();
        const bool sampleImpact = impact_.valid();
        impactCells_.clear();
        // HullSurface marks cells outside its search box with GREAT.
        const scalar outsideMask = .5*(1-Foam::tanh(GREAT/thickness_));
        forAll(distances, cell) {
            const scalar distance = distances[cell];
            const scalar value = distance == GREAT
                ? outsideMask : .5*(1-Foam::tanh(distance/thickness_));
            mask[cell] = value;
            // Sample the fluid side of the smoothed interface, in cell order.
            if (sampleImpact && value > .001 && value < .5) impactCells_.append(cell);
        }
        mask_.correctBoundaryConditions();
        shieldedPatches_.setSize(mesh_.boundaryMesh().size());
        for (auto& patch : shieldedPatches_) patch.valid = false;
        if (sampleImpact) distances_.transfer(distances);
    }

    tmp<volScalarField> inversePenalisationTime() const { return mask_/penaltyTime_; }

    void shieldBed(label patch, vectorField& shields) const {
        auto& cached = shieldedPatches_[patch];
        if (!cached.valid) {
            const auto& centres = mesh_.Cf().boundaryField()[patch].primitiveField();
            const auto distances = surface_.signedDistances(centres);
            cached.faces.clear();
            forAll(distances, face) {
                if (distances[face] < -thickness_) cached.faces.append(face);
            }
            cached.valid = true;
        }
        for (const label face : cached.faces) shields[face] = Zero;
    }

    void writeLoads(const volVectorField& U,const volScalarField& Cs,const volVectorField& Ws,
                    scalar rhoWater,scalar rhoSediment) {
        const tmp<volVectorField> maskGradient = fvc::grad(mask_);
        const auto& gradient = maskGradient().primitiveField();
        vector force(Zero), moment(Zero), incidentMomentum(Zero);
        scalar incidentMass = 0, coreSediment = 0;
        forAll(gradient,i) {
            const scalar volume = mesh_.V()[i];
            const auto water = rhoWater*mask_[i]*volume/penaltyTime_.value()*U[i];
            force += water;
            moment += (mesh_.C()[i]-centre_) ^ water;
            if (mask_[i]>.95) coreSediment += Cs[i]*volume;
            // -grad(mask) points into the water. Use only the fluid-side shell.
            // Doubling the Heaviside half-gradient recovers unit surface weight
            // in the continuum limit; the finite-grid result is a diagnostic.
            if (mask_[i]>=.5 || mask_[i]<=.001) continue;
            const vector particleVelocity = U[i]+Ws[i];
            const scalar arrival = max(scalar(0),particleVelocity & gradient[i]);
            const scalar massRate = 2*rhoSediment*Cs[i]*arrival*volume;
            incidentMass += massRate;
            incidentMomentum += massRate*particleVelocity;
        }
        reduce(force,sumOp<vector>());
        reduce(moment,sumOp<vector>());
        reduce(incidentMomentum,sumOp<vector>());
        reduce(incidentMass,sumOp<scalar>());
        reduce(coreSediment,sumOp<scalar>());
        const label inlet = mesh_.boundaryMesh().findPatchID("inlet");
        const auto& area = mesh_.magSf().boundaryField()[inlet];
        const scalar totalArea = gSum(area);
        const scalar meanSpeed = gSum(area*U.boundaryField()[inlet].component(vector::X))/totalArea;
        const scalar meanCs = gSum(area*Cs.boundaryField()[inlet])/totalArea;
        if (Pstream::master()) {
            auto& log = loads_();
            log << mesh_.time().value() << ',' << meanSpeed << ',' << meanCs;
            for (const vector& v : {force,moment}) log << ',' << v.x() << ',' << v.y() << ',' << v.z();
            log << ',' << incidentMass << ',' << incidentMomentum.x() << ',' << incidentMomentum.y()
                << ',' << incidentMomentum.z() << ',' << coreSediment << endl;
        }
        if (impact_.valid()) impact_().sample(Cs,U,Ws,impactCells_,gradient,distances_);
    }
};

} // namespace Foam
#endif
