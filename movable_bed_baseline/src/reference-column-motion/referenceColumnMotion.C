// Reference-coordinate column interpolation for the vertically structured bed.
// Unlike a walk on the current geometry, a held boundary target is idempotent.
#include "displacementMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedValuePointPatchFields.H"
#include "pointConstraints.H"
#include "HashTable.H"
#include "Pstream.H"
#include <cmath>

namespace Foam
{
class referenceColumnMotion : public displacementMotionSolver
{
    label bedPatch_;
    scalar top_;
    labelList bedIndex_;
    scalarField weight_;
    word key(const point& p) const
    {
        return Foam::name(label(std::llround(p.x()*1e8)))+"_"
            +Foam::name(label(std::llround(p.y()*1e8)));
    }
public:
    TypeName("referenceColumnMotion");
    referenceColumnMotion(const polyMesh& mesh,const IOdictionary& dict)
    : displacementMotionSolver(mesh,dict,typeName),
      bedPatch_(mesh.boundaryMesh().findPatchID("bed")),
      top_(coeffDict().getOrDefault<scalar>("fixedTop",.2)),
      bedIndex_(mesh.nPoints(),-1),weight_(mesh.nPoints(),0)
    {
        if(bedPatch_<0 || top_<=0)FatalErrorInFunction<<"Invalid bed/top"<<exit(FatalError);
        const labelList& labels=mesh.boundaryMesh()[bedPatch_].meshPoints();
        List<pointField> allFeet(Pstream::nProcs());
        allFeet[Pstream::myProcNo()]=pointField(points0(),labels);
        Pstream::gatherList(allFeet);Pstream::scatterList(allFeet);
        label count=0;for(const auto& feet:allFeet)count+=feet.size();
        pointField feet(count);label cursor=0;
        for(const auto& chunk:allFeet)for(const point& p:chunk)feet[cursor++]=p;
        HashTable<label> columns(2*count);
        forAll(feet,i)
        {
            const point& p=feet[i];
            if(mag(p.z())>1e-8)
                FatalErrorInFunction<<"Bed is not a flat reference surface"<<exit(FatalError);
            // Processor boundaries duplicate physical bed points. Choose the
            // first owner deterministically and verify targets in solve().
            if(!columns.found(key(p)))columns.insert(key(p),i);
        }
        forAll(bedIndex_,i)
        {
            const point& p=points0()[i];
            if(p.z()>=top_-1e-10)continue;
            auto found=columns.find(key(p));
            if(found==columns.end())
                FatalErrorInFunction<<"No global bed column at "<<p<<exit(FatalError);
            const label b=found();
            const point& foot=feet[b];
            if(mag(p.x()-foot.x())>1e-8 || mag(p.y()-foot.y())>1e-8)
                FatalErrorInFunction<<"Column coordinate mismatch"<<exit(FatalError);
            bedIndex_[i]=b;weight_[i]=max(scalar(0),min(scalar(1),1-p.z()/top_));
        }
    }
    virtual tmp<pointField> curPoints() const override
    {
        return tmp<pointField>::New(points0()+pointDisplacement_.primitiveField());
    }
    virtual void solve() override
    {
        // Read the prescribed fixed-value bed target directly. All point
        // constraints are imposed below in reference coordinates; evaluating
        // current-normal slip fields here is unnecessary and reintroduces the
        // very projection feedback this solver avoids.
        const vectorField localTarget(refCast<const fixedValuePointPatchVectorField>
            (pointDisplacement_.boundaryField()[bedPatch_]));
        List<vectorField> allTargets(Pstream::nProcs());
        allTargets[Pstream::myProcNo()]=localTarget;
        Pstream::gatherList(allTargets);Pstream::scatterList(allTargets);
        label count=0;for(const auto& chunk:allTargets)count+=chunk.size();
        vectorField target(count);label cursor=0;
        for(const auto& chunk:allTargets)for(const vector& v:chunk)target[cursor++]=v;
        const labelList& bedLabels=mesh().boundaryMesh()[bedPatch_].meshPoints();
        forAll(bedLabels,i)
            if(mag(localTarget[i]-target[bedIndex_[bedLabels[i]]])>1e-12)
                FatalErrorInFunction<<"Conflicting shared bed target"<<exit(FatalError);
        forAll(bedIndex_,i)
        {
            const label b=bedIndex_[i];
            pointDisplacement_[i]=b<0?vector::zero:vector(0,0,weight_[i]*target[b].z());
        }
        // The prescribed motion already satisfies the reference geometry:
        // x/y fixed, z=top fixed, vertical wall tangency. Re-projecting onto
        // normals calculated from the current wall introduces lateral roundoff
        // and feeds it back at the next update. Keep the exact vertical map.
        forAll(pointDisplacement_.boundaryField(),patchi)
        {
            if(patchi==bedPatch_)continue;
            if(isA<fixedValuePointPatchVectorField>(pointDisplacement_.boundaryField()[patchi]))
            {
                const labelList& ids=mesh().boundaryMesh()[patchi].meshPoints();
                const auto& prescribed=refCast<const fixedValuePointPatchVectorField>
                    (pointDisplacement_.boundaryField()[patchi]);
                forAll(ids,j)
                    if(mag(pointDisplacement_[ids[j]]-prescribed[j])>1e-12)
                        FatalErrorInFunction<<"Vertical map conflicts with fixed patch "
                            <<mesh().boundaryMesh()[patchi].name()<<exit(FatalError);
            }
        }
    }
};
defineTypeNameAndDebug(referenceColumnMotion,0);
addToRunTimeSelectionTable(motionSolver,referenceColumnMotion,dictionary);
}
