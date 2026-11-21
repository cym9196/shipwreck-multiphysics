// Isolate dynamic-mesh interpolation under a prescribed, smooth real bed target.
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include <cmath>

int main(int argc,char *argv[])
{
    argList::addNote("Mesh-only replay: scale the saved bed target linearly to t/0.2");
    argList::addBoolOption("held","Keep the target constant and verify every point update");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    pointVectorField& displacement=const_cast<pointVectorField&>
        (mesh.lookupObject<pointVectorField>("pointDisplacement"));
    const label patch=mesh.boundaryMesh().findPatchID("bed");
    if(patch<0)FatalErrorInFunction<<"Missing bed patch"<<exit(FatalError);
    const vectorField target(refCast<const fixedValuePointPatchVectorField>
        (displacement.boundaryField()[patch]));
    scalar worstCo=0,worstStep=0;
    const bool held=args.found("held");
    scalar minimum=GREAT;
    // A held-target experiment starts at that prescribed geometry. Its initial
    // placement is not a physical .001 s motion from the flat bed.
    if(held)
    {
        mesh.update();
    }
    while(runTime.run())
    {
        ++runTime;
        displacement.boundaryFieldRef()[patch]==target*(held?scalar(1):runTime.value()/.2);
        const pointField before(mesh.points());
        mesh.update();
        const scalar step=gMax(mag(mesh.points()-before));
        scalar co=0;
        if(!held || runTime.timeIndex()==1 || step>1e-12)
        {
            minimum=gMin(mesh.V().field());
            tmp<volScalarField> sumPhi=fvc::surfaceSum(mag(mesh.phi()));
            co=.5*runTime.deltaTValue()*gMax(sumPhi().primitiveField()/mesh.V().field());
        }
        worstCo=max(worstCo,co);worstStep=max(worstStep,step);
        Info<<"REPLAY time="<<runTime.value()<<" maxPointStep="<<step
            <<" minV="<<minimum<<" meshCo="<<co<<endl;
        if(minimum<=0 || !std::isfinite(co) || co>1 || (held && runTime.timeIndex()>1 && step>1e-12))
        {
            mesh.write();displacement.write();
            Info<<"REPLAY_FAILURE mesh-only instability"<<endl;
            return 1;
        }
        runTime.write();
    }
    Info<<"REPLAY_RESULT maxMeshCo="<<worstCo<<" maxPointStep="<<worstStep<<endl;
    return 0;
}
