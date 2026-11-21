if (writeSedimentBudget)
{
    const faMesh& bm = bed.aMesh();
    const scalarField& projectedArea = bed.aProjMesh().S();
    const vectorField& pointsNow = bm.patch().localPoints();
    const pointVectorField& pd = mesh.lookupObject<pointVectorField>("pointDisplacement");
    const fixedValuePointPatchVectorField& displacement =
        refCast<const fixedValuePointPatchVectorField>
        (
            pd.boundaryField()[bed.bedPatchesID()[0]]
        );
    scalar pending = 0;
    forAll(bm.faces(), facei)
    {
        const face& f = bm.faces()[facei];
        scalar meanPending = 0;
        forAll(f, pi)
        {
            const label pointi = f[pi];
            const vector remaining = budgetBedPoints0()[pointi]
                + displacement[pointi] - budgetBedDisplacement0()[pointi]
                - pointsNow[pointi];
            meanPending += remaining & (-eg);
        }
        // Rectangular projected faces with vertical-only mesh motion.
        pending += projectedArea[facei]*meanPending/f.size();
    }
    pending = CsMax.value()*returnReduce(pending, sumOp<scalar>());
    scalar suspended = 0;
    if (switchSuspension == "on")
    {
        suspended = gSum(CsPtr.ref().primitiveField()*mesh.V().field());
    }
    const scalar actualBed = CsMax.value()*(budgetVolume0-gSum(mesh.V().field()));
    scalar bedOut = 0;
    scalar bedGross = 0;
    forAll(bm.boundary(), patchi)
    {
        if (!qbPtr.ref().boundaryField()[patchi].coupled()
            && bm.boundary()[patchi].type() != "empty")
        {
            const faPatch& patch = bm.boundary()[patchi];
            forAll(patch, ei)
            {
                const scalar flux = phiqbProjPtr.ref()[patch.start()+ei];
                bedOut += flux;
                bedGross += mag(flux);
            }
        }
    }
    bedOut = returnReduce(bedOut, sumOp<scalar>());
    bedGross = returnReduce(bedGross, sumOp<scalar>());
    budgetCumulativeOut += runTime.deltaTValue()*(bedOut+budgetSuspendedBoundaryRate);
    budgetCumulativeGrossBoundary += runTime.deltaTValue()
        *(bedGross+budgetSuspendedBoundaryGross);
    const scalar actualResidual = suspended-budgetSuspended0
        + actualBed+budgetCumulativeOut;
    const scalar dt = runTime.deltaTValue();
    scalar exchangeRate = 0;
    if (switchSuspension == "on")
    {
        exchangeRate = gSum(DEptr.ref().primitiveField()*bm.S().field());
    }
    budgetCumulativeSuspendedEquation += suspended-budgetSuspendedPrevious
        + dt*(budgetSuspendedBoundaryRate+budgetBedSuspendedRate-budgetErosionRate);
    budgetSuspendedPrevious = suspended;
    budgetCumulativeExchange += dt*(exchangeRate+budgetErosionRate-budgetBedSuspendedRate);
    budgetCumulativeExner += budgetFaceStep+dt*(bedOut-exchangeRate);
    budgetCumulativeInterpolation += budgetPointStepBeforeFixed-budgetFaceStep;
    budgetCumulativeFixed += budgetPointStepAfterFixed-budgetPointStepBeforeFixed;
    budgetCumulativePointSteps += budgetPointStepAfterFixed;
    const scalar geometry = actualBed+pending-budgetCumulativePointSteps;
    const scalar explained = budgetCumulativeSuspendedEquation+budgetCumulativeExchange
        + budgetCumulativeExner+budgetCumulativeInterpolation+budgetCumulativeFixed+geometry;
    if (Pstream::master())
    {
        budgetStream() << runTime.value() << ',' << suspended << ',' << actualBed
            << ',' << pending << ',' << budgetCumulativeOut << ','
            << budgetCumulativeGrossBoundary << ',' << actualResidual << ','
            << actualResidual+pending << ',' << budgetCumulativeSuspendedEquation
            << ',' << budgetCumulativeExchange << ',' << budgetCumulativeExner
            << ',' << budgetCumulativeInterpolation << ',' << budgetCumulativeFixed
            << ',' << geometry << ',' << actualResidual+pending-explained << endl;
    }

    if (runTime.writeTime())
    {
        localIOdictionary checkpoint(IOobject("sedimentBudgetState",runTime.timeName(),"uniform",mesh,IOobject::NO_READ,IOobject::NO_WRITE));
        checkpoint.add("savedTime",runTime.value());
        checkpoint.add("volume0",budgetVolume0);
        checkpoint.add("suspended0",budgetSuspended0);
        checkpoint.add("suspendedPrevious",budgetSuspendedPrevious);
        checkpoint.add("cumulativeOut",budgetCumulativeOut);
        checkpoint.add("cumulativeGrossBoundary",budgetCumulativeGrossBoundary);
        checkpoint.add("cumulativeSuspendedEquation",budgetCumulativeSuspendedEquation);
        checkpoint.add("cumulativeExchange",budgetCumulativeExchange);
        checkpoint.add("cumulativeExner",budgetCumulativeExner);
        checkpoint.add("cumulativeInterpolation",budgetCumulativeInterpolation);
        checkpoint.add("cumulativeFixed",budgetCumulativeFixed);
        checkpoint.add("cumulativePointSteps",budgetCumulativePointSteps);
        checkpoint.add("initialBedPoints",budgetBedPoints0());
        checkpoint.add("initialBedDisplacement",budgetBedDisplacement0());
        if (!checkpoint.writeObject(IOstreamOption(IOstreamOption::ASCII),true))
            FatalErrorInFunction<<"Could not persist sediment inventory"<<exit(FatalError);
    }
}
