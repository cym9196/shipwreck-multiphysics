// Optional inventory diagnostics for a fixed outer domain and one moving bed.
// Keep actual geometry and next-step boundary targets on separate time layers.
const bool writeSedimentBudget =
    bedloadProperties.lookupOrDefault<bool>("writeSedimentBudget", false);
scalar budgetVolume0 = 0;
scalar budgetSuspended0 = 0;
scalar budgetSuspendedBoundaryRate = 0;
scalar budgetSuspendedBoundaryGross = 0;
scalar budgetCumulativeOut = 0;
scalar budgetCumulativeGrossBoundary = 0;
scalar budgetSuspendedPrevious = 0;
scalar budgetBedSuspendedRate = 0;
scalar budgetErosionRate = 0;
scalar budgetFaceStep = 0;
scalar budgetPointStepBeforeFixed = 0;
scalar budgetPointStepAfterFixed = 0;
scalar budgetCumulativeSuspendedEquation = 0;
scalar budgetCumulativeExchange = 0;
scalar budgetCumulativeExner = 0;
scalar budgetCumulativeInterpolation = 0;
scalar budgetCumulativeFixed = 0;
scalar budgetCumulativePointSteps = 0;
autoPtr<vectorField> budgetBedPoints0;
autoPtr<vectorField> budgetBedDisplacement0;
autoPtr<OFstream> budgetStream;
if (writeSedimentBudget)
{
    if (!bed.exist() || !bed.bedMotion() || bed.bedPatchesID().size() != 1
        || morphoAccFactor != 1 || ABorder != 0 || switchMULES != "off")
    {
        FatalErrorInFunction
            << "Inventory diagnostics require one moving bed, morphoAccFactor=1, "
            << "ABorder=0 and MULESformulation=off" << abort(FatalError);
    }
    budgetVolume0 = gSum(mesh.V().field());
    if (switchSuspension == "on")
    {
        budgetSuspended0 = gSum(CsPtr.ref().primitiveField()*mesh.V().field());
    }
    budgetSuspendedPrevious = budgetSuspended0;
    budgetBedPoints0.reset(new vectorField(bed.aMesh().patch().localPoints()));
    const pointVectorField& pd = mesh.lookupObject<pointVectorField>("pointDisplacement");
    budgetBedDisplacement0.reset
    (
        new vectorField
        (
            refCast<const fixedValuePointPatchVectorField>
            (
                pd.boundaryField()[bed.bedPatchesID()[0]]
            )
        )
    );

    IOobject budgetRestartIO("sedimentBudgetState",runTime.timeName(),"uniform",mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE);
    if (budgetRestartIO.typeHeaderOk<localIOdictionary>(true))
    {
        localIOdictionary restored(budgetRestartIO);
        if (mag(restored.get<scalar>("savedTime")-runTime.value())>1e-10)
            FatalErrorInFunction<<"Budget checkpoint clock mismatch"<<exit(FatalError);
        budgetVolume0 = restored.get<scalar>("volume0");
        budgetSuspended0 = restored.get<scalar>("suspended0");
        budgetSuspendedPrevious = restored.get<scalar>("suspendedPrevious");
        budgetCumulativeOut = restored.get<scalar>("cumulativeOut");
        budgetCumulativeGrossBoundary = restored.get<scalar>("cumulativeGrossBoundary");
        budgetCumulativeSuspendedEquation = restored.get<scalar>("cumulativeSuspendedEquation");
        budgetCumulativeExchange = restored.get<scalar>("cumulativeExchange");
        budgetCumulativeExner = restored.get<scalar>("cumulativeExner");
        budgetCumulativeInterpolation = restored.get<scalar>("cumulativeInterpolation");
        budgetCumulativeFixed = restored.get<scalar>("cumulativeFixed");
        budgetCumulativePointSteps = restored.get<scalar>("cumulativePointSteps");
        budgetBedPoints0.ref() = vectorField(restored.lookup("initialBedPoints"));
        budgetBedDisplacement0.ref() = vectorField(restored.lookup("initialBedDisplacement"));
        if (budgetBedPoints0().size()!=bed.aMesh().nPoints() || budgetBedDisplacement0().size()!=bed.aMesh().nPoints())
            FatalErrorInFunction<<"Budget checkpoint partition mismatch"<<exit(FatalError);
        Info<<"Restored complete sediment inventory at "<<runTime.timeName()<<endl;
    }
    else if (runTime.value()>SMALL)
        FatalErrorInFunction<<"Nonzero restart requires sedimentBudgetState"<<exit(FatalError);

    if (Pstream::master())
    {
        const fileName directory = runTime.globalPath()/"postProcessing"/"sedimentBudget";
        mkDir(directory);
        budgetStream.reset(new OFstream(directory/"budget.csv"));
        budgetStream().precision(16);
        budgetStream() << "time_s,suspended_solid_m3,actual_bed_change_solid_m3,"
            << "pending_point_target_solid_m3,cumulative_external_out_solid_m3,"
            << "cumulative_gross_boundary_solid_m3,actual_raw_residual_solid_m3,"
            << "target_raw_residual_solid_m3,"
            << "suspended_equation_cumulative_m3,exchange_cumulative_m3,"
            << "exner_cumulative_m3,interpolation_cumulative_m3,"
            << "fixed_boundary_cumulative_m3,geometry_quadrature_m3,"
            << "decomposition_closure_m3" << endl;
    }
}
