/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

vectorField& dispVals = dispValsPtr.ref();
if (bed.bedPatchesID().size() != 1)
{
    FatalErrorInFunction
        << "This displacement mapping currently requires one bed polyPatch"
        << abort(FatalError);
}
for (const label patchID : bed.bedPatchesID())
{
    // Accumulate from the prescribed boundary target, not from an
    // interpolated volume-internal displacement near the bed.
    const fixedValuePointPatchVectorField& bedDisplacement =
        refCast<const fixedValuePointPatchVectorField>
        (
            pointDisplacement.boundaryField()[patchID]
        );
    if (bedDisplacement.size() != dispVals.size())
    {
        FatalErrorInFunction << "Bed point mapping size mismatch" << abort(FatalError);
    }
    dispVals = bedDisplacement;
}


scalarField& dHpoints = dHpointsPtr.ref();
dHpoints = Zero;

//interpolate bed level increment from face center to vertices
const labelList& fixedPointsLabels = fixedPatchPointLabelsPtr.ref();
bed.interpFaceToVertices(dH, dHpoints, fixedPointsLabels);

auto integratedBedPointStep = [&]()
{
    scalar volume = 0;
    const scalarField& area = bed.aProjMesh().S();
    forAll(bed.aMesh().faces(), facei)
    {
        const face& f = bed.aMesh().faces()[facei];
        scalar height = 0;
        forAll(f, pi) height += dHpoints[f[pi]];
        volume += area[facei]*height/f.size();
    }
    return CsMax.value()*returnReduce(volume, sumOp<scalar>());
};
if (writeSedimentBudget) budgetPointStepBeforeFixed = integratedBedPointStep();

// zero displacements for faPatches

for (const label i : fixedPointsLabels)
{
    dHpoints[i] = 0.;
}
if (writeSedimentBudget) budgetPointStepAfterFixed = integratedBedPointStep();

dispVals -= eg * dHpoints;

// interpolate displacement on points, vertical orientation
//dispVals -= patchInterp.faceToPointInterpolate(dH) * eg;

for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << max(mag(dispVals));
    Info << " m , min:" << min(mag(dispVals));
    Info << " m , mean:" << average(mag(dispVals)) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID]
        == tmp<vectorField>::New(dispVals);
}
pointDisplacement.correctBoundaryConditions();
