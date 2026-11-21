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

#include "sedimentBed.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedimentBed::sedimentBed
(
    dictionary& dict,
    Foam::fvMesh& mesh,
    const meshObjects::gravity& g
)
:
    dict_(dict),
    bedExist_(false),
    bedMotion_(false),
    rigidBed_(false),
    mesh_(mesh),
    bedloadModel_(nullptr),
    g_(g),
    eg_((g_ / Foam::mag(g)).value()),
    betaPtr_(nullptr),
    slopeDirPtr_(nullptr)
{
    checkBedExistence();
    if (bedExist_)
    {
        Info << "create finite-area mesh" << endl;
        aMesh_.reset(new faMesh(mesh_));
        vsm.reset(new volSurfaceMapping(aMesh_));
        getPatchesID();
        // check if bed inclination consistent with gravity
        checkFaMeshOrientation();
        // instantiate projected finite-area mesh
        aProjMesh_.reset(new projectedFaMesh(aMesh_.ref(), eg_));
        // instantiate bedloadModel
        getBedloadModel();
        // instantiate criticalShieldsModel
        getCritShieldsModel();

        // activate/deactivate presence of a rigid bed below sediment layer
        word switchRigidBed =
            dict_.lookupOrDefault<word>("rigidBed", "off");
        if (switchRigidBed=="on")
        {
            rigidBed_ = true;
        }
        
        Info << "faMesh patches ID : " << bedPatchesID_ << endl;
        Info << "faMesh patches names : " << bedPatchesNames_ << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedimentBed::~sedimentBed()
{
    deleteDemandDrivenData(betaPtr_);
    deleteDemandDrivenData(slopeDirPtr_);
}

// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

const faMesh& Foam::sedimentBed::aMesh()
{
    return aMesh_.ref();
}

const Foam::projectedFaMesh& Foam::sedimentBed::aProjMesh()
{
    return aProjMesh_.ref();
}

const Foam::bedloadModels::bedloadModel&
Foam::sedimentBed::bedloadModel()
{
    return bedloadModel_.ref();
}

Foam::criticalShieldsModels::criticalShieldsModel&
Foam::sedimentBed::critShieldsModel()
{
    return critShieldsModel_.ref();
}

labelList Foam::sedimentBed::bedPatchesID()
{
    return bedPatchesID_.clone();
}

bool Foam::sedimentBed::exist() const
{
    if (bedExist_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::sedimentBed::bedMotion() const
{
    if (bedMotion_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::sedimentBed::rigidBed() const
{
    if (rigidBed_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

const Foam::scalarField& Foam::sedimentBed::beta() const
{
    if (betaPtr_==nullptr)
    {
        calcBeta();
    }

    return *betaPtr_;
}

const Foam::vectorField& Foam::sedimentBed::slopeDir() const
{
    if (slopeDirPtr_==nullptr)
    {
        calcSlopeDir();
    }

    return *slopeDirPtr_;
}

void Foam::sedimentBed::interpFaceToVertices
(
    scalarField& dHfaces,
    scalarField& dHpoints,
    const labelList& fixedPoints
)
{
    // Access to face centers in projected horizontal plane
    const vectorField& xFacesProj = aProjMesh_->areaCentres();

    // Access to points coordinates in projected horizontal plane
    const vectorField& xPointsProj = aProjMesh_->pointCoords();

    // normalise interpolation weight
    scalarField normWeight(aMesh_->nPoints());
    normWeight = Zero;
    const bool conservativeFixed =
        dict_.lookupOrDefault<bool>("conservativeFixedPointProjection", false);
    const labelList& meshPoints = aMesh_->patch().meshPoints();
    scalarField fixed(mesh_.nPoints(), Zero);
    if (conservativeFixed)
    {
        for (const label pointi : fixedPoints) fixed[meshPoints[pointi]] = 1;
        syncTools::syncPointList(mesh_, fixed, maxEqOp<scalar>(), scalar(0));
    }

    // loop over all faces
    for (label facei=0; facei < aMesh_->nFaces(); facei++)
    {
        const face& f = aMesh_->faces()[facei];
        const vector& xf = xFacesProj[facei];
        // number of vertices defining face f
        const label nVerts = f.size();
        scalarField weights(nVerts, Zero);

        // loop over points belonging to face
        for (label pi=0; pi < nVerts; pi++)
        {
            // label of point pi, global
            const label piGlob = f.thisLabel(pi);

            scalar areaPF = 0;

            // vertices coordinates
            const vector& xv = xPointsProj[piGlob];
            const vector& xn = xPointsProj[f.nextLabel(pi)];
            const vector& xp = xPointsProj[f.prevLabel(pi)];

            // vector from vertice to face center
            vector vp = xf - xv;
            // vector from vertice to right edge center
            vector vredg = 0.5 * (xn - xv);
            // vector from vertice to left edge center
            vector vledg = 0.5 * (xp - xv);

            areaPF +=
                0.5 * Foam::mag(vredg ^ vp)
                + 0.5 * Foam::mag(vp ^ vledg);

            weights[pi] = areaPF;
            normWeight[piGlob] += areaPF;
        }
        scalar freeWeight = 0;
        const scalar allWeight = sum(weights);
        forAll(f, pi)
        {
            if (!conservativeFixed || fixed[meshPoints[f[pi]]] == 0)
                freeWeight += weights[pi];
        }
        if (freeWeight <= VSMALL && mag(dHfaces[facei]) > SMALL)
        {
            FatalErrorInFunction << "Sediment exchange on a fully constrained face"
                << abort(FatalError);
        }
        forAll(f, pi)
        {
            if (!conservativeFixed || fixed[meshPoints[f[pi]]] == 0)
            {
                // Distribute each face's solid volume only to free vertices;
                // retain the full geometric quadrature weight in the denominator.
                dHpoints[f[pi]] += weights[pi]*allWeight/max(freeWeight,VSMALL)*dHfaces[facei];
            }
        }
    }
    // Sum the numerator and weight before normalising shared bed vertices.
    // Averaging already-normalised local displacements loses contributions
    // when a processor boundary cuts faces of different projected areas.
    scalarField numerator(mesh_.nPoints(), Zero);
    scalarField denominator(mesh_.nPoints(), Zero);
    forAll(meshPoints, pointi)
    {
        numerator[meshPoints[pointi]] = dHpoints[pointi];
        denominator[meshPoints[pointi]] = normWeight[pointi];
    }
    syncTools::syncPointList(mesh_, numerator, plusEqOp<scalar>(), scalar(0));
    syncTools::syncPointList(mesh_, denominator, plusEqOp<scalar>(), scalar(0));
    forAll(meshPoints, pointi)
    {
        const label meshPointi = meshPoints[pointi];
        if (denominator[meshPointi] <= VSMALL)
        {
            FatalErrorInFunction << "Zero bed interpolation weight" << abort(FatalError);
        }
        dHpoints[pointi] = numerator[meshPointi]/denominator[meshPointi];
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sedimentBed::checkBedExistence()
{
    word isBed(dict_.lookup("sedimentBed"));
    if (isBed=="on")
    {
        bedExist_ = true;
        checkBedMotion();
        //checkAvalancheModel_();
    }
    else if (isBed=="off")
    {
        bedExist_ = false;
    }
    else
    {
        FatalError << "wrong keyword for entry sedimentBed!"
            << " possible options are: on off" << endl;
        Info << abort(FatalError) << endl;
    }
}

void Foam::sedimentBed::checkFaMeshOrientation() const
{
    if (not bedExist_)
    {
        return;
    }
    Info << "g : " << g_.value() << endl;
    forAll(aMesh_->faceLabels(), facei)
    {
        double res = g_.value() & aMesh_->faceAreaNormals()[facei];
        if (res <= 0)
        {
            FatalError
                << "sediment bed and gravity orientation "
                << "are not consistent" << endl
                << "g and patch faces dot product "
                << "should be positive" << endl;
            Info << abort(FatalError) << endl;
        }
    }
}

void Foam::sedimentBed::getPatchesID()
{
    List<word> bedPatchesNames(dict_.lookup("sedimentBedPatches"));
    forAll(bedPatchesNames, i)
    {
        word patchName = bedPatchesNames[i];
        label patchID = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchID==-1)
        {
            FatalError << "bedPatch " << patchName
                << " does not exist" << endl
                << "existing patches are:" << endl
                << mesh_.boundaryMesh().names() << endl;
            Info << abort(FatalError) << endl;
        }
        bedPatchesNames_.append(patchName);
        bedPatchesID_.append(patchID);
    }
}

void Foam::sedimentBed::checkBedMotion()
{
    if (dict_.found("bedMotion"))
    {
        word bedMotionState(dict_.lookup("bedMotion"));
        if (bedMotionState=="on")
        {
            bedMotion_ = true;
            Info << "bed motion activated" << endl;
        }
        else if (bedMotionState=="off")
        {
            bedMotion_ = false;
            Info << "no bed motion, "
                << "bedload and shields are computed "
                << "but bed will not move" << endl;
        }
        else
        {
            FatalError << "wrong keyword for entry bedMotion!"
                << " possible options are: on off,"
                << " default is on" << endl;
            Info << abort(FatalError) << endl;
        }
    }
    else
    {
        bedMotion_ = true;
    }
}

void Foam::sedimentBed::getBedloadModel()
{
    // default dict in case bedloadModel not specified
    dictionary defaultBedloadDict;
    defaultBedloadDict.add(keyType("type"), "MeyerPeter");

    if (dict_.found("bedloadModel"))
    {
        const dictionary bedloadModelDict =
            dict_.subDict("bedloadModel");

        bedloadModel_.reset
        (
            bedloadModels::bedloadModel::New
            (
                bedloadModelDict
            )
        );
    }
    else
    {
        bedloadModel_.reset(
            bedloadModels::bedloadModel::New
            (
                defaultBedloadDict
            )
        );
    }
}

void Foam::sedimentBed::getCritShieldsModel()
{
    // default dict in case criticalShieldsModel not specified
    dictionary defaultCritShieldsDict;
    defaultCritShieldsDict.add(keyType("type"), "fixedValue");
    defaultCritShieldsDict.add(keyType("value"), 0.047);
    if (dict_.found("criticalShieldsModel"))
    {
        const dictionary& critShieldsDict =
            dict_.subDict("criticalShieldsModel");
        word modelType(critShieldsDict.lookup("type"));
        if (modelType != "readFromFile")
        {
            critShieldsModel_.reset
                (
                    criticalShieldsModels::criticalShieldsModel::New
                    (
                        critShieldsDict
                    )
                );
        }
    }
    else
    {
        critShieldsModel_.reset
        (
            criticalShieldsModels::criticalShieldsModel::New
            (
                defaultCritShieldsDict
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
