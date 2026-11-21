#ifndef ShipwreckHullSurface_H
#define ShipwreckHullSurface_H

#include "triSurfaceMesh.H"

namespace Foam {

class HullSurface {
    const triSurfaceMesh surface_;
    const scalar thickness_;
    boundBox searchBox_;

public:
    HullSurface(const fvMesh& mesh, const word& geometry, scalar thickness)
    : surface_(IOobject(geometry, mesh.time().constant(), "triSurface", mesh,
                        IOobject::MUST_READ, IOobject::NO_WRITE), triSurfaceMesh::localOrGlobal),
      thickness_(thickness), searchBox_(surface_.localPoints(), false)
    {
        if (!(thickness_ > 0) || !surface_.hasVolumeType()) {
            FatalErrorInFunction << "A closed hull and positive interface thickness are required" << exit(FatalError);
        }
        const vector padding = vector::one * (8 * thickness_);
        searchBox_.min() -= padding;
        searchBox_.max() += padding;
    }

    scalarField signedDistances(const pointField& positions) const {
        DynamicList<label> indices;
        DynamicList<point> samples;
        forAll(positions, i) if (searchBox_.contains(positions[i])) {
            indices.append(i);
            samples.append(positions[i]);
        }
        scalarField distances(positions.size(), GREAT);
        if (samples.empty()) return distances;

        const pointField query(std::move(samples));
        List<pointIndexHit> hits;
        List<volumeType> types;
        surface_.findNearest(query, scalarField(query.size(), sqr(searchBox_.mag())), hits);
        surface_.getVolumeType(query, types);
        forAll(query, i) {
            if (!hits[i].hit()) FatalErrorInFunction << "Hull surface distance failed" << exit(FatalError);
            const label cell = indices[i];
            const scalar distance = mag(query[i] - hits[i].hitPoint());
            if (distance <= SMALL) distances[cell] = 0;
            else if (types[i] == volumeType::INSIDE) distances[cell] = -distance;
            else if (types[i] == volumeType::OUTSIDE) distances[cell] = distance;
            else FatalErrorInFunction << "Ambiguous hull surface classification at " << query[i] << exit(FatalError);
        }
        return distances;
    }
};

} // namespace Foam
#endif
