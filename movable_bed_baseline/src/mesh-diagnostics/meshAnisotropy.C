// Reproduce every native determinant flag from current face geometry.
#include "fvCFD.H"
#include "primitiveMeshTools.H"
#include "syncTools.H"
#include "timeSelector.H"
#include <fstream>
#include <iomanip>
#include <cmath>

int main(int argc,char *argv[])
{
    Foam::timeSelector::addOptions();
    Foam::argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList times=Foam::timeSelector::select0(runTime,args);
    forAll(times,ti)
    {
        runTime.setTime(times[ti],ti);
        Foam::fvMesh mesh(Foam::IOobject(Foam::fvMesh::defaultRegion,runTime.timeName(),runTime,Foam::IOobject::MUST_READ));
        forAll(mesh.boundaryMesh(),patchi)
            if(mesh.boundaryMesh()[patchi].coupled())
                FatalErrorInFunction<<"Reconstruct the mesh before attribution"<<Foam::exit(Foam::FatalError);
        const auto& sf=mesh.faceAreas();
        const auto& owner=mesh.faceOwner();const auto& neighbour=mesh.faceNeighbour();
        Foam::List<Foam::symmTensor> moment(mesh.nCells(),Foam::symmTensor::zero);
        Foam::scalarField areaSum(mesh.nCells(),0);Foam::labelList count(mesh.nCells(),0);
        forAll(neighbour,facei)
        {
            const Foam::symmTensor outer=Foam::sqr(sf[facei]);
            const Foam::scalar area=Foam::mag(sf[facei]);
            const Foam::label a=owner[facei],b=neighbour[facei];
            moment[a]+=outer;moment[b]+=outer;areaSum[a]+=area;areaSum[b]+=area;++count[a];++count[b];
        }
        Foam::tmp<Foam::scalarField> native=Foam::primitiveMeshTools::cellDeterminant(
            mesh,mesh.geometricD(),sf,Foam::syncTools::getInternalOrCoupledFaces(mesh));
        const Foam::fileName outDir=runTime.path()/"quality_attribution"/runTime.timeName();Foam::mkDir(outDir);
        std::ofstream cells((outDir/"flagged_cells.f64").c_str(),std::ios::binary);
        Foam::label flagged=0,unexplained=0;double maxRelative=0,minNative=Foam::GREAT;
        forAll(count,celli)
        {
            if(count[celli]==0)FatalErrorInFunction<<"No active internal face"<<Foam::exit(Foam::FatalError);
            const double average=areaSum[celli]/count[celli];
            const double value=Foam::det(moment[celli])/(8*std::pow(average,6));
            const double reference=native()[celli];
            const double error=std::abs(value-reference)/Foam::max(std::abs(reference),1e-100);
            maxRelative=Foam::max(maxRelative,error);minNative=Foam::min(minNative,reference);
            if(reference<.001)
            {
                ++flagged;
                const Foam::vector eig=Foam::eigenValues(moment[celli]);
                const double lo=Foam::min(eig.x(),Foam::min(eig.y(),eig.z()));
                const double hi=Foam::max(eig.x(),Foam::max(eig.y(),eig.z()));
                if(!std::isfinite(value)||error>1e-8||lo<=0)++unexplained;
                const double row[7]={double(celli),reference,value,lo,hi,double(count[celli]),average};
                cells.write(reinterpret_cast<const char*>(row),sizeof(row));
            }
        }
        cells.close();if(!cells)FatalErrorInFunction<<"Attribution output failed"<<Foam::exit(Foam::FatalError);
        std::ofstream report((outDir/"summary.json").c_str());
        report<<std::setprecision(17)<<"{\"cells\":"<<mesh.nCells()<<",\"flagged_cells\":"<<flagged
              <<",\"unexplained_flags\":"<<unexplained<<",\"max_relative_formula_difference\":"<<maxRelative
              <<",\"minimum_native_determinant\":"<<minNative<<",\"threshold_unchanged\":0.001,"
              <<"\"record_format\":\"little-endian float64 rows\",\"record_columns\":[\"cell_id\",\"native_determinant\",\"reproduced_determinant\",\"face_tensor_min_eigenvalue\",\"face_tensor_max_eigenvalue\",\"internal_face_count\",\"mean_internal_face_area_m2\"],"
              <<"\"mechanism\":\"Current face-area tensor anisotropy and internal-face connectivity; reproduction explains flags but does not establish numerical harmlessness.\",\"formal_validation_passed\":false}\n";
        if(unexplained)FatalErrorInFunction<<"Unexplained determinant flags"<<Foam::exit(Foam::FatalError);
        Foam::Info<<"ATTRIBUTION_COMPLETE "<<runTime.timeName()<<" flagged="<<flagged<<Foam::endl;
    }
    return 0;
}
