"""Validate and simplify the provided binary STL as a five-metre inverted ship."""
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


ROOT = Path(sys.argv[1]).resolve()
GEOM = ROOT/'geometry'
SHA = 'f0209313129d4a36b47ccd97bdd99a1483b40ced701fd3178ba0e85a0c2f0cb4'
LENGTH = 5.0
# The central hull spreads to its full beam at original elevations 12–14.
# Adopt 14.2 as the original deck reference; it becomes z=0 after placement.
DECK_Z_NATIVE = 14.2


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def measure(poly):
    points = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    polys = poly.GetPolys()
    offsets = vtk_to_numpy(polys.GetOffsetsArray())
    if not np.all(np.diff(offsets)==3):
        raise ValueError('Expected triangles throughout the STL')
    faces = vtk_to_numpy(polys.GetConnectivityArray()).reshape(-1,3)
    xyz = points[faces]
    crosses = np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
    if not np.isfinite(xyz).all() or np.any(np.linalg.norm(crosses,axis=1)==0):
        raise ValueError('Nonfinite or degenerate facets')
    bounds = np.array([points.min(0),points.max(0)])
    volume = np.einsum('ij,ij->',xyz[:,0]-bounds.mean(0),crosses)/6
    return points, faces, bounds, volume


def manifold(faces):
    edges = np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    order = np.lexsort((edges.max(1),edges.min(1)))
    sorted_edges = edges[order]
    return (len(sorted_edges)%2==0 and np.array_equal(sorted_edges[::2],sorted_edges[1::2,::-1])
            and not np.any(np.all(np.sort(sorted_edges[1:-1:2],axis=1)==np.sort(sorted_edges[2::2],axis=1),axis=1)))


def build_locator(surface):
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()
    return locator


def sample_error(locator, sample, count=12000):
    points, faces, _, _ = sample
    take = faces[np.linspace(0,len(faces)-1,min(count,len(faces)),dtype=int)]
    coords = np.concatenate((points[take[:,0]],points[take].mean(1)))
    dist = vtk.mutable(0.)
    cell = vtk.mutable(0)
    sub = vtk.mutable(0)
    closest = [0.,0.,0.]
    errors = []
    for xyz in coords:
        locator.FindClosestPoint(xyz,closest,cell,sub,dist)
        errors.append(float(dist)**.5)
    return {'max_m':max(errors),'p95_m':float(np.percentile(errors,95)),'samples':len(errors)}


def save_binary(path, poly):
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToBinary()
    writer.SetInputData(poly)
    if writer.Write()!=1:
        raise RuntimeError('Failed to write candidate STL')


def voxel_surface(pose,bounds,spacing):
    # Rasterize the solid, extract its half-voxel isosurface, then smooth the
    # stair steps before decimation. This separates coincident collapsed edges.
    start = bounds[0]-3*spacing
    dims = np.ceil((bounds[1]-start)/spacing).astype(int)+4
    image = vtk.vtkImageData()
    image.SetDimensions(*map(int,dims))
    image.SetOrigin(*start)
    image.SetSpacing(spacing,spacing,spacing)
    image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)
    vtk_to_numpy(image.GetPointData().GetScalars()).fill(1)
    stencil = vtk.vtkPolyDataToImageStencil()
    stencil.SetInputData(pose)
    stencil.SetOutputOrigin(*start)
    stencil.SetOutputSpacing(spacing,spacing,spacing)
    stencil.SetOutputWholeExtent(image.GetExtent())
    convert = vtk.vtkImageStencil()
    convert.SetInputData(image)
    convert.SetStencilConnection(stencil.GetOutputPort())
    convert.SetBackgroundValue(0)
    convert.Update()
    from scipy.ndimage import gaussian_filter
    binary = vtk_to_numpy(convert.GetOutput().GetPointData().GetScalars()).reshape(tuple(dims[::-1]))
    smooth_grid = gaussian_filter(binary.astype(np.float32),sigma=.82,mode='constant')
    signed_grid = vtk.vtkImageData()
    signed_grid.CopyStructure(image)
    signed_grid.GetPointData().SetScalars(numpy_to_vtk(smooth_grid.ravel(),deep=True))
    surface = vtk.vtkFlyingEdges3D()
    surface.SetInputData(signed_grid)
    surface.SetValue(0,.50037)
    surface.ComputeNormalsOff()
    smooth = vtk.vtkWindowedSincPolyDataFilter()
    smooth.SetInputConnection(surface.GetOutputPort())
    smooth.SetNumberOfIterations(8)
    smooth.SetPassBand(.2)
    smooth.NormalizeCoordinatesOn()
    smooth.BoundarySmoothingOff()
    smooth.NonManifoldSmoothingOff()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(smooth.GetOutputPort())
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.Update()
    result = vtk.vtkPolyData()
    result.ShallowCopy(normals.GetOutput())
    return result


def main():
    original = GEOM/'input.original.stl'
    if not original.exists():
        with TemporaryDirectory(prefix='.hull-original-', dir=GEOM) as staging:
            tmp = Path(staging)/original.name
            with gzip.open(GEOM/'input.stl.gz','rb') as source, tmp.open('wb') as target:
                shutil.copyfileobj(source,target)
            if sha256_file(tmp)!=SHA:
                raise RuntimeError('Uploaded original hash does not match')
            tmp.replace(original)
    elif sha256_file(original)!=SHA:
        raise RuntimeError('Original hash does not match')
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(original))
    reader.MergingOn()
    reader.Update()
    raw = reader.GetOutput()
    _,_,original_bounds,_ = measure(raw)
    scale = LENGTH/(original_bounds[1,0]-original_bounds[0,0])
    centre = original_bounds.mean(0)
    # Proper 180-degree roll about X (determinant +1, winding is preserved).
    transform = vtk.vtkTransform()
    transform.Scale(scale,-scale,-scale)
    transform.Translate(-centre[0],-centre[1],-DECK_Z_NATIVE)
    filt = vtk.vtkTransformPolyDataFilter()
    filt.SetTransform(transform)
    filt.SetInputData(raw)
    filt.Update()
    pose = filt.GetOutput()
    reference = measure(pose)
    if reference[3]<=0 or not manifold(reference[1]):
        raise RuntimeError('Input is not a consistently outward, closed oriented manifold')
    report = {'original_sha256':SHA,'original_faces':int(pose.GetNumberOfPolys()),'length_m':LENGTH,
              'roll_deg':180,'yaw_deg_for_case':90,'deck_reference_native_z':DECK_Z_NATIVE,
              'original_bounds_native':original_bounds.tolist(),'reference_bounds_m':reference[2].tolist(),
              'reference_signed_volume_m3':reference[3],'candidates':[]}
    method = 'direct' if '--direct-only' in sys.argv else 'voxel_filtered'
    if method=='voxel_filtered':
        spacing = float(sys.argv[2]) if len(sys.argv)>2 else .0125
        reduced_input = voxel_surface(pose,reference[2],spacing)
        report['voxel_spacing_m'] = spacing
        report['voxel_isosurface_faces'] = reduced_input.GetNumberOfPolys()
    else:
        reduced_input = pose
    report['method'] = method
    reference_locator = build_locator(pose)
    for target in (30000,60000,120000,240000):
        path = GEOM/f'hull_{method}_{target}.stl'
        decimate = vtk.vtkQuadricDecimation()
        decimate.SetInputData(reduced_input)
        decimate.SetTargetReduction(max(0,1-target/reduced_input.GetNumberOfPolys()))
        decimate.VolumePreservationOn()
        decimate.Update()
        save_binary(path,decimate.GetOutput())
        reread = vtk.vtkSTLReader()
        reread.SetFileName(str(path))
        reread.MergingOn()
        reread.Update()
        candidate = reread.GetOutput()
        metrics = {'target_faces':target,'actual_faces':candidate.GetNumberOfPolys(),'file':str(path)}
        report['candidates'].append(metrics)
        try:
            actual = measure(candidate)
            if not manifold(actual[1]) or actual[3]<=0:
                raise ValueError('Candidate edge incidence or winding failed')
            metrics.update(bounds_m=actual[2].tolist(),volume_m3=actual[3],
                           max_bounds_difference_m=float(abs(actual[2]-reference[2]).max()),
                           volume_relative_difference=abs(actual[3]/reference[3]-1),
                           original_to_candidate=sample_error(build_locator(candidate),reference),
                           candidate_to_original=sample_error(reference_locator,actual))
            if (metrics['max_bounds_difference_m']>.03 or metrics['volume_relative_difference']>.02 or
                max(metrics[n]['max_m'] for n in ['original_to_candidate','candidate_to_original'])>.075):
                raise ValueError('Shape threshold exceeded')
            log = ROOT/'logs'/f'hull_{method}_{target}.surfaceCheck.log'
            with log.open('w') as out:
                result = subprocess.run(['surfaceCheck','-checkSelfIntersection',str(path)],cwd=GEOM,
                                        stdout=out,stderr=subprocess.STDOUT)
            native = log.read_text(errors='replace')
            if result.returncode or 'Surface is closed' not in native or 'Surface is not self-intersecting' not in native:
                raise ValueError('Native closure or self-intersection check failed')
            metrics['passed'] = True
            report.update(status='passed',selected=str(path),sha256=sha256_file(path),
                          actual_faces=actual[1].shape[0],burial_m=float(-actual[2][0,2]*LENGTH/np.ptp(actual[2][:,0])),
                          scale_native_to_m=scale)
            break
        except ValueError as exc:
            metrics.update(passed=False,error=str(exc))
        finally:
            (ROOT/'reports/hull_geometry.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if report.get('status')!='passed':
        raise RuntimeError('No simplified geometry passed. Original and candidate logs preserved.')
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
