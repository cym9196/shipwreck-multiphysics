"""Collapse one documented edge in an independent binary STL copy."""
import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--face', type=int, required=True)
    parser.add_argument('--from-vertex', type=int, choices=range(3), required=True)
    parser.add_argument('--to-vertex', type=int, choices=range(3), required=True)
    parser.add_argument('--expected-sha256', required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise ValueError('Output must be a new file')
    original_hash = file_hash(args.source)
    if original_hash != args.expected_sha256:
        raise ValueError('Selected face index is only valid for the specified input hash')
    dtype = np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attribute','<u2')])
    with args.source.open('rb') as stream:
        header = stream.read(84)
        count = struct.unpack_from('<I',header,80)[0]
        if args.source.stat().st_size != 84+50*count or not 0<=args.face<count:
            raise ValueError('Unexpected triangle count or binary STL size')
        faces = np.fromfile(stream,dtype=dtype,count=count)
    vertices = faces['vertices']
    old = vertices[args.face,args.from_vertex].copy()
    new = vertices[args.face,args.to_vertex].copy()
    if np.array_equal(old,new):
        raise ValueError('Choose distinct edge endpoints')
    using_old = (vertices==old).all(axis=2)
    using_new = (vertices==new).all(axis=2)
    old_star = using_old.any(axis=1)
    new_star = using_new.any(axis=1)
    removed = old_star & new_star
    # The link condition preserves the local manifold when this edge is collapsed.
    old_neighbors = {tuple(v) for v in vertices[old_star].reshape(-1,3)}-{tuple(old),tuple(new)}
    new_neighbors = {tuple(v) for v in vertices[new_star].reshape(-1,3)}-{tuple(old),tuple(new)}
    opposite = {tuple(v) for v in vertices[removed].reshape(-1,3)}-{tuple(old),tuple(new)}
    if int(removed.sum()) != 2 or len(opposite)!=2 or old_neighbors & new_neighbors != opposite:
        raise ValueError('One interior manifold edge with a valid link is required')
    extent = np.ptp(vertices,axis=(0,1))
    incident = old_star & ~removed
    earlier = vertices[incident].astype(float)
    vertices[using_old] = new
    later = vertices[incident].astype(float)
    normal = lambda v: np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])
    new_normals = normal(later)
    magnitudes = np.linalg.norm(new_normals,axis=1)
    if (magnitudes<=0).any() or (np.sum(normal(earlier)*new_normals,axis=1)<=0).any():
        raise ValueError('Collapse would degenerate or invert an adjacent face')
    faces['normal'][incident] = new_normals/magnitudes[:,None]
    surviving = faces[~removed]
    with args.output.open('xb') as stream:
        stream.write(header[:80]+struct.pack('<I',len(surviving)))
        surviving.tofile(stream)
    summary = {'source_sha256':original_hash,'candidate_sha256':file_hash(args.output),
               'selected_source_face_zero_based':args.face,'from_vertex':old.tolist(),'to_vertex':new.tolist(),
               'removed_triangles':np.flatnonzero(removed).tolist(),'modified_triangles':np.flatnonzero(incident).tolist(),
               'endpoint_move_fraction_of_hull_length':float(np.linalg.norm(new-old)/extent[0]),
               'native_surface_validation_required':True}
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
