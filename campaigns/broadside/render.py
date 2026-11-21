"""Render actual OpenFOAM samples with fixed scales and export a ParaView timeline."""
from pathlib import Path
import csv
import json
import shutil
import sys
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
import matplotlib.tri as tri
import numpy as np
from PIL import GifImagePlugin, Image, ImageChops
from scipy.interpolate import LinearNDInterpolator
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROI = (-4.5,8.,-4.5,4.5)
SURFACES = ('flow', 'nearBed', 'bed')


def vtk_sample(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    converter = vtk.vtkTriangleFilter()
    converter.SetInputConnection(reader.GetOutputPort())
    converter.PassVertsOff()
    converter.PassLinesOff()
    converter.Update()
    poly = converter.GetOutput()
    if not poly.GetNumberOfPolys(): raise ValueError(f'No surface triangles at {path}')
    pts = vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
    faces = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1,3).copy()
    data = poly.GetCellData()
    values = {}
    for name in ('U','Cs','shipwreckMask'):
        arr = data.GetArray(name)
        if arr is None: raise ValueError(f'Missing cell field {name} in {path}')
        a = vtk_to_numpy(arr).astype(float)
        if len(a)!=len(faces) or not np.isfinite(a).all(): raise ValueError(f'Invalid sampled {name} in {path}')
        values[name] = a
    return pts,faces,values


def section(hull, elevation):
    plane = vtk.vtkPlane()
    plane.SetOrigin(0,0,elevation)
    plane.SetNormal(0,0,1)
    cutter = vtk.vtkCutter()
    cutter.SetInputData(hull)
    cutter.SetCutFunction(plane)
    cutter.Update()
    output = cutter.GetOutput()
    if output.GetNumberOfLines()==0: return []
    pts = vtk_to_numpy(output.GetPoints().GetData())
    edges = output.GetLines()
    offsets = vtk_to_numpy(edges.GetOffsetsArray())
    conn = vtk_to_numpy(edges.GetConnectivityArray())
    return [pts[conn[a:b],:2].copy() for a,b in zip(offsets,offsets[1:])]


def xy_area(xyz):
    u,v = xyz[...,1,:2]-xyz[...,0,:2],xyz[...,2,:2]-xyz[...,0,:2]
    return .5*abs(u[...,0]*v[...,1]-u[...,1]*v[...,0])


def integral_positive(zones):
    """Integrate max(z,0) over the bed using linear triangle heights."""
    z = zones[:,:,2]
    area = xy_area(zones)
    all_positive = (z>=0).all(1)
    total = float(np.sum(area[all_positive]*z[all_positive].mean(1)))
    for vertices in zones[(z.min(1)<0)&(z.max(1)>0)]:
        poly = []
        for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
            if a[2]>=0: poly.append(a)
            if (a[2]<0)!=(b[2]<0):
                intersection = a + (b-a)*(-a[2]/(b[2]-a[2]))
                intersection[2] = 0
                poly.append(intersection)
        for i in range(1,len(poly)-1):
            xyz = np.asarray((poly[0],poly[i],poly[i+1]))
            total += float(xy_area(xyz)*xyz[:,2].mean())
    return total


def bed_frame(sample):
    pts,faces,data = sample
    xyz = pts[faces]
    outside = data['shipwreckMask']<.5
    if not outside.any(): raise ValueError('No exposed bed faces')
    valid = xyz[outside]
    scour = valid.copy();scour[:,:,2] *= -1
    return {'max_scour_mm':max(0.,-float(valid[:,:,2].min())*1000),
            'max_deposition_mm':max(0.,float(valid[:,:,2].max())*1000),
            'scour_bulk_m3':integral_positive(scour),'deposition_bulk_m3':integral_positive(valid)}


def draw_surface(ax,sample,values,mask,cmap,norm):
    pts,faces,_ = sample
    raster = tri.Triangulation(pts[:,0],pts[:,1],faces)
    raster.set_mask(mask | (xy_area(pts[faces])<=1e-14))
    return ax.tripcolor(raster,facecolors=values,cmap=cmap,norm=norm,shading='flat',rasterized=True)


def flow_lines(ax,sample,x,y):
    pts,faces,data = sample
    centres = pts[faces].mean(1)
    idx = (centres[:,0]>=ROI[0]-.5)&(centres[:,0]<=ROI[1]+.5)&(abs(centres[:,1])<=ROI[3]+.5)
    if np.count_nonzero(idx)<4: return
    fields = np.column_stack((data['U'][:,:2],data['shipwreckMask']))
    values = LinearNDInterpolator(centres[idx,:2],fields[idx])(x,y)
    mask = (values[:,:,2]>=.5) | ~np.isfinite(values[:,:,2])
    ax.streamplot(x[0,:],y[:,0],np.ma.masked_where(mask,values[:,:,0]),np.ma.masked_where(mask,values[:,:,1]),
                  color='#f7f5ed',linewidth=.48,density=.88,arrowsize=.6)


def sample_loads(path, times, field):
    with path.open(newline='') as stream:
        values = np.asarray([(float(row['time_s']), float(row[field]))
                             for row in csv.DictReader(stream)], dtype=float)
    if not values.size or not np.isfinite(values).all():
        raise ValueError(f'Empty or nonfinite load data in {path}')
    log_times, first_rows = np.unique(values[:, 0], return_index=True)
    insertion = np.searchsorted(log_times, times)
    before = (insertion-1).clip(0, len(log_times)-1)
    after = insertion.clip(0, len(log_times)-1)
    before_error = np.abs(log_times[before]-times)
    after_error = np.abs(log_times[after]-times)
    take_before = (before_error < after_error) | (
        (before_error == after_error) & (first_rows[before] < first_rows[after])
    )
    nearest = np.where(take_before, before, after)
    error = np.abs(log_times[nearest]-times)
    if (error > 1e-7).any():
        raise ValueError(f'Sample time {times[error.argmax()]:g} has no matching load time in {path}')
    return values[first_rows[nearest], 1]


def write_gif_frame(stream, palette, previous, duration):
    encoding = {'duration': duration, 'disposal': 1}
    if previous is None:
        header, _ = GifImagePlugin.getheader(palette, info={'loop': 0})
        stream.writelines(header)
    current = palette.convert('RGB')
    if previous is not None:
        difference = ImageChops.difference(previous, current)
        # A one-pixel transparent frame preserves the duration of a duplicate.
        bounds = difference.getbbox() or (0, 0, 1, 1)
        palette = palette.crop(bounds)
        red, green, blue = difference.crop(bounds).split()
        changed = ImageChops.lighter(ImageChops.lighter(red, green), blue)
        changed = changed.point([0] + [255]*255)
        # A transparent index must not be used by any pixel that changed.
        for index, count in enumerate(palette.histogram(mask=changed)):
            if count == 0:
                palette.paste(index, mask=ImageChops.invert(changed))
                encoding['transparency'] = index
                break
        encoding.update(offset=bounds[:2], include_color_table=True)
    stream.writelines(GifImagePlugin.getdata(palette, **encoding))
    return current


def write_gif(fig, draw_frame, count, path, final_path, fps):
    # Keep only adjacent frames to find changes and write the GIF incrementally.
    fig.set_dpi(91)
    duration = int(1000/fps)
    previous = None
    with path.open('wb') as stream:
        for index in range(count):
            draw_frame(index)
            if index == count-1:
                fig.savefig(final_path, dpi=140)
            pixels, size = fig.canvas.print_to_buffer()
            image = Image.frombuffer('RGBA', size, pixels, 'raw', 'RGBA', 0, 1)
            palette = image.convert('RGB').quantize()
            previous = write_gif_frame(stream, palette, previous, duration)
        stream.write(b';')


def render(case,out):
    meta = json.loads((case/'scenario.json').read_text())
    root = case/'postProcessing/cinema'
    dirs = sorted([p for p in root.iterdir() if p.is_dir() and (p/'bed.vtp').is_file()],key=lambda p:float(p.name))
    if not dirs: raise RuntimeError('No sampled times found')
    times = np.asarray([float(p.name) for p in dirs])
    if not np.isfinite(times).all():
        raise ValueError('Nonfinite cinema sample time')
    for folder in dirs:
        for name in SURFACES:
            path = folder/f'{name}.vtp'
            if not path.is_file():
                raise FileNotFoundError(path)
    hullpath = case/'constant/triSurface/shipwreck.stl'
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(hullpath))
    reader.Update()
    hull = reader.GetOutput()
    sections = [section(hull,h) for h in (.45,.2,0)]
    shutil.copyfile(hullpath,out/'shipwreck.stl')
    pvd = ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian')
    collection = ET.SubElement(pvd,'Collection')
    loadpath = case/'postProcessing/shipwreckLoads/0/loads.csv'
    impactpath = case/'postProcessing/shipwreckImpact/0/impact.csv'
    for name,src in [('shipwreck_loads.csv',loadpath),('sediment_budget.csv',case/'postProcessing/sedimentBudget/budget.csv'),
                     ('suspended_impact.csv',impactpath),('scenario.json',case/'scenario.json')]:
        shutil.copyfile(src,out/name)
    loads = sample_loads(loadpath, times, 'water_Fx_N')
    impact = sample_loads(impactpath, times, 'impact_Fx_N')
    time_history = []
    limits = {'z':0.,'cs':0.,'u':0.}
    for p,t,water_force,impact_force in zip(dirs,times,loads,impact):
        samples = [vtk_sample(p/f'{name}.vtp') for name in SURFACES]
        row = {'time_s':float(t),**bed_frame(samples[2]),
               'water_Fx_N':float(water_force),'suspended_impact_Fx_N':float(impact_force)}
        time_history.append(row)
        for key,s,values in [('u',samples[0],np.linalg.norm(samples[0][2]['U'],axis=1)),
                             ('cs',samples[1],samples[1][2]['Cs']*meta['assumed_sediment_density_kg_m3']*1000)]:
            valid = s[2]['shipwreckMask']<.5
            limits[key] = max(limits[key],float(values[valid].max()))
        limits['z'] = max(limits['z'],row['max_scour_mm'],row['max_deposition_mm'])
        destination = Path('surfaces')/p.name
        (out/destination).mkdir(parents=True)
        for part,name in enumerate(SURFACES):
            relative = destination/f'{name}.vtp'
            shutil.copyfile(p/f'{name}.vtp',out/relative)
            ET.SubElement(collection,'DataSet',timestep=p.name,group=name,part=str(part),file=relative.as_posix())
    with (out/'shipwreck_series.pvd').open('wb') as stream:
        ET.ElementTree(pvd).write(stream, encoding='utf-8', xml_declaration=True)
    with (out/'scour_history.csv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(time_history[0]))
        writer.writeheader();writer.writerows(time_history)

    fig,axes = plt.subplots(2,2,figsize=(12,7.2),layout='constrained')
    axes = axes.flatten()
    flow_norm = Normalize(0,max(meta['inlet_speed_m_s'],limits['u']))
    cs_norm = Normalize(0,max(1.,limits['cs']))
    z = max(.001,limits['z'])
    bed_norm = TwoSlopeNorm(vmin=-z,vcenter=0,vmax=z)
    norms = (flow_norm,cs_norm,bed_norm)
    cmaps = ('cividis','magma','RdBu_r')
    titles = ('Water speed and streamlines at z = 0.45 m','Suspended sand at z = 0.20 m','Bed change relative to flat start')
    units = ('m/s','mg/L','mm')
    for ax,norm,cmap,unit in zip(axes,norms,cmaps,units):
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,label=unit,shrink=.8)
    px,py = np.meshgrid(np.linspace(ROI[0],ROI[1],75),np.linspace(ROI[2],ROI[3],55))
    tvals = np.asarray([r['time_s'] for r in time_history])
    scour = [r['max_scour_mm'] for r in time_history];deposition = [r['max_deposition_mm'] for r in time_history]

    def frame(n):
        p = dirs[n]; row = time_history[n]
        frames = [vtk_sample(p/f'{name}.vtp') for name in SURFACES]
        values = [np.linalg.norm(frames[0][2]['U'],axis=1),
                  frames[1][2]['Cs']*meta['assumed_sediment_density_kg_m3']*1000,
                  frames[2][0][frames[2][1],2].mean(1)*1000]
        for ax,sample,value,outline,cmap,norm,title in zip(axes,frames,values,sections,cmaps,norms,titles):
            ax.clear()
            draw_surface(ax,sample,value,sample[2]['shipwreckMask']>=.5,cmap,norm)
            ax.add_collection(LineCollection(outline,linewidths=1.0,colors='#53d5eb' if ax!=axes[2] else '#222222'))
            ax.set(xlim=ROI[:2],ylim=ROI[2:],aspect='equal',xlabel='Downstream x (m)',ylabel='Cross-stream y (m)',title=title)
        flow_lines(axes[0],frames[0],px,py)
        ax=axes[3];ax.clear()
        ax.plot(tvals[:n+1],scour[:n+1],color='#1d60a5',label='Maximum scour depth')
        ax.plot(tvals[:n+1],deposition[:n+1],color='#ba4b30',label='Maximum deposition')
        ax.set(xlim=(0,tvals[-1]),ylim=(0,max(z*1.08,.001)),xlabel='Physical time (s)',ylabel='Depth or height (mm)')
        ax.grid(alpha=.25);ax.legend(loc='upper left',fontsize=8)
        ax.set_title(f"Water Fx {row['water_Fx_N']:.3g} N | Sand impact Fx {row['suspended_impact_Fx_N']:.3g} N")
        fig.suptitle(f"5 m inverted ship | 90° to current | t = {row['time_s']:.3g} s",fontweight='bold')

    fps = 5
    try:
        write_gif(fig, frame, len(dirs), out/'shipwreck_scour.gif', out/'final_frame.png', fps)
    finally:
        plt.close(fig)
    status = {'physical_end_s':float(tvals[-1]),'frames':len(dirs),'display_fps':fps,'colour_scales':{'water_m_s':flow_norm.vmax,
              'suspended_sand_mg_L':cs_norm.vmax,'bed_limits_mm':[-z,z]},'initial_bed_z_m':0.,
              'last':time_history[-1],'geometry_sha256':meta['solver_geometry_sha256'],
              'volume_scope':'Exposed bed triangles (shipwreckMask < 0.5) in computational domain; volume is bulk bed change'}
    summary = json.dumps(status, indent=2, allow_nan=False)+'\n'
    (out/'animation_summary.json').write_bytes(summary.encode('utf-8'))
    return status


def main(case,out):
    if out.exists():
        raise FileExistsError(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f'.{out.name}.', dir=out.parent) as directory:
        staging = Path(directory)/out.name
        staging.mkdir()
        status = render(case, staging)
        if out.exists():
            raise FileExistsError(out)
        staging.rename(out)
    print(json.dumps(status,indent=2,allow_nan=False))


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
