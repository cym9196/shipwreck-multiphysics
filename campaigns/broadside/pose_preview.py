"""Technical views of the exact placed computational hull, independent of CFD."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def main(case,output):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(case/'constant/triSurface/shipwreck.stl'))
    reader.Update()
    poly = reader.GetOutput()
    pts = vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
    faces = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1,3)
    xyz = pts[faces]
    normals = np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
    light = normals@np.array([-.48,-.30,.823]) / np.linalg.norm(normals,axis=1)
    bright = .58+.35*np.clip(light,0,1)
    above = xyz[:,:,2].mean(1)>=0
    base = np.where(above[:,None],np.array([.19,.61,.74]),np.array([.53,.54,.50]))
    rgb = base*bright[:,None]+.15
    fig,axes = plt.subplots(1,2,figsize=(10.3,4),layout='constrained')
    specs = [(0,1,2,(-2.1,2.1),(-2.8,2.8),'Plan view'),
             (1,2,0,(-2.8,2.8),(-.75,1.02),'View along the incoming flow')]
    for ax,(a,b,d,xlim,ylim,title) in zip(axes,specs):
        order = np.argsort(xyz[:,:,d].mean(1))
        if d==0: order=order[::-1]
        ax.add_collection(PolyCollection(xyz[order][:,:,[a,b]],facecolors=rgb[order],edgecolors='none',rasterized=True))
        ax.set(aspect='equal',xlim=xlim,ylim=ylim,title=title,xlabel='xyz'[a]+' (m)',ylabel='xyz'[b]+' (m)')
        ax.grid(alpha=.2)
    axes[0].annotate('Flow toward +x',xy=(-.99,0),xytext=(-2.06,0),
                     arrowprops={'arrowstyle':'->','color':'#183c68','lw':2},color='#183c68',fontsize=9)
    axes[1].axhline(0,color='#594d41',linestyle='--',linewidth=1)
    axes[1].annotate('Initial bed',xy=(1.9,0),xytext=(1.45,-.28),arrowprops={'arrowstyle':'->'},fontsize=9)
    fig.suptitle('Computational hull placement | Length 5 m | Upside down | Broadside')
    fig.savefig(output,dpi=160)
    plt.close(fig)


if __name__=='__main__':
    main(Path(sys.argv[1]),Path(sys.argv[2]))
