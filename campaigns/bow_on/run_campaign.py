"""Run the independent bow-on cases in order after the current CFD job ends."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time


parser = argparse.ArgumentParser()
parser.add_argument('root',type=Path)
parser.add_argument('--model',type=Path,required=True)
parser.add_argument('--baseline',type=Path,required=True)
parser.add_argument('--prepare-only',action='store_true')
args = parser.parse_args()
ROOT = args.root.resolve()
MODEL = args.model.resolve()
BASE = args.baseline.resolve()
REPORT = ROOT/'reports/campaign.json'
U = 4.28*0.3048
DEPTH = 31.5*0.3048
CS = (645*0.15/1000)/2650
X = [-10.,-2.5,2.5,15.]
Y = [-7.5,-3.,3.,7.5]


def atomic(path, text):
    tmp = path.with_name('.'+path.name+'.tmp')
    tmp.write_text(text,encoding='utf-8')
    tmp.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def active_solvers():
    names = ('pimpleFoam','sedExnerFoam','shipwreckSedEx')
    pids = []
    for p in Path('/proc').iterdir():
        if p.name.isdecimal():
            try:
                if (p/'comm').read_text().strip().startswith(names) and (p/'stat').read_text().rsplit(')',1)[1].split()[0]!='Z':
                    pids.append(int(p.name))
            except (OSError,IndexError):
                pass
    return sorted(pids)


def update_mesh(path):
    text = path.read_text()
    coords = ''.join(f'({x:.16g} {y:.16g} {z:.16g})\n' for z in (0,DEPTH) for y in Y for x in X)
    text,n = re.subn(r'vertices \(\n.*?\n\);',lambda m:'vertices (\n'+coords+');',text,count=1,flags=re.S)
    if n!=1: raise RuntimeError('Unexpected blockMeshDict vertices')
    matches = list(re.finditer(r'^hex \([^)]*\) \([^)]*\) simpleGrading \([^)]*\)$',text,re.M))
    if len(matches)!=9: raise RuntimeError('Expected 9 conformal blocks from the existing generator')
    nx,ny = [10,24,16],[8,32,8]
    gx,gy = [.42,1,3.8],[.42,1,2.4]
    for n,m in reversed(list(enumerate(matches))):
        j,i = divmod(n,3)
        indices = m.group().split(')')[0]+')'
        block = f'{indices} ({nx[i]} {ny[j]} 24) simpleGrading ({gx[i]} {gy[j]} 25)'
        text = text[:m.start()]+block+text[m.end():]
    atomic(path,text)


def prepare(tag,end,checkpoint,sampling,hull):
    case = ROOT/'cases'/tag
    if case.exists(): raise FileExistsError(case)
    flow = ROOT/'reports'/f'{tag}_inflow.csv'
    atomic(flow,f'time_s,speed_m_s,Cs_volume_fraction\n0,{U:.16g},{CS:.16g}\n{end},{U:.16g},{CS:.16g}\n')
    subprocess.run([str(MODEL/'build/bin/shipwreck-case'),'--case',str(case),'--template',
        str(BASE/'cases/caisson_prism_cpu_002_20261118'),'--stl',hull['selected'],'--flow',str(flow),
        '--length','5','--yaw','0','--burial',repr(hull['burial_m']),'--depth',repr(DEPTH),
        '--domain-span','6','--interface-thickness','.16','--end',str(end),'--dt','.004',
        '--write-interval',str(checkpoint),'--d50','.00042','--initial-cs',repr(CS),
        '--impact-model','impulse','--restitution','.5','--friction','.3'],check=True)
    update_mesh(case/'system/blockMeshDict')
    control = case/'system/controlDict'
    original = control.read_text()
    text = original[:original.index('functions')].replace('maxDeltaT 0.004;','maxDeltaT 0.02;').replace('purgeWrite 2;','purgeWrite 0;')
    if 'maxDeltaT 0.02;' not in text or 'purgeWrite 0;' not in text:
        raise RuntimeError('Prepared controlDict disagrees with output schedule')
    text += f'''functions
{{
    bedShear {{ type wallShearStress; libs (fieldFunctionObjects); patches (bed); writeControl writeTime; }}
    cinema
    {{
        type surfaces;
        libs (sampling);
        writeControl adjustableRunTime;
        writeInterval {sampling};
        fields (U Cs shipwreckMask);
        surfaceFormat vtk;
        formatOptions {{ vtk {{ format binary; legacy false; }} }}
        sampleScheme cell;
        interpolationScheme cellPoint;
        surfaces
        {{
            bed {{ type patch; patches (bed); triangulate true; interpolate false; }}
            nearBed {{ type cuttingPlane; point (0 0 .20); normal (0 0 1); triangulate true; interpolate false; }}
            flow {{ type cuttingPlane; point (0 0 .45); normal (0 0 1); triangulate true; interpolate false; }}
        }}
    }}
}}
'''
    atomic(control,text)
    meta = json.loads((case/'scenario.json').read_text())
    meta.update(cells=57600, length_m=5,roll_deg=180,bow_direction=[-1,0,0],inflow_direction=[1,0,0],
                physical_time_s=end,max_dt_s=.02,display_fps=1,
                source_geometry_sha256=hull['original_sha256'],solver_geometry_sha256=hashlib.sha256((case/'constant/triSurface/shipwreck.stl').read_bytes()).hexdigest(),
                simplified_geometry_sha256=hull['sha256'],domain_m={'x':X,'y':Y,'z':[0,DEPTH]},
                checkpoint_interval_s=checkpoint,sampling_interval_s=sampling,
                inlet_speed_m_s=U,inlet_Cs_volume_fraction=CS,mean_depth_m=DEPTH,
                measured_source='https://pubs.usgs.gov/wsp/1819j/report.pdf',measurement_date='1961-04-21',
                measurement_tables=[1,7,8],bed_material_d50_m=.00042,measurement_total_suspended_ppm=645,
                measured_suspended_sand_mass_fraction=.15,assumed_sediment_density_kg_m3=2650,
                geometry_role='5 m hypothetical fixed hull with deck reference on a flat sand bed')
    atomic(case/'scenario.json',json.dumps(meta,indent=2,allow_nan=False)+'\n')
    return case


def monitor(case,state,save):
    if active_solvers(): raise RuntimeError('The six-core solver slot is occupied')
    if shutil.disk_usage(case).free < 8*2**30: raise RuntimeError('Less than 8 GiB free on the case drive')
    with (ROOT/'logs'/f'{case.name}.controller.log').open('x') as log:
        process = subprocess.Popen(['bash',str(MODEL/'run.sh'),str(BASE),str(case),'6'],stdout=log,stderr=subprocess.STDOUT)
        state.update(worker_pid=process.pid,phase='computing',case=str(case),last_time_s=0.)
        save()
        while process.poll() is None:
            time.sleep(5)
            solver = case/'log.solver'
            if not solver.is_file(): continue
            with solver.open('rb') as buffer:
                buffer.seek(max(0,buffer.seek(0,2)-100000))
                tail = buffer.read().decode(errors='replace')
            times = re.findall(r'^Time = (\S+)\s*$',tail,re.M)
            if times: state['last_time_s'] = float(times[-1])
            co = re.findall(r'^Courant Number[^\n]* max: (\S+)',tail,re.M)
            fatal = re.search(r'FOAM FATAL|segmentation fault|\b(?:nan|inf)\b',tail,re.I)
            if not state.get('stop_reason') and (shutil.disk_usage(case).free<4*2**30 or fatal or any(float(c)>1 for c in co)):
                state['stop_reason']='Resource or numerical alarm; write checkpoint and diagnose'
                control = case/'system/controlDict'
                atomic(control,re.sub(r'\bstopAt\s+[^;]+;','stopAt writeNow;',control.read_text(),count=1))
            save()
        if process.returncode or state.get('stop_reason'):
            raise RuntimeError(f'{case.name} stopped; see preserved solver and controller logs')
    summary = json.loads((case/'run_summary.json').read_text())
    if not summary['numerical_checks_passed']:
        raise RuntimeError(f'{case.name} failed numerical checks')
    view = case/'postProcessing/cinema'
    present = sorted(float(p.name) for p in view.iterdir() if p.is_dir() and (p/'bed.vtp').is_file()
                     and (p/'flow.vtp').is_file() and (p/'nearBed.vtp').is_file())
    meta = json.loads((case/'scenario.json').read_text())
    expected = [i*meta['sampling_interval_s'] for i in range(1,round(meta['physical_time_s']/meta['sampling_interval_s'])+1)]
    if len(present)!=len(expected) or any(abs(x-y)>1e-7 for x,y in zip(present,expected)):
        raise RuntimeError(f'{case.name} has incomplete or at unexpected times cinema sampling: {present}')
    summary['sample_times_s'] = present
    return summary


def main():
    with (ROOT/'reports/.campaign.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if REPORT.exists(): raise FileExistsError('An existing campaign requires diagnosis before reuse')
        def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
        state = {'status':'running','phase':'preparing','controller_pid':os.getpid(),'started_utc':now(),'completed':{},
                 'provenance':{str(p):sha(p) for p in [MODEL/'build/bin/shipwreckSedExnerFoam',MODEL/'build/bin/shipwreck-case',
                     MODEL/'build/lib/libshipwreck_impact.a',MODEL/'run.sh',Path(__file__).resolve()]}}
        def save():
            state['updated_utc'] = now()
            atomic(REPORT,json.dumps(state,indent=2,allow_nan=False)+'\n')
        save()
        try:
            hull = json.loads((ROOT/'reports/hull_geometry.json').read_text())
            if hull['status']!='passed' or sha(Path(hull['selected']))!=hull['sha256']:
                raise RuntimeError('Verified hull required')
            schedule = [('gate_0p1',.1,.05,.05),('preview_2s',2.,.5,1.),('bow_on_600s',600.,2.,1.)]
            for tag,end,checkpoint,sampling in schedule:
                case = prepare(tag,end,checkpoint,sampling,hull)
                state['phase'] = 'prepared';state['case'] = str(case);save()
                if tag=='gate_0p1':
                    subprocess.run([sys.executable,str(ROOT/'pose_preview.py'),str(case),str(ROOT/'deliverables/hull_pose.png')],check=True)
                if args.prepare_only: continue
                while active_solvers():
                    state['phase'] = 'waiting_for_existing_six_core_CFD'
                    state['waiting_on_pids'] = active_solvers();save();time.sleep(15)
                state.pop('waiting_on_pids',None)
                for path,digest in state['provenance'].items():
                    if sha(Path(path))!=digest: raise RuntimeError(f'Provenance changed: {path}')
                state['completed'][tag] = monitor(case,state,save);save()
                if tag in ('preview_2s','bow_on_600s'):
                    state['phase']='rendering';save()
                    subprocess.run([sys.executable,str(ROOT/'render.py'),str(case),str(ROOT/'deliverables'/tag)],check=True)
            state.update(status='prepared' if args.prepare_only else 'completed',phase='done');save()
        except Exception as exc:
            state.update(status='needs_diagnostic',error=repr(exc));save();raise


if __name__=='__main__':
    main()
