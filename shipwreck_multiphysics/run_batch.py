"""Prepare and run independent hull cases from a declarative JSON batch."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess


REQUIRED = {'stl', 'flow', 'length', 'depth', 'end'}
OPTIONS = REQUIRED | {'beam','height','yaw','burial','x','y','dt','write_interval','d50',
                      'initial_cs','refine','domain_span','interface_thickness','penalty_time',
                      'impact_model','restitution','friction'}


def command_options(row, directory):
    if set(row)-OPTIONS or REQUIRED-set(row):
        raise ValueError(f'Unknown or missing case options: {set(row)-OPTIONS}, {REQUIRED-set(row)}')
    arguments = []
    for key, value in row.items():
        if key in ('stl','flow'):
            if not isinstance(value,str):
                raise ValueError(f'{key} must be a file path')
            value = (directory/value).resolve(strict=True)
            if not value.is_file():
                raise ValueError(f'{key} must be a file')
        elif key == 'impact_model':
            if value not in ('off','impulse'):
                raise ValueError('impact_model must be off or impulse')
        elif type(value) not in (int,float) or not math.isfinite(value):
            raise ValueError(f'{key} must be a finite number')
        arguments.extend(['--'+key.replace('_','-'),str(value)])
    return arguments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    parser.add_argument('--base',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    base = args.base.resolve(strict=True)
    config = args.config.resolve(strict=True)
    output = args.output.resolve()
    plan = json.loads(config.read_text(encoding='utf-8'))
    if (not isinstance(plan, dict) or set(plan) != {'schema', 'template_case', 'ranks', 'defaults', 'cases'}
            or type(plan['schema']) is not int or plan['schema'] != 1):
        raise ValueError('Expected batch schema 1 with template_case, ranks, defaults and cases')
    if not isinstance(plan['template_case'], str) or not isinstance(plan['defaults'], dict):
        raise ValueError('Expected a template_case path and a defaults object')
    if (type(plan['ranks']) is not int or plan['ranks'] not in (1, 6)
            or not isinstance(plan['cases'], list) or not plan['cases']):
        raise ValueError('Provide a nonempty cases array and choose 1 or 6 ranks')
    template = (base/plan['template_case']).resolve(strict=True)
    if not template.is_relative_to(base) or not template.is_dir():
        raise ValueError('Template must be a directory in the baseline project')
    if output.exists():
        raise FileExistsError('Batch output directory already exists')
    expected = []
    reserved_names = {'batch_summary.json'}
    for item in plan['cases']:
        name = item.get('name') if isinstance(item, dict) else None
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,63}', name):
            raise ValueError(f'Invalid case name {name!r}')
        case_paths = {name.casefold(), f'{name}.prepare.log'.casefold(), f'{name}.solve.log'.casefold()}
        if case_paths & reserved_names:
            raise ValueError(f'Case name {name!r} collides with another case or batch file')
        reserved_names.update(case_paths)
        settings = plan['defaults'] | {k: v for k, v in item.items() if k != 'name'}
        expected.append((name, command_options(settings, config.parent)))
    output.mkdir(parents=True,exist_ok=False)
    state = {'config':str(config),'prepare_only':args.prepare_only,'status':'preparing','cases':{}}

    def save():
        path = output/'.batch_summary.tmp'
        path.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        path.replace(output/'batch_summary.json')

    def step(name,phase,argv):
        state['status'] = phase
        state['active_case'] = name
        save()
        with (output/f'{name}.{phase}.log').open('x') as log:
            return subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)

    try:
        for name,parameters in expected:
            directory = output/name
            step(name,'prepare',[str(project/'build/bin/shipwreck-case'),'--case',str(directory),
                                 '--template',str(template),*parameters]).check_returncode()
            state['cases'][name] = {'prepared':True}
            save()
        if not args.prepare_only:
            for name,_ in expected:
                directory = output/name
                process = step(name, 'solve', ['bash', str(project/'run.sh'), str(base), str(directory), str(plan['ranks'])])
                report = directory/'run_summary.json'
                if report.is_file():
                    state['cases'][name].update(json.loads(report.read_text(encoding='utf-8')))
                process.check_returncode()
                if not report.is_file():
                    raise FileNotFoundError(f'Solver finished without writing {report}')
                save()
        state.update(status='prepared' if args.prepare_only else 'completed',active_case=None)
        save()
    except (Exception, KeyboardInterrupt) as exc:
        message = 'Interrupted by user' if isinstance(exc, KeyboardInterrupt) else str(exc)
        state.update(status='needs_diagnostic', error=message)
        save()
        raise
    print(output/'batch_summary.json')


if __name__=='__main__':
    main()
