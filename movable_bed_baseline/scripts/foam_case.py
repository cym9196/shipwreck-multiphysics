"""Small OpenFOAM file writers and checked command execution."""
from pathlib import Path
import os
import re
import subprocess


def header(name, cls='dictionary'):
    return f'FoamFile {{version 2.0; format ascii; class {cls}; object {name};}}\n'


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix=='.json':
        temporary=path.with_name('.'+path.name+f'.{os.getpid()}.tmp')
        temporary.write_text(text,encoding='utf-8')
        os.replace(temporary,path)
    else:
        path.write_text(text, encoding='utf-8')


def field(name, dimensions, boundary, internal='uniform 0', cls='volScalarField'):
    return (header(name, cls) + f'dimensions {dimensions};\ninternalField {internal};\n'
            + 'boundaryField\n{\n' + ''.join(f'{k} {{ {v} }}\n' for k, v in boundary.items()) + '}\n')


def set_entry(text, name, value):
    result, count = re.subn(r'(\b' + re.escape(name) + r'\s+)[^;]+;',
                            lambda m: m.group(1) + str(value) + ';', text, count=1)
    if count != 1:
        raise ValueError(f'Missing entry: {name}')
    return result


def run(case, args, label=None, timeout=300):
    case = Path(case)
    logfile = case / ('log.' + (label or args[0]))
    with logfile.open('w') as stream:
        done = subprocess.run(args, cwd=case, env={**os.environ, 'PWD': str(case)},
                              stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
    if done.returncode:
        raise RuntimeError(f'{args[0]} failed ({done.returncode}); see {logfile}')
    return logfile
