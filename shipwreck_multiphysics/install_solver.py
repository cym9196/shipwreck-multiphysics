"""Assemble the independent C++ solver with guarded edits to the baseline."""
from hashlib import sha256
from pathlib import Path
import json
import sys
from tempfile import TemporaryDirectory


def assemble(base, project):
    source = base / 'src/restart_candidate_20261115'
    paths = sorted((*source.glob('*.H'), *source.glob('*.C'), source / 'Make/options'))
    originals = {path.relative_to(source): path.read_bytes() for path in paths}
    files = dict(originals)
    files.update({Path(path.name): path.read_bytes() for path in sorted((project / 'src').glob('*.H'))})

    def edit(name, old, new):
        path = Path(name)
        if path not in files:
            raise ValueError(f'Baseline is missing {name}')
        text = files[path].decode('utf-8').replace('\r\n', '\n')
        if text.count(old) != 1:
            raise ValueError(f'{name}: incompatible baseline near {old!r}')
        files[path] = text.replace(old, new).encode('utf-8')

    files[Path('Make/files')] = b'sedExnerFoam.C\n\nEXE = $(FOAM_USER_APPBIN)/shipwreckSedExnerFoam\n'
    edit('Make/options', '    -I$(LIB_SRC)/sampling/lnInclude',
         '    -I$(LIB_SRC)/surfMesh/lnInclude \\\n    -I$(LIB_SRC)/sampling/lnInclude')
    edit('Make/options', '    -lmeshTools', '    -lsurfMesh \\\n    -lmeshTools')
    edit('Make/options', '    -lsettlingModels',
         '    -L$(SHIPWRECK_IMPACT_LIB) -lshipwreck_impact -ldl -lpthread \\\n    -lsettlingModels')
    edit('sedExnerFoam.C', '#include "OFstream.H"', '#include "OFstream.H"\n#include "Shipwreck.H"')
    edit('sedExnerFoam.C', '    #include "createSedimentBudget.H"', """    #include "createSedimentBudget.H"
    if (switchSuspension != "on")
        FatalErrorInFunction << "Shipwreck load output requires suspension on" << exit(FatalError);
    Shipwreck shipwreck(mesh,rhoS.value(),dS.value());""")
    edit('sedExnerFoam.C', '                    MRF.update();',
         '                    MRF.update();\n                    shipwreck.update();')
    edit('sedExnerFoam.C', '        #include "writeSedimentBudget.H"', """        #include "writeSedimentBudget.H"
        shipwreck.writeLoads(U, CsPtr(), WsPtr(), rhoF.value(), rhoS.value());""")
    edit('UEqn.H', '  + turbulence->divDevReff(U)',
         '  + turbulence->divDevReff(U)\n  + fvm::Sp(shipwreck.inversePenalisationTime(), U)')
    edit('bedShearStress.H', '    if (filterShields)',
         '    shipwreck.shieldBed(patchID, ssp);\n\n    if (filterShields)')

    manifest = {path.as_posix(): sha256(content).hexdigest() for path, content in originals.items()}
    outputs = {Path('solver') / path: content for path, content in files.items()}
    outputs[Path('baseline_sources.json')] = (json.dumps(manifest, indent=2) + '\n').encode('utf-8')
    publish(project / 'build', outputs)


def publish(build, outputs):
    build.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='.solver-', dir=build) as stage_name:
        stage = Path(stage_name)
        changed = []
        for path, content in outputs.items():
            destination = build / path
            if destination.is_file() and destination.read_bytes() == content:
                continue
            staged = stage / path
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.write_bytes(content)
            changed.append(path)

        solver = build / 'solver'
        for path in (*solver.glob('*.H'), *solver.glob('*.C')):
            if path.relative_to(build) not in outputs:
                path.unlink()
        for path in changed:
            destination = build / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            (stage / path).replace(destination)


if __name__ == '__main__':
    assemble(*(Path(arg).resolve() for arg in sys.argv[1:]))
