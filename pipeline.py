import yaml
from rdkit_utils import EmbedConformers, CleanConformers
from job_runner import TurboOpt, TurboSP, TurboDesc
from rdkit import Chem
from rdkit.Chem import AllChem


def run_though_pipeline(obj, config):
    for c in config:
        k, v = list(c.items())[0]
        print('Running ', k)
        print(obj.GetNumConformers())
        print([x.GetId() for x in obj.GetConformers()])
        print(obj.__dict__.get('energies'))
        obj = eval(k)(obj, **v)
    return obj

if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.description = 'Conformer generation pipeline'
    parser.add_argument('--infile', '-i', type=argparse.FileType('r'),
                        help='Input SDF file')
    parser.add_argument('--outfile', '-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output SDF file (default STDOUT)')
    parser.add_argument('--config', '-c', nargs='?', type=argparse.FileType('r'), default='pipeline.yml',
                        help='Pipeline config (default pipeline.yml)')

    ns = parser.parse_args()

    config = yaml.load(ns.config.read())

    writer = Chem.SDWriter(ns.outfile)
    for mol in Chem.SDMolSupplier(ns.infile.name):
        mol = run_though_pipeline(mol, config)
        writer.write(mol)
