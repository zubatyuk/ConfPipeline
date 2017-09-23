import yaml
from rdkit_utils import EmbedConformers, CleanConformers
from job_runner import TurboOpt, TurboSP, TurboDesc
from rdkit import Chem
from rdkit.Chem import AllChem


def run_though_pipeline(obj, config):
    for k, v in config:
        obj = eval(k)(obj, **v)
    return obj

if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.description = 'Morse encoder/decoder, uses extended table (all printable ASCII characters)'
    parser.add_argument('--infile', '-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='Input SDF file (default STDIN)')
    parser.add_argument('--outfile', '-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output SDF file (default STDOUT)')
    parser.add_argument('--config', '-c', nargs='?', type=argparse.FileType('r'), default='pipeline.yml',
                        help='Pipeline config (default pipeline.yml)')

    ns = parser.parse_args()

    config = yaml.load(ns.config.read())

    for mol in Chem.SDMolSupplier(ns.infile):
        ns.outfile.write(run_though_pipeline(mol))
