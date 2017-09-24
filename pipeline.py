import yaml
from rdkit_utils import EmbedConformers, CleanConformers
from job_runner import TurboOpt, TurboSP, TurboDesc
from rdkit import Chem
from rdkit.Chem import AllChem
from mongo_uploader import MongoUploader


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
    import pickle

    parser = argparse.ArgumentParser()
    parser.description = 'Conformer generation pipeline'
    parser.add_argument('--infile', '-i', type=argparse.FileType('rb'),
                        help='Input SDF file')
    parser.add_argument('--config', '-c', nargs='?', type=argparse.FileType('r'), default='pipeline.yml',
                        help='Pipeline config (default pipeline.yml)')

    ns = parser.parse_args()

    config = yaml.load(ns.config.read())
    mol = pickle.load(ns.infile)
    run_though_pipeline(mol, config)
