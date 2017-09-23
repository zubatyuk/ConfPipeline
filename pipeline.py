import yaml
from rdkit_utils import EmbedConformers, CleanConformers
from job_runner import TurboOpt, TurboSP, TurboDesc
from rdkit import Chem
from rdkit.Chem import AllChem

def pipeline(mol, config):
    imax = len(config)
    for i in range(imax):
        k, v = list(config[i].items())[0]
        mol = eval(k)(mol, **v)
        print('>>', len(mol.energies), mol.energies)
    return mol


if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.description = 'Morse encoder/decoder, uses extended table (all printable ASCII characters)'
    parser.add_argument('--infile', '-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='Input SDF file (default STDIN)')
    parser.add_argument('--outfile', '-o', nargs='?', type=argparse.FileType('r'), default=sys.stdout,
                        help='Output SDF file (default STDOUT)')
    parser.add_argument('--decode', '-d', action='store_true', default=False,
                        help='Run decoder (encoder by default)')

    config = yaml.load(open('pipeline.yml').read())


    mol = Chem.MolFromSmiles('COCO')
    mol = pipeline(mol, config)
    print(len(mol.energies), mol.energies)
    writer = Chem.SDWriter('out.sdf')
    for c in mol.GetConformers():
        writer.write(mol, c.GetId())
    writer.close()

