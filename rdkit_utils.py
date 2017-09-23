from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def CleanConformers(mol, max_keep=1, rmsd_thresh=0.5, max_energy=10.0, max_delta_e_compare=0.5):
    assert isinstance(mol, Chem.Mol)
    N = mol.GetNumConformers()
    if not N:
        # nothing to do
        return
    conformers = mol.GetConformers()
    has_energies = hasattr(mol, 'energies')
    # 1000 is the safe max for conformer energy difference
    energies = mol.energies if has_energies else [i * 1000 for i in range(N)]
    d_energies = np.subtract(energies, min(energies))
    keep = [0] * len(conformers) 
    # sorted by energy
    for iconf in np.argsort(energies):
        _id = conformers[iconf].GetId()
        # max number of conformers reached -> delete
        if d_energies[iconf] > max_energy:
            mol.RemoveConformer(_id)
            continue
        # make list of conformer ids in keep within max_delta_e_compare window
        cids = [conformers[x].GetId() for x in range(len(keep)) if keep[x] and (energies[iconf] - energies[x]) < max_delta_e_compare]
        if sum(keep) != max_keep and is_uniq_conformer(mol, _id, cids, rmsd_thresh):
            keep[iconf] = 1
        else:
            mol.RemoveConformer(_id)
    # update energies
    if has_energies:
        mol.energies = [e for e, i in zip(energies, keep) if i]
    return mol


def OptimizeConformersFF(mol, ff='uff'):
    assert ff in ('uff', 'mmff')
    mol.energies = list()
    for conf in mol.GetConformers():
        if ff == 'uff':
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())
        else:
            AllChem.MMFFSanitizeMolecule(mol)
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf.GetId())
        ff.Minimize()
        mol.energies.append(ff.CalcEnergy())
    return mol

def EmbedConformers(mol, max_out=-1, nrot_min=4, nrot_max=10, pool_mult=50, out_mult=4, emax=20.0):
    nrot = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    nrot = max(nrot_min, min(nrot_max, nrot))
    nout = max(max_out, nrot * out_mult)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=nrot*pool_mult, pruneRmsThresh=-1.)
    OptimizeConformersFF(mol)
    CleanConformers(mol, max_keep=nout, max_energy=emax)
    return mol


def is_uniq_conformer(mol, ref_cid, target_cids, rmsd_threshold=0.5):
    """
    Check if reference conformer is close to any of target conformers
    """
    uniq = True
    for t in target_cids:
        rmsd = AllChem.GetBestRMS(mol, mol, int(ref_cid), int(t))
        if rmsd < rmsd_threshold:
            uniq = False
            break
    return uniq


def save_mol(mol, filename):
    with open(filename, 'w') as f:
        f.write(mol.ToBinary())

def load_mol(mol, filename):
    with open(filename, 'w') as f:
        f.write(mol.ToBinary())



if __name__ == '__main__':
    mol = Chem.MolFromSmiles('COOCOCOON')
    Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol)
    OptimizeConformersFF(mol)
    print(mol.energies)
    CleanConformers(mol, max_keep=10)
    print(mol.energies)
