from pymongo import MongoClient
import os
import numpy as np
import pybel
from rdkit import Chem
import pickle
from tqdm import tqdm
from io import StringIO
from job_runner import update_conf_positions


def MongoUploader(mol, coll_name='politzer'):
    MONGO_URL = os.environ.get('MONGO_URL')
    with StringIO() as f:
        w = Chem.SDWriter(f)
        w.write(mol)
        w.close()
        sdf = f.getvalue()
    assert hasattr(mol, '_id')
    mol.sdf = sdf
    with MongoClient(MONGO_URL) as conn:
        coll = conn['bigcsd'][coll_name]
        coll.update_one({'_id': mol._id}, {'$set': mol.__dict__}, upsert=True)


def MongoPuller(coll_name='conf3', all_converged=True):
    MONGO_URL = os.environ.get('MONGO_URL')
    with MongoClient(MONGO_URL) as conn:
        coll = conn['bigcsd'][coll_name]
        ids = coll.distinct('inchikey', {'converged': True})
        for _id in tqdm(ids, smoothing=0):
            dlist = list(coll.find({'inchikey': _id}, projection=['conf_id', 'labels', 'xyz1', 'converged', 'energy']))
            if all_converged and not all(d.get('converged') for d in dlist):
                continue
            for i in range(min(10, len(dlist))):
                mol = lapo2rdmol(dlist[i]['labels'], dlist[i]['xyz1'])
                if mol is not None:
                    break
            if not mol:
                print('Failed:', _id)
                continue
            mol.energies = [d['energy'] for d in dlist]
            mol._id = _id
            for d in dlist:
                c = Chem.Conformer(len(d['labels']))
                c.SetId(d['conf_id'])
                mol.AddConformer(c)
                update_conf_positions(mol, d['conf_id'], np.asarray(d['xyz1']).reshape(-1, 3))
            with open('mols/' + _id, 'wb') as f:
                pickle.dump(mol, f)


def lapo2rdmol(labels, positions):
    xyz = '{}\n\n'.format(len(labels))
    positions = np.asarray(positions).reshape(-1, 3)
    for l, p in zip(labels, positions):
        xyz += '{} {} {} {}\n'.format(l, *p)
    pmol = pybel.readstring('xyz', xyz)
    rdmol = Chem.MolFromPDBBlock(pmol.write('pdb'), sanitize=False, removeHs=False)
    if not rdmol:
        return None
    rdmol.RemoveAllConformers()
    return rdmol
