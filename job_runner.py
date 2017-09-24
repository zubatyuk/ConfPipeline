from tempdir import in_tempdir
import sys
import os
import numpy as np
from joblib import Parallel, delayed
from surface_desc import politzer
from rdkit_utils import RemoveConformers


def TurboOpt(mol, script, n_jobs=1):
    conf_ids = [c.GetId() for c in mol.GetConformers()]
    results = Parallel(n_jobs=n_jobs)(delayed(run_calc)(
        *[
            (write_xyz, {'mol': mol, 'conf_id': conf_id, 'fil': 'tm.xyz'}),
            (os.system, [script]),
            (get_tm_xyz_energy, {})
        ]
    ) for conf_id in conf_ids)
    # !!! remove conformer if calculation failed
    RemoveConformers(mol, [not res[2][1] for res in results])
    for res, conf_id in zip(results, conf_ids):
        if res[2][1]:
            update_conf_positions(mol, conf_id, res[2][0])
            update_conf_energy(mol, conf_id, res[2][1])
    return mol


def TurboDesc(mol, script_gas, script_cosmo=None, n_jobs=1):
    conf_ids = [c.GetId() for c in mol.GetConformers()]
    results = Parallel(n_jobs=n_jobs)(delayed(run_calc)(
        *[
            (write_xyz, {'mol': mol, 'conf_id': conf_id, 'fil': 'tm.xyz'}),
            (os.system, [script_gas]),
            (parse_cosmo, {}),
            (os.system if script_cosmo else noop, [script_cosmo]),
            (get_tm_energy if script_cosmo else noop, {}),
            (parse_cosmo, {}),
        ]
    ) for conf_id in conf_ids)
    # !!! remove conformer if calculation failed
    RemoveConformers(mol, [not res[2][0] for res in results])
    mol.politzer = list()
    for res in results:
        if res[2][0]:
            desc = dict(Vol=res[2][0], S_tot=np.sum(res[2][2]), **politzer(res[2][2], res[2][3]))
            if script_cosmo and res[5][0]:
                desc = dict(**desc,
                            **politzer(res[5][2], res[5][3], '_cosmo'),
                            **politzer(res[5][2], res[5][3] - res[2][3], '_diff')
                            )
            mol.politzer.append(desc)
    if script_cosmo:
        mol.energies_cosmo = [x[4] for x in results if x[2][0]]
    return mol


def TurboSP(mol, script, n_jobs=1):
    conf_ids = [c.GetId() for c in mol.GetConformers()]
    results = Parallel(n_jobs=n_jobs)(delayed(run_calc)(
        *[
            (write_xyz, {'mol': mol, 'conf_id': conf_id, 'fil': 'tm.xyz'}),
            (os.system, [script]),
            (get_tm_energy, {}),
        ]
    ) for conf_id in conf_ids)
    for res, conf_id in zip(results, conf_ids):
        update_conf_energy(mol, conf_id, res[2])
    return mol


def run_calc(*funcargs, basedir='.'):
    '''
    :param funcargs: of tuples of function, **kwargs 
    '''
    results = list() 
    with in_tempdir(basedir=basedir):
        for f, a in funcargs:
            if isinstance(a, dict):
                results.append(f(**a))
            else:
                results.append(f(*a))
    return results


def update_conf_positions(mol, conf_id, positions):
    conf = mol.GetConformer(conf_id)
    for i, x in enumerate(positions):
        conf.SetAtomPosition(i, x)


def update_conf_energy(mol, conf_id, energy):
    if not hasattr(mol, 'energies'):
        mol.energies = [None] * mol.GetNumConformers()
    for i, conf in enumerate(mol.GetConformers()):
        if conf.GetId() == conf_id:
            mol.energies[i] = energy


def write_xyz(mol, conf_id, fil=sys.stdout):
    la, po = conf2lapo(mol, conf_id)
    res = '{}\n\n'.format(len(la))
    for l, p in zip(la, po):
        res += '{} {} {} {}\n'.format(l, *p)
    if hasattr(fil, 'write'):
        fil.write(res)
    else:
        with open(fil, 'w') as f:
            f.write(res)


def conf2lapo(mol, conf_id):
    la = [a.GetSymbol() for a in mol.GetAtoms()]
    po = mol.GetConformer(conf_id).GetPositions()
    return la, po


def get_tm_xyz_energy():
    xyz = None
    ene = None
    try:
        with open('gradient') as fil:
            for line in fil:
                if line.startswith('  cycle ='):
                    xyz = list()
                    ene = float(line[32:52])
                    for line in fil:
                        ll = line.split()
                        if len(ll) != 4:
                            break
                        xyz.append([float(i) * 0.529177 for i in ll[:3]])
        result = np.asarray(xyz), ene
    except:
        result = (None, None)
    return result


def get_tm_energy():
    try:
        with open('energy') as fil:
            for line in fil:
                ll = line.split()
                if len(ll) == 4:
                    try:
                        ene = float(ll[1])
                    except:
                        pass
        result = ene
    except:
        # !!! failed
        result = None
    return result


def parse_cosmo():
    vol = None
    N = None
    try:
        with open('out.cosmo') as fil:
            for line in fil:
                if not vol and line.startswith('  volume='):
                    vol = float(line.split()[1])
                elif not N and line.startswith('  nps='):
                    N = int(line.split()[1])
                    coord = np.empty((N, 3))
                    area = np.empty(N)
                    potential = np.empty(N)
                elif line.startswith('$segment_information'):
                    for _ in range(10):
                        next(fil)
                    for i in range(N):
                        line = next(fil)
                        ll = [float(x) for x in line.split()]
                        coord[i] = ll[2:5]
                        area[i] = ll[6]
                        potential[i] = ll[8]
        result = vol, coord, area, potential
    except:
        # !!! failed
        result = (None, None, None, None)
    return result


def noop(*args, **kwargs):
    return None
