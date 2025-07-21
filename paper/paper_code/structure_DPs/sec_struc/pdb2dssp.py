import os

from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import ss_to_index

import functools
from collections import Counter
import pandas as pd
import numpy as np

from multiprocessing import Pool


cpus = 80

structure_path = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_h/'
structures = os.listdir(structure_path)

def dssp_to_dict(structure, path = structure_path):

    try:
        print(structure[:-4])
        p = PDBParser(QUIET=True)
        model = p.get_structure("", path + structure)[0]
        dssp = DSSP(model, path + structure)

        data = np.array(list(dict(dssp).values()))

        init_c = Counter({
            'H': 0,
            'B': 0,
            'E': 0,
            'G': 0,
            'I': 0,
            'T': 0,
            'S': 0,
            '-': 0
        })

        get_c = Counter(data[..., 2])

        ss_c = dict(functools.reduce(lambda a, b: a.update(b) or a, [init_c, get_c], Counter()))

        values = np.mean(data[..., 3:].astype(np.float32), axis=0)

        ss_d = {}
        keys = ['relative_ASA', 'phi', 'psi', 'NH_O_1_relidx', 
                'NH_O_1_energy', 'O_NH_1_relidx', 'O_NH_1_energy',
                'NH_O_2_relidx', 'NH_O_2_energy', 'O_NH_2_relidx', 
                'O_NH_2_energy']

        ss_d.update(ss_c)

        for n, key in enumerate(keys):
            ss_d.update({key: values[n]})
        
        return ss_d

    except Exception:
        print('Error:', structure[:-4])
        ss_d = {
            'H': 0,
            'B': 0,
            'E': 0,
            'G': 0,
            'I': 0,
            'T': 0,
            'S': 0,
            '-': 0,
            'relative_ASA': 0, 
            'phi': 0, 
            'psi': 0, 
            'NH_O_1_relidx': 0, 
            'NH_O_1_energy': 0,
            'O_NH_1_relidx': 0, 
            'O_NH_1_energy': 0,
            'NH_O_2_relidx': 0,
            'NH_O_2_energy': 0, 
            'O_NH_2_relidx': 0, 
            'O_NH_2_energy': 0
        }
        return ss_d


       
pool = Pool(processes = cpus)
results = pool.map(dssp_to_dict, structures)
pool.close()
pool.join()

print(len(results))
#new_results = []
#for result in results:
#    _ = {k: {k2: 0 if v2 == 'NaN' else v2 for k2, v2 in v.items()} for k, v in result.items()}

#new_results = []
#for i in results:
#    if i is not None:
#        new_results.append(i)
#    elif i is None:
#        new_results.append({H,B,E,G,I,T,S,-,relative_ASA,phi,psi,NH_O_1_relidx,NH_O_1_energy,O_NH_1_relidx,O_NH_1_energy,NH_O_2_relidx,NH_O_2_energy,O_NH_2_relidx,O_NH_2_energy})

df = pd.DataFrame(results)
#print(df.shape)

df.insert(0, 'ID', [s[:-4] for s in structures])
df.to_csv('./ABB_paired_dssp_2.csv', index=None)
#print(df.shape)

