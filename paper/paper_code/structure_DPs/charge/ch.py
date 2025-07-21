import os
import pandas as pd
import prody as prd
import numpy as np
from collections import defaultdict, Counter
from multiprocessing import Pool


def get_charge_dispersion(pdb):
    
    try:

        prot = prd.parsePDB(pdb, QUIET=True)
        positive_res = ['HIS', 'ARG', 'LYS']
        negative_res = ['ASP', 'GLU']

        central_atom = prd.pickCentralAtom(prot)
        prd.moveAtoms(prot, to=np.zeros(3))

        positive_distances = [prd.calcDistance(central_atom, prd.pickCentralAtom(res)) for res in prot.iterResidues() if res.getResname() in positive_res]
        negative_distances = [prd.calcDistance(central_atom, prd.pickCentralAtom(res)) for res in prot.iterResidues() if res.getResname() in negative_res]
        positive_dispersion = np.std(positive_distances)
        negative_dispersion = np.std(negative_distances)
        
    except Exception:
        positive_dispersion, negative_dispersion = 0, 0
    
    return positive_dispersion, negative_dispersion

cpus = 80

path = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_h/'
files = [path + f for f in os.listdir(path) if f.endswith('.pdb')]


pool = Pool(processes = cpus)
results = pool.map(get_charge_dispersion, files)
pool.close()
pool.join()
    
print(len(results))

df = pd.DataFrame(results)
df.columns = ['positive_charge_heterogeneity', 'negative_charge_heterogeneity']
df.insert(0, 'ID', [s.split('/')[-1][:-4] for s in files])
df.to_csv('./ABB_paired_charge_heterogeneity_2.csv', index=None)
