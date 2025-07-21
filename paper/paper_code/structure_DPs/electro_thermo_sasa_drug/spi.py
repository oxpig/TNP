import os
import numpy as np
from multiprocessing import Pool
import pandas as pd

cpus = 70

path = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_propka/'
files = [path + f for f in os.listdir(path) if f.endswith('.pka')] 
print(len(files))

def run(file):
    try:
        with open(file) as rf:
            pI_line = [line for line in rf.readlines() if line.startswith('The pI is')][0]
            return pI_line.split()[3], pI_line.split()[6]
    except Exception:
        return None, None


pool = Pool(processes = cpus)
results = pool.map(run, files)
pool.close()
pool.join()

print(len(results))

df = pd.DataFrame(results)
df.columns = ['folded_pI', 'unfolded_pI']
df.insert(0, 'ID', [s.split('/')[-1][:-4] for s in files])
df.to_csv('./ABB_paired_structural_pI_2.csv', index=None)
