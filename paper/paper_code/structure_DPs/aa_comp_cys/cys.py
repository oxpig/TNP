import os
import prody as prd
import math
from collections import defaultdict
import itertools
import pandas as pd
from multiprocessing import Pool

def distance(distances_1, distances_2):
    x1, y1, z1 = distances_1
    x2, y2, z2 = distances_2
    distance = math.sqrt(math.pow(x2 - x1, 2) +
               math.pow(y2 - y1, 2) +
               math.pow(z2 - z1, 2) * 1.0)
    return distance

def ss_bonds(pdb):
    
    try:
        cysteines = defaultdict(list)
        free_cys = 0 
        cys_bridges = 0

        prot = prd.parsePDB(pdb, QUIET=True)
        for residue in prot.iterResidues():
            name_id = str(residue).split()
            name_id. insert(0, str(residue.getChain()).split()[1])
            if name_id[1] == 'CYS':
                cysteines['_'.join(name_id)].extend(list(residue.getAtom('SG').getCoords()))
                #print('_'.join(name_id))

        ss_bonds = []

        for pair in itertools.combinations(cysteines.keys(), r=2):
            if pair[0] != pair[1]:
                d = distance(cysteines[pair[0]], cysteines[pair[1]])
                if d <= 4.0:
                    ss_bonds.append((pair) + (d, ))
                    cys_bridges += 1
                else:
                    free_cys +=1
                    
    except Exception:
        free_cys, cys_bridges = 0, 0

    return free_cys, cys_bridges


cpus = 70

path = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_h/'
files = [path + f for f in os.listdir(path) if f.endswith('.pdb')]

pool = Pool(processes = cpus)
results = pool.map(ss_bonds, files)
pool.close()
pool.join()
    
print(len(results))

df = pd.DataFrame(results)
df.columns = ['free_cys', 'cys_bridges']
df.insert(0, 'ID', [s.split('/')[-1][:-4] for s in files])
df.to_csv('./ABB_paired_CYS_2.csv', index=None)
