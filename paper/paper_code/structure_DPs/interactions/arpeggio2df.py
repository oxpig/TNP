import os
import json
from operator import itemgetter
from collections import defaultdict, Counter
import numpy as np
from multiprocessing import Pool
import pandas as pd

cpus = 80

path = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_arpeggio/'

files = [path + f for f in os.listdir(path) if f.endswith('.json')]

def flatten(t):
    return [item for sublist in t for item in sublist]


def run(f):
    
    print(f)
    d = {}
    with open(f) as rf:
        json_d = json.load(rf)

        contacts = Counter(flatten([j['contact'] for j in json_d]))

        interactions = [itemgetter('distance', 'interacting_entities', 'type')(j) for j in json_d]
        interacting_entities = Counter([i[1] for i in interactions])
        types = Counter([i[2] for i in interactions])

        d.update(contacts)
        d.update(interacting_entities)
        d.update(types)
        d.update({'distance': np.mean([i[0] for i in interactions])})
        
    return d


pool = Pool(processes = cpus)
results = pool.map(run, files)
pool.close()
pool.join()

print(len(results))

d_protint = defaultdict(dict)

for n, (f, r) in enumerate(zip(files, results)):
    #print(n)
    name = f.split('/')[-1].split('.')[0]
    d_protint[name] = r
    
df = pd.DataFrame(d_protint.values(), index=d_protint.keys())

df.to_csv('./ABB_paired_arpeggio_int_2.csv') #, index=None)
