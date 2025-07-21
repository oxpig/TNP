import os
import json
from operator import itemgetter
from collections import defaultdict, Counter
import numpy as np
from multiprocessing import Pool
import pandas as pd
import freesasa
import copy

cpus = 70

path_1 = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_propka/'
path_2 = '/storage/evagsm/nobackup/crystal_dataset/structures/ABB_paired_h/'

files = [f for f in os.listdir(path_1) if f.endswith('.pka')] 

files_1 = [path_1 + f for f in files]
files_2 = [path_2 + f.replace('pka', 'pdb') for f in files]

#files_1 = [path_1 + f for f in os.listdir(path_1) if f.endswith('.pka')]
#files_2 = [path_2 + f for f in os.listdir(path_2) if f.endswith('.pdb')]

print(len(files_1), len(files_2))
print(files_1[:3])
print(files_2[:3])


# Black, S. D., & Mould, D. R. (1991). Development of hydrophobicity parameters to analyze proteins which bear post- or cotranslational modifications. Analytical Biochemistry, 193(1), 72–82. doi:10.1016/0003-2697(91)90045-u 
hydrophobicity = {'ALA': 0.616,
                  'CYS': 0.680,
                  'ASP': 0.028,
                  'GLU': 0.043,
                  'PHE': 1.000,
                  'GLY': 0.501,
                  'HIS': 0.165,
                  'ILE': 0.943,
                  'LYS': 0.283,
                  'LEU': 0.943,
                  'MET': 0.738,
                  'ASN': 0.236,
                  'PRO': 0.711,
                  'GLN': 0.251,
                  'ARG': 0.000,
                  'SER': 0.359,
                  'THR': 0.450,
                  'VAL': 0.825,
                  'TRP': 0.878,
                  'TYR': 0.880}

def propka_data(f):
     
    #files = [f for f in os.listdir(path) if f.endswith('.pka')]
    
    try:
        
        ph_dicts = defaultdict(dict)
        #res_dicts = defaultdict(dict)
    
        with open(f) as rf:
            all_lines = [l.strip() for l in rf.readlines()] 

        sp = [n for n, l in enumerate(all_lines) if l.startswith('SUMMARY OF THIS PREDICTION')]
        ffe = [n for n, l in enumerate(all_lines) if l.startswith('Free energy of')]
        pc = [n for n, l in enumerate(all_lines) if l.startswith('Protein charge of folded')]

        summary_prediction = [i for i in all_lines[sp[0]+2:ffe[0]] if not i.startswith('-')]
        folding_free_energy = [i for i in all_lines[ffe[0]+1:pc[0]] if i and i[0].isnumeric()]
        protein_charge = [i for i in all_lines[pc[0]+2:-1]]

        ph_dict = defaultdict(list)   # key: pH; values: energy of folding (kcal/mol), protein charge of folded state, protein charge of unfolded state
        #res_dict = defaultdict(list)   # key: residue name_residue id; values: chain, pKa, model pKa

        for i, j in zip(folding_free_energy, protein_charge):
            new_i, new_j = i.split(), j.split()
            if new_i[0] == new_j[0]:
                ph_dict[float(new_i[0])].extend([float(new_i[1]), float(new_j[1]), float(new_j[2])])
        
        #for i in summary_prediction:
        #    new_i = i.split()
        #    key = '_'.join(new_i[:3])
        #    res_dict[key].extend([float(new_i[3]), float(new_i[4])])
        #
        name = f.split('/')[-1][:-4]
        ph_dicts[name].update(ph_dict)
        #res_dicts[name].update(res_dict)
        
    except Exception:
        print('PROPKA:', f.split('/')[-1][:-4]) 
        ph_dict = {0.0: [np.inf, np.inf, np.inf], 
                   1.0: [np.inf, np.inf, np.inf], 
                   2.0: [np.inf, np.inf, np.inf], 
                   3.0: [np.inf, np.inf, np.inf], 
                   4.0: [np.inf, np.inf, np.inf], 
                   5.0: [np.inf, np.inf, np.inf], 
                   6.0: [np.inf, np.inf, np.inf], 
                   7.0: [np.inf, np.inf, np.inf], 
                   8.0: [np.inf, np.inf, np.inf], 
                   9.0: [np.inf, np.inf, np.inf], 
                   10.0: [np.inf, np.inf, np.inf], 
                   11.0: [np.inf, np.inf, np.inf], 
                   12.0: [np.inf, np.inf, np.inf], 
                   13.0: [np.inf, np.inf, np.inf], 
                   14.0: [np.inf, np.inf, np.inf]}
        
    return f.split('/')[-1][:-4], ph_dict
    
    #return ph_dicts, res_dicts

def exposed_sasas():
    
    backbone = ['C', 'N', 'CA', 'O']
    
    exposed_sasas = defaultdict(dict)
    
    tripeptides = ['A_A_A.pdb', 'A_C_A.pdb', 'A_D_A.pdb', 'A_E_A.pdb', 
                   'A_F_A.pdb', 'A_G_A.pdb', 'A_H_A.pdb', 'A_I_A.pdb',
                   'A_K_A.pdb', 'A_L_A.pdb', 'A_M_A.pdb', 'A_N_A.pdb',
                   'A_P_A.pdb', 'A_Q_A.pdb', 'A_R_A.pdb', 'A_S_A.pdb',
                   'A_T_A.pdb', 'A_V_A.pdb', 'A_W_A.pdb', 'A_Y_A.pdb']
    
    for peptide in tripeptides:
        structure = freesasa.Structure(peptide)
        sasa = freesasa.calc(structure)
        
        for atom in range(structure.nAtoms()):
            if structure.residueNumber(atom).strip() == '3':
                data = structure.residueName(atom).strip(), structure.atomName(atom).strip(), sasa.atomArea(atom)
                exposed_sasas[data[0]].update({data[1]: data[2]})
                
    exposed_side_chain_sasas = {k: 0 for k in exposed_sasas.keys()}
    
    for aa in exposed_sasas:
        for atom in exposed_sasas[aa].keys():
            if not atom in backbone:
                exposed_side_chain_sasas[aa] += exposed_sasas[aa][atom]

    return exposed_side_chain_sasas


exposed_side_chain_sasas = exposed_sasas()


def sap(filename, exposed_side_chain_sasas = exposed_side_chain_sasas, hydrophobicity = hydrophobicity):
    
    def sasa(filename):

        backbone = ['C', 'N', 'CA', 'O']

        sasas = defaultdict(dict)
        sasa_int = 0

        structure = freesasa.Structure(filename)
        sasa = freesasa.calc(structure)

        for atom in range(structure.nAtoms()):
            data = structure.residueName(atom).strip(), structure.residueNumber(atom).strip(), structure.chainLabel(atom).strip(), structure.atomName(atom).strip(), sasa.atomArea(atom)
            idx = '_'.join(data[:3])
            sasas[idx].update({data[3]: data[4]})
            sasa_int += data[4]

        return sasas, sasa_int
    
    def side_chain_sasa(sasas):

        backbone = ['C', 'N', 'CA', 'O']

        side_chain_sasas = {k: 0 for k in sasas.keys()}
        side_chain_sasa_int = 0

        for res in sasas:
            for atom in sasas[res].keys():
                if not atom in backbone:
                    side_chain_sasas[res] += sasas[res][atom]
                    side_chain_sasa_int += sasas[res][atom]

        return side_chain_sasas, side_chain_sasa_int

    try:
    
        sasas = sasa(filename)
        side_chain_sasas_ = side_chain_sasa(sasas[0])
        side_chain_sasas = side_chain_sasas_[0]

        sap_values = {}

        for i in side_chain_sasas:
            for j in exposed_side_chain_sasas:
                if i.startswith(j):
                    try:
                        #print(i, j, side_chain_sasas[i] / exposed_side_chain_sasas[j])
                        sasa = side_chain_sasas[i] / exposed_side_chain_sasas[j]
                    except ZeroDivisionError:
                        if i.startswith('GLY'):
                            #print(i, j, 0)
                            sasa = 0 
                    sap = sasa * hydrophobicity[j]
                    sap_values[i] = sap

        sap_score = sum(sap_values.values())
        return filename.split('/')[-1][:-4], sasas[1], side_chain_sasas_[1], sap_values, sap_score
        
    except Exception:
        print('SASA:', filename.split('/')[-1][:-4])
        sap_values = {'SER_2_L': np.inf, 'VAL_3_L': np.inf, 'LEU_4_L': np.inf}
        
        return filename.split('/')[-1][:-4], np.inf, np.inf, sap_values, np.inf    

    
pool = Pool(processes = cpus)

propka_dicts = np.array(pool.map(propka_data, files_1))
saps = np.array(pool.map(sap, files_2))

names_1 = propka_dicts[..., 0]
names_2 = saps[..., 0]
#names = np.intersect1d(names_1, names_2)

print(len(names_1), len(names_2))

saps_ = np.array(saps[..., -1])
propka = np.array([propka[7.0] for propka in propka_dicts[..., 1]])

#def DI(charges, sap_score, beta = 0.0815):
#    
#    di = {}
#    for k, v in charges.items():
#        di[k] = sap_score - beta * v[1]
#        
#    return di

beta = 0.0815

dis = saps_ - beta * propka[..., 2]

pool.close()
pool.join()

print(len(propka_dicts))

df = pd.DataFrame(propka)
df.columns = ['folding_free_energy', 'unfolded_charge', 'folded_charge']
df[['sasas', 'side_chain_sasas']] = saps[..., 1:3]
df['SAP'] = saps_
df['DI'] = dis
df.insert(0, 'ID', [s.split('/')[-1][:-4] for s in files_1])
df.replace([np.inf, -np.inf], 0, inplace=True)
df.to_csv('./ABB_paired_SASA_SAP_ID_2.csv', index=None)
