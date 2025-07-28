import numpy as np
import os
import json
import anarci
from Bio.PDB import PDBParser
from Bio import PDB
from ABDB.AbPDB.AntibodyParser import AntibodyParser
from ABDB.AbPDB.Select import select_all
#from ABDB.AB_Utils.sequence_liabilities import get_liabilities
from scripts.sequence_liabilities import get_liabilities # modified from sabdab version
from ABDB.AB_Utils.region_definitions import annotate_regions
from ImmuneBuilder.sequence_checks import number_sequences

# dict to convert three letter code to one letter code
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

available_schemes = ["imgt", "chothia", "kabat", "martin"]

# https://github.com/biopython/biopython/blob/master/Bio/PDB/DSSP.py
# Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
SASA_max = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLN": 225.0,
    "GLU": 223.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}

def get_sequences_from_pdb(input_pdb):
    # Run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', input_pdb)    

    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    model_sequences = []
    for model in structure:
        for chain in model:
            sequence = [d3to1[residue.resname] for residue in chain]
            model_sequences.append((chain.id, ''.join(sequence)))

    return model_sequences

def renumber_pdb(input_pdb, scheme="imgt"):
    """
    Renumbers an input structure into the IMGT numbering scheme.
    """
    schemes_list = available_schemes if scheme == "all" else [scheme]

    for scheme in schemes_list:
        # Define output filepaths
        output_dir = os.path.dirname(input_pdb)
        output_filename = ''.join([os.path.basename(input_pdb).split('.')[0],'_',scheme, '.pdb'])
        output_structure_file = os.path.join(output_dir, output_filename)

        # Renumber the prediction with the selected scheme
        p = AntibodyParser(QUIET=True)
        p.set_numbering_scheme(scheme=scheme)
        p.set_numbering_method("anarci")
        s = p.get_antibody_structure('modeller_out', input_pdb)

        for r in s.get_residues():
            r.parent.id = r.chain_type

        # Output structure file
        with open(output_structure_file, 'w') as out:
            out.write(s._get_output_string(select_all(), 1)[0])

def get_liabilities_from_pdb(input_pdb, scheme="imgt", save=False, restrict_species=True):
    """
    Get numberings for any or all available schemes
    """
    schemes_list = available_schemes if scheme == "all" else [scheme]
    
    liabilities = {}
    for scheme in schemes_list:
        # Define output filepaths
        output_dir = os.path.dirname(input_pdb)
        output_filename = ''.join([os.path.basename(input_pdb).split('.')[0],'_',scheme, '.pdb'])
        output_structure_file = os.path.join(output_dir, output_filename)

        # Run parser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', input_pdb)    

        # Iterate each model, chain, and residue
        # Extract chain sequences 

        model_sequences = []
        for model in structure:
            for chain in model:
                sequence = [d3to1[residue.resname] for residue in chain]
                model_sequences.append((chain.id, ''.join(sequence)))

        # Extract liabilities for extracted sequences
        # GG updated to deal with nanobodies only VHH
        try:
            heavy_sequence, light_sequence = model_sequences[0][1], model_sequences[1][1]
        except:
            heavy_sequence, light_sequence = model_sequences[0][1], ''
        liabilities[scheme], n_liabilities = get_liabilities(heavy_sequence,
                                                            light_sequence,
                                                            scheme=scheme,
                                                            save=save,
                                                            outfile= os.path.join(output_dir,"sequence_liabilities.csv"),
                                                            restrict_species=restrict_species)

    return liabilities, n_liabilities

def get_bfactors_from_pdb(input_pdb):
    """
    Extract B-factor values per residue for given input PDB
    """
    # Run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', input_pdb)    

    # Iterate each model, chain, and residue
    # Extract chain sequences 
    bfactors = {}
    for model in structure:
        for chain in model:
            bfactors[chain.id] = []
            for residue in chain:
                resnum, resname, bfactor = (residue.get_id(), d3to1[residue.get_resname()], list(residue)[0].get_bfactor())
                bfactors[chain.id].append((resnum, resname, bfactor))

    return bfactors

def get_range( residues, chain, scheme, definition):
    """
    Get the indices which correspond to CDRs for a particular definition.
    """
    regions = annotate_regions(residues, chain, numbering_scheme=scheme, definition=definition)       
    ranges = []
    framework, cdr = True, False
    for i in range( len(residues) ):
        reg = regions[i][2]
        if "cdr" in reg and framework:
            ranges.append(i)
            framework, cdr = False, True
        elif "fw" in reg and cdr:
            ranges.append(i-1)
            framework, cdr = True, False

    return ranges

def get_cdr_ranges_from_pdb(input_pdb, type = "", scheme = "imgt", definition="imgt"):
    """
    Get CDR index ranges corresponding to a particular defintion given an input PDB
    """
    sequences = dict(get_sequences_from_pdb(input_pdb))

    cdr_ranges = {}
    if type == "antibody":
        available_definitions = ["kabat", "chothia", "imgt", "north", "contact"]
        definitions_list = available_definitions if definition == "all" else [definition]
        
        numberings = number_sequences(sequences, scheme=scheme)
        heavy_numbering = numberings['H']
        light_numbering = numberings['L']
        for definition in definitions_list:
            cdr_ranges[f"{definition}_CDRH_ranges"] = get_range(heavy_numbering, chain="H", scheme=scheme, definition=definition)
            cdr_ranges[f"{definition}_CDRL_ranges"] = get_range(light_numbering, chain="L", scheme=scheme, definition=definition)

    elif type == "nanobody":
        available_definitions = ["kabat", "chothia", "imgt", "north", "contact"]
        definitions_list = available_definitions if definition == "all" else [definition]
        
        numberings = number_sequences(sequences, scheme=scheme)
        heavy_numbering = numberings['H']
        for definition in definitions_list:
            cdr_ranges[f"{definition}_CDRH_ranges"] = get_range(heavy_numbering, chain="H", scheme=scheme, definition=definition)

    elif type == "tcr":
        numberings = number_sequences(sequences, scheme=scheme)
        beta_numbering = numberings['B']
        alpha_numbering = numberings['A']
        cdr_ranges[f"{definition}_CDRB_ranges"] = get_range(beta_numbering, chain="H", scheme=scheme, definition="imgt")
        cdr_ranges[f"{definition}_CDRA_ranges"] = get_range(alpha_numbering, chain="H", scheme=scheme, definition="imgt")

    else:
        print("ERROR: Not a valid model 'type'. Options are: 'antibody', 'nanobody', 'tcr'")

    return cdr_ranges

def exposed_buried(input_pdb, cutoff=0.075, probe_radius=1.4, n_points=100):
    """
    Calculate whether a resiude is exposed of buried
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('struct', input_pdb)    

    sr = PDB.SASA.ShrakeRupley(probe_radius=probe_radius, n_points=n_points)

    fab_structure = structure[0]

    sr.compute(fab_structure, level="R")

    rel_sasas = []
    for res in fab_structure.get_residues():
        rel_sasas.append((res.sasa / SASA_max[res.resname]) > cutoff)

    return rel_sasas

def run_dssp(structurefile, type=''):
    """
    Run dssp on the structurefile. This will calculate the solvent exposure and the phi/psi angles etc.
    Requires DSSP to be in the path.
    @param structurefile: The file to calculate the surface area with on
    @return: A dictionary containing H and L as keys and a list of residue identifiers and dictionaries containing the properties of each residue.
    """

    from Bio.PDB.DSSP import DSSP
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import three_to_one

    if type == 'tcr':
        dssp_properties = {"B":[], "A":[]}
        # Parse the model to get a structure object
        s = PDBParser(QUIET=True).get_structure("model", structurefile)

        # Run DSSP on the model to annotate residues
        try:
            _ = DSSP( s[0], structurefile) 
        except OSError: 
            # dssp is not installed so cannot be run.
            return {"B":[], "A":[]}, None, "DSSP is not in the path"

    else:
        dssp_properties = {"H":[], "L":[]}

        # Parse the model to get a structure object
        s = PDBParser(QUIET=True).get_structure("model", structurefile)

        # Run DSSP on the model to annotate residues
        try:
            _ = DSSP( s[0], structurefile) 
        except OSError: 
            # dssp is not installed so cannot be run.
            return {"H":[], "L":[]}, None, "DSSP is not in the path"

    # These are not required as Biopython kindly implements it for us :)
    # Standard solvent accessibilities for residues taken from the reference:
    # C. Chothia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14. 
#    standards={'CYS': 135, 'ASP': 150, 'SER': 115, 'ASN': 160, 'GLN': 180, 'LYS': 200,
#     'THR': 140, 'PRO': 145, 'HIS': 195, 'PHE': 210, 'ALA': 115, 'GLY': 75, 'ILE': 175,
#     'LEU': 170, 'ARG': 225, 'TRP': 255, 'VAL': 155, 'GLU': 190, 'TYR': 230, 'MET': 185}


    # Populate a dictionary of properties from dssp. These are stored in the xtra dictionary in the residue object.
    # They loop like:
#    {'PSI_DSSP': 142.5, 'EXP_DSSP_ASA': 88, 'EXP_DSSP_RASA': 0.44, 'SS_DSSP': 'E', 'PHI_DSSP': -148.40}

    outputfile = os.path.splitext( structurefile )[0]+"_annotations.csv"
    message=""
    with open(outputfile , 'w' ) as fout:
        print(",".join(["Chain","Index","Insertion","Amino Acid","Absolute Surface Area",
                                   "Relative Surface Area","Secondary Structure","Phi angle","Psi angle"]), file=fout)
        for chain in s[0]:
            for r in chain:
                if r.id[0] !="W":
                    dssp_properties[chain.id].append( [ r.id[1:], r.resname, r.xtra ] )
                    try:
                        print(",".join( [ chain.id, str(r.id[1]), r.id[2], three_to_one(r.get_resname()), "%.2f"%r.xtra['EXP_DSSP_ASA'], "%.2f"%r.xtra['EXP_DSSP_RASA'],
                                r.xtra['SS_DSSP'],"%.2f"%r.xtra['PHI_DSSP'],"%.2f"%r.xtra['PSI_DSSP'] ]), file=fout)
                    except KeyError:
                        print(",".join([chain.id, str(r.id[1]),r.id[2],"-","-","-","-","-","-"]), file=fout)
                        message = "DSSP annotations could not be calculated for at least one residue"

    return {"dssp_properties": dssp_properties, "annotation_file": outputfile, "dssp_message": message}

def get_bfactors_from_pdb_per_region(input_pdb, type= '', scheme='imgt', definition='imgt'):
    
    pdb_cdr_ranges = get_cdr_ranges_from_pdb(input_pdb, type=type, scheme=scheme, definition=definition)
    pdb_bfactors = get_bfactors_from_pdb(input_pdb)

    chain_regions_names = {}
    chain_regions_indices = {}
    for chain in pdb_bfactors:
        first_indx, last_indx = pdb_bfactors[chain][0][0][1], pdb_bfactors[chain][-1][0][1]
        chain_regions_names[chain] = [f'fw{chain.lower()}1',
                                      f'cdr{chain.lower()}1', 
                                      f'fw{chain.lower()}2',
                                      f'cdr{chain.lower()}2',
                                      f'fw{chain.lower()}3',
                                      f'cdr{chain.lower()}3',
                                      f'fw{chain.lower()}4']
        x = pdb_cdr_ranges[f'imgt_CDR{chain}_ranges']
        chain_regions_indices[chain] = list(zip([first_indx]+x, x+[last_indx]))

    bfactors_per_region = {}
    for chain in chain_regions_indices:
        for i in range(len(chain_regions_indices[chain])):
            bfactor_values = [x[-1] for x in pdb_bfactors[chain] if x[0][1] in range(*chain_regions_indices[chain][i])]
            prediction_score_mean = np.mean(bfactor_values)
            prediction_score_std = np.std(bfactor_values)
            bfactors_per_region[chain_regions_names[chain][i]] = [list(chain_regions_indices[chain][i]), prediction_score_mean, prediction_score_std]

    return bfactors_per_region

def get_modelling_details(input_pdb, type="", scheme="imgt", save=False):
    """
    Compile a modelling_details dict from an input PDB. This can be saved into JSON file
    """

    pdb_sequences = dict(get_sequences_from_pdb(input_pdb))
    pdb_numberings = number_sequences(pdb_sequences, scheme=scheme)
    pdb_bfactors = get_bfactors_from_pdb(input_pdb)

    modelling_details = {'scheme': scheme,
                         'sequences': pdb_sequences,
                         'prediction_scores': pdb_bfactors}

    if type == "antibody":
        pdb_liabilities, n_liabilities = get_liabilities_from_pdb(input_pdb, scheme=scheme, save=True, restrict_species=True)
        heavy_numbering = [''.join(map(str,x[0])) for x in pdb_numberings['H']]  
        light_numbering = [''.join(map(str,x[0])) for x in pdb_numberings['L']]
        pdb_cdr_ranges = get_cdr_ranges_from_pdb(input_pdb, type=type, scheme=scheme, definition='all')
        pdb_bfactors_per_region = get_bfactors_from_pdb_per_region(input_pdb, type='antibody')
        pdb_dssp_output = run_dssp(input_pdb)

        modelling_details.update({'heavy_numbering': heavy_numbering, 
                                  'light_numbering': light_numbering,
                                  'sequence_liabilities': pdb_liabilities[scheme],
                                  'n_possible_seq_liabilities': n_liabilities})

    elif type == "nanobody":
        # GG include liabilities for nanobodies 
        pdb_liabilities, n_liabilities = get_liabilities_from_pdb(input_pdb, scheme=scheme, save=True, restrict_species=True)
        heavy_numbering = [''.join(map(str,x[0])) for x in pdb_numberings['H']]
        pdb_cdr_ranges = get_cdr_ranges_from_pdb(input_pdb, type=type, scheme=scheme, definition='all')
        pdb_bfactors_per_region = get_bfactors_from_pdb_per_region(input_pdb, type='nanobody')
        pdb_dssp_output = run_dssp(input_pdb)

        # NOTE updated here to make sure liabilities included in json for nanobodies
        modelling_details.update({'heavy_numbering': heavy_numbering,
                                  'sequence_liabilities': pdb_liabilities[scheme],
                                  'n_possible_seq_liabilities': n_liabilities})

    elif type == "tcr":
        beta_numbering = [''.join(map(str,x[0])) for x in pdb_numberings['B']]
        alpha_numbering = [''.join(map(str,x[0])) for x in pdb_numberings['A']]
        pdb_cdr_ranges = get_cdr_ranges_from_pdb(input_pdb, type=type, scheme=scheme, definition='imgt')
        pdb_bfactors_per_region = get_bfactors_from_pdb_per_region(input_pdb, type='tcr')
        pdb_dssp_output = run_dssp(input_pdb, type='tcr')

        modelling_details.update({'beta_numbering': beta_numbering, 
                                  'alpha_numbering': alpha_numbering})

    else:
        print("ERROR: Not a valid model 'type'. Options are: 'antibody', 'nanobody', 'tcr'")
                              
    modelling_details.update(pdb_cdr_ranges)
    modelling_details.update(pdb_dssp_output)
    modelling_details.update({'prediction_scores_per_region': pdb_bfactors_per_region})

    if save:
        output_dir = os.path.dirname(input_pdb)

        # Save ANARCI numbering
        allowed_species = None
        anarci.anarci(list(pdb_sequences.items()),
                      scheme=scheme, output=True,
                      outfile=os.path.join(output_dir,"numbering.txt"),
                      allowed_species=allowed_species)

        # Define output filepaths
        output_json = os.path.join(output_dir, 'modelling_details.jsonp')

        with open(output_json,'w') as fp:
            json.dump(modelling_details, fp)

    return modelling_details