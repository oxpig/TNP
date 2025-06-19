# -*- coding: utf-8 -*-
import numpy as np
from os import listdir
import os,sys,tempfile
from os.path import isfile, join
from .Common.PDBUtils import PDBchain
import pickle

local_path = os.path.dirname(os.path.realpath(__file__))

##################################################################################
#FUNCTION to REMOVE INSERTION CODES from NUMBERING
#
def stripInsertionCode(res):
        insertion = ""
        icode = ['Z','X','C','V','B','N','M','L','K','J','H','G','F','D','S','A','W','Q','E','R','T','Y','U','I','O','P'] #eg. 100A -> 100, 100B -> 100...
        for q in icode:
            if res.endswith(q):        #if last letter is non-numical
                #print res
                res,insertion = res[0:(len(res)-1)],q    #strip it
                #print res,insertion
        return res,insertion

##################################################################################
#FUNCTION to CHECK if a RESIDUE is a CDR RESIDUE (+/- 2) in Chothia/IMGT

def is_CDR(res,deff):
      
        _res,ins = stripInsertionCode(res)
        for CDR in definitions[deff]:
                if _res in definitions[deff][CDR]:
                        return CDR
        return False

##################################################################################
#FUNCTION to CHECK if a RESIDUE is a SPECIFIC CDR RESIDUE (+/-2) in Chothia/IMGT

def is_partic_CDR(res,deff,cdr_id):
        
        _res,ins = stripInsertionCode(res)
        if _res in definitions[deff][str(cdr_id)]:
                return cdr_id
        return False

##################################################################################
#PREREQUISITES#

#Classifying residues by hydrophobicity
Naive_HH = {    
'ALA': 'B',
'ARG': 'C',
'ASN': 'U',
'ASP': 'C',
'CYS': 'S',
'GLU': 'C',
'GLN': 'U',
'GLY': 'S',
'HIS': 'C',
'ILE': 'H',
'LEU': 'H',
'LYS': 'C',
'MET': 'B',
'PHE': 'H',
'PRO': 'S',
'SER': 'U',
'THR': 'U',
'TRP': 'H',
'TYR': 'H',
'VAL': 'H'}

coloring = {'H':100,'C':25,'U':50,'S':75}

#Three letter:One letter translation table
AAtable = {
'ALA': 'A',
'ARG': 'R',
'ASN': 'N',
'ASP': 'D',
'CYS': 'C',
'GLU': 'E',
'GLN': 'Q',
'GLY': 'G',
'HIS': 'H',
'ILE': 'I',
'LEU': 'L',
'LYS': 'K',
'MET': 'M',
'PHE': 'F',
'PRO': 'P',
'SER': 'S',
'THR': 'T',
'TRP': 'W',
'TYR': 'Y',
'VAL': 'V'}

#CDR Definitions - Chothia/IMGT +2 on either side...


definitions = {                 
                "chothia" : {"H1" : ["H24","H25","H26", "H27", "H28", "H29", "H30", "H31", "H32","H33","H34"], 
                                        "H2" : ["H50","H51","H52", "H53", "H54", "H55", "H56","H57","H58"] ,
                                        "H3" : ["H93","H94","H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102","H103","H104"]},
                "imgt" : {"H1" : ["H25","H26","H27","H28","H29","H30","H31","H32","H33","H34","H35","H36","H37","H38","H39","H40"],
                                        "H2" : ["H54","H55","H56","H57","H58","H59","H60","H61","H62","H63","H64","H65","H66", "H67"],
                                        "H3" : ["H103","H104","H105","H106","H107","H108","H109","H110","H111","H112","H113","H114","H115","H116","H117","H118","H119"]
                        }
                }

##############################################################################################
#FETCHES  the TOTAL AREA of the SIDE-CHAIN

def parsePSAArea(psafile):

        dicto = dict()
        curr_chain = 'H'
        dicto['H'] = dict()
        dicto['L'] = dict()     #Initialise an empty dictionary, dicto = {'H':{},'L':{}}
        last = -1

        for line in open(psafile):
                        line =        line.strip()    #Creat a list, items of which were separated by \n
                        
                        if line[0:6] =='ACCESS':
                                sid_int = int(line[6:11])
                                if sid_int<last:
                                        curr_chain = 'L'
                                last = sid_int
                                sid = line[6:12].replace(" ","")
                                dicto[curr_chain][sid] = float(line[55:61])
                                
        return dicto

#############################################################################################
#PARSES the PSA OUTPUT, Chains should be HEAVY FIRST, then LIGHT.

def parsePSA(psafile, verbose=True):

        dicto = dict()
        curr_chain = 'H'
        dicto['H'] = dict()
        dicto['L'] = dict()     #Initialise an empty dictionary, dicto = {'H':{},'L':{}}
        last = -1
        
        for line in open(psafile):
                        line =        line.strip()    #Create a list, items of which were separated by \n
                        # if verbose:
                        #     print(line)
                        if line[0:6] =='ACCESS':
                                sid_int = int(line[6:11])
                                if sid_int<last:
                                        curr_chain = 'L'
                                last = sid_int
                                sid = line[6:12].replace(" ","")
                                dicto[curr_chain][sid] = float(line[61:67])     #This code populates the dictionary with the appropriate info...
        return dicto

#############################################################################################
#STARTS the PSA SOFTWARE - VERSION DEPENDS ON OS

def runPSA(pdb_file,where_to,verbose=True):
        
        extension = ""

        #If we are on a mac, run the mac version.
        if (sys.platform=='darwin'):
                extension = "_mac"
        os.system("psa"+extension+" -t "+pdb_file+" > "+where_to)      #the bash code to run PSA

#############################################################################################
#CREATES a TEMPORARY FOLDER in the USER DIRECTORY. Returns as a variable#
def create_temp_folder():
        dirpath = tempfile.mkdtemp()
        return dirpath

#############################################################################################
#NORMALISES HYDROPHOBICITY to a Z-SCORE#

def normalize(matrix,index):

        vals = []

        for residue in matrix:
                vals.append(matrix[residue][index]) #adds the appropriate hydrophobicity value H(R,S), unnormalised hydrophobicity R in scoring scheme S
        
        annotation = dict()
        
        mu = np.mean(vals)      #compute mean of the list
        sigma = np.std(vals)    #compute sd of the list
        
        for residue in matrix:

                #Normalizing to between 1 and 2
                annotation[residue]= ((matrix[residue][index]+(0-min(vals)))/((max(vals)-min(vals))))+1.0
        
        return annotation

############################################################################################
#FUNCTION to EVALUATE CHARGE of RESIDUES#


#######################################################################################
#MORE PREREQUISITES

#Which of the three classes do we belong to: positive, negative, or hydrophobic (or none)?
residue_classes = ['hydrophobic','positive','negative']

#######################################################################################
#PRIMARY ANNOTATING FUNCTION#

def CreateAnnotation(annotation_index,which_ph,input_file,chains,numbering_scheme,verbose=True):
        
        #Load In All Hydrophobicity Information from hydrophobics.txt

        hydrophobics = dict()
        
        for line in open(join(local_path, 'dat', 'Hydrophobics.txt')):
                line = line.strip().split('\t')
                hydrophobics[line[0].upper()] = (float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]))    #key = ALA,ILE etc. value = tuple with 5 separate scores
        
        #Normalize the hydrophobics

        hydrophobicity_annotation = normalize(hydrophobics,annotation_index)    #calls up a function to normalise the scores within requested scoring system
        
        #Set up a temporary directory

        temp_dir = create_temp_folder()
        
        os.mkdir(join(temp_dir,'surface'))      #Creates a temporary directory called surface

        if verbose:
            print("Temporary results stored in", temp_dir)

        chothia_file = join(temp_dir,'input.pdb')

        os.system('cp '+input_file+' '+chothia_file)    #note variable is called chothia_file, but may be numbered in imgt

        
        #STEP 1. Get the surface-exposed residues by running PSA

        surface_file = join(temp_dir,'surface','result.psa')    #generate a new temporary file in temp directory 'surface'

        runPSA(chothia_file,surface_file)  #Starts up the PSA software, depending on your current OS

        asa = parsePSA(surface_file, verbose=verbose)    #Returns a dictionary with PSA values for each H chain
        
        asa_area = parsePSAArea(surface_file)   #Returns a dictionary with total areas for each side chain within each H chain
        

        #STEP 2: Evaluate Hydrophobicity-Scaled Accessible Surface Areas (HASAs) and Charges
        
        nb_structure = {'H':PDBchain(chothia_file,'H')}  #Load in the ab heavy and light chains in following format
        heavy_fv_charge = 0
        
        h1_charge = 0
        h2_charge = 0
        h3_charge = 0

        num_cdr_charged = 0
        num_heavy_fv_charged = 0

        residues = dict()
        
        for chain in ['H']:
                for elem in sorted(nb_structure[chain].residues):
                        
                        typ = nb_structure[chain].residues[elem].res_type   #Assign the residue a type

                        r_elem = {'res':nb_structure[chain].residues[elem],'typ':typ,'charge':0}   #store this type, add future properties
                        bfact = 0
                        res_asa = 0
                        res_asa_area = 0

                        try:
                                res_asa = asa[chain][str(elem[0])+elem[1]]      #Looks up in the returned dictionaries from PSA software for asa['H']['4H'], for example
                                res_asa_area = asa_area[chain][str(elem[0])+elem[1]]    #Looks up as above, but in the area dictionary
                        
                        except KeyError:
                                #IN PROGRESS: Sometimes ASA skips residues -- especially the terminal ones.
                                sys.stderr.write("ASA ERROR: %s" % input_file)  #Write these messages to standard error for diagnosis
                                
                        if res_asa < 7.5:       #If residue is buried (accessible surface area of less than 7.5, ignore.
                                continue

                        r_elem['asa'] = res_asa             #Populates residue dictionary with asa score
                        r_elem['asa_sum'] = res_asa_area    #Populates residue diction with asa area score
                        if np.isnan(r_elem['asa_sum']) == True:
                            if verbose: print("NAN Residue Detected: "+r_elem['typ'],chain,elem[0]+elem[1])
                            break
                        r_elem['hydrophobic'] = float(hydrophobicity_annotation[typ]) #Looks up normalised hydrophobicity score for this residue, populates residue dictionary.        
                
                        
                        residues[(chain,elem)] = r_elem     #Finally collates all residues in the ab together, key is chain and number and value is their property dictionary
                

                #Define CDR neighborhood- so as not to calculate CDRs only...
                cdr_vicinity = []

                #Salt Bridge List, for second iteration...
                sb_list = []
                num_saltbridges = 0

                #Calculate cdrs and their neighbors.
                number_cdr_residues = 0
                number_heavy_fv_residues = 0       
        
                for r1 in residues:
                        numbering = r1[0]+str(r1[1][0])          #residues dictionary has tuple keys: (chain,elem); so eg. H111, no insertion code
                        if is_CDR(numbering,numbering_scheme)==False:
                                if r1[0] == "H":

                                        if AAtable[str(residues[r1]['typ'])] == "D" or AAtable[str(residues[r1]['typ'])] == "E":
                                                heavy_fv_charge = heavy_fv_charge - 1
                                                num_heavy_fv_charged += 1
                                                residues[r1]['negative'] = -1
                                                residues[r1]['charge'] = -1
                                        elif AAtable[str(residues[r1]['typ'])] == "R" or AAtable[str(residues[r1]['typ'])] == "K":
                                                heavy_fv_charge = heavy_fv_charge + 1
                                                num_heavy_fv_charged += 1
                                                residues[r1]['positive'] = 1
                                                residues[r1]['charge'] = 1
                                        elif AAtable[str(residues[r1]['typ'])] == "H":
                                                heavy_fv_charge = heavy_fv_charge + 0.1        #appropriate at pH 7.4 according to Henderson-Hasselbalch. Could generalise method to any given pH.
                                                num_heavy_fv_charged += 1
                                                residues[r1]['positive'] = 0.1
                                                residues[r1]['charge'] = 0.1
                                
                                        number_heavy_fv_residues += 1

                                        continue

                        else:
                                if r1[0] == "H":
                                        if is_partic_CDR(numbering,numbering_scheme,"H1")=="H1":
                                                if AAtable[str(residues[r1]['typ'])] == "D" or AAtable[str(residues[r1]['typ'])] == "E":
                                                        h1_charge = h1_charge - 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['negative'] = -1
                                                        residues[r1]['charge'] = -1
                                                elif AAtable[str(residues[r1]['typ'])] == "R" or AAtable[str(residues[r1]['typ'])] == "K":
                                                        h1_charge = h1_charge + 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 1
                                                        residues[r1]['charge'] = 1
                                                elif AAtable[str(residues[r1]['typ'])] == "H":
                                                        h1_charge = h1_charge + 0.1        #appropriate at pH 7.4 according to Henderson-Hasselbalch. Could generalise method to any given pH by evaluating this equation.
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 0.1
                                                        residues[r1]['charge'] = 0.1

                                        if is_partic_CDR(numbering,numbering_scheme,"H2")=="H2":
                                                if AAtable[str(residues[r1]['typ'])] == "D" or AAtable[str(residues[r1]['typ'])] == "E":
                                                        h2_charge = h2_charge - 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['negative'] = -1
                                                        residues[r1]['charge'] = -1
                                                elif AAtable[str(residues[r1]['typ'])] == "R" or AAtable[str(residues[r1]['typ'])] == "K":
                                                        h2_charge = h2_charge + 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 1
                                                        residues[r1]['charge'] = 1
                                                elif AAtable[str(residues[r1]['typ'])] == "H":
                                                        h2_charge = h2_charge + 0.1        #appropriate at pH 7.4 according to Henderson-Hasselbalch. Could generalise method to any given pH by evaluating this equation.
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 0.1
                                                        residues[r1]['charge'] = 0.1

                                        if is_partic_CDR(numbering,numbering_scheme,"H3")=="H3":
                                                if AAtable[str(residues[r1]['typ'])] == "D" or AAtable[str(residues[r1]['typ'])] == "E":
                                                        h3_charge = h3_charge - 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['negative'] = -1
                                                        residues[r1]['charge'] = -1
                                                elif AAtable[str(residues[r1]['typ'])] == "R" or AAtable[str(residues[r1]['typ'])] == "K":
                                                        h3_charge = h3_charge + 1
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 1
                                                        residues[r1]['charge'] = 1
                                                elif AAtable[str(residues[r1]['typ'])] == "H":
                                                        h3_charge = h3_charge + 0.1        #appropriate at pH 7.4 according to Henderson-Hasselbalch. Could generalise method to any given pH by evaluating this equation.
                                                        num_cdr_charged += 1
                                                        residues[r1]['positive'] = 0.1
                                                        residues[r1]['charge'] = 0.1
                                        

                                        if AAtable[str(residues[r1]['typ'])] == "D" or AAtable[str(residues[r1]['typ'])] == "E":
                                                heavy_fv_charge = heavy_fv_charge - 1
                                                num_heavy_fv_charged += 1
                                                residues[r1]['negative'] = -1
                                                residues[r1]['charge'] = -1
                                        elif AAtable[str(residues[r1]['typ'])] == "R" or AAtable[str(residues[r1]['typ'])] == "K":
                                                heavy_fv_charge = heavy_fv_charge + 1
                                                num_heavy_fv_charged += 1
                                                residues[r1]['positive'] = 1
                                                residues[r1]['charge'] = 1
                                        elif AAtable[str(residues[r1]['typ'])] == "H":
                                                heavy_fv_charge = heavy_fv_charge + 0.1        #appropriate at pH 7.4 according to Henderson-Hasselbalch. Could generalise method to any given pH.
                                                num_heavy_fv_charged += 1
                                                residues[r1]['positive'] = 0.1
                                                residues[r1]['charge'] = 0.1
                                        number_heavy_fv_residues += 1


                                number_cdr_residues += 1

                                if r1 not in cdr_vicinity:
                                        cdr_vicinity.append(r1)     #build up cdr vicinity list, populated with resiudes that comprise the CDR +/- 2

                        for r2 in residues:

                                #Skip if we are the same residue or if the other is a CDR as well!
                                if r1 == r2 or r2 in cdr_vicinity:
                                        continue
                        
                                #If we are here, we know r1 is CDR and r2 is not. So let's measure whether r2 is also in the neighborhood...
                                res1 = residues[r1]['res']
                                res2 = residues[r2]['res']
                                d = res2.distance(res1,res2)    #Works out the lowest pairwise distance between all atoms in both residues

                                if d < 4.0:       #Arbitrary cutoff distance, here 4.0A, to assign non-CDR residue as being in the "CDR neighbourhood".
                                        if r2 not in cdr_vicinity:
                                                cdr_vicinity.append(r2)     #Add this one if not already in the list!
                                                
                for r1 in residues:
                        for r2 in residues:
                                if r1 == r2:
                                        continue
                        
                        res1 = residues[r1]['res']
                        res2 = residues[r2]['res']
                        
                        donors = ["LYS","ARG"]      #not using histidine as only infrequently positively charged at pH 7.4
                        acceptors = ["ASP","GLU"]

                        if res1.res_type in donors:
                                if res2.res_type in acceptors:
                                        d = res1.sb_distance(res1,res2)
                                        if d < 3.2:     #Literature recommended cutoff for salt bridge formation
                                                residues[r1]['Positive'] = 0
                                                residues[r1]['Charge'] = 0
                                                residues[r1]['hydrophobic'] = hydrophobicity_annotation["GLY"]
                                                residues[r2]['Negative'] = 0
                                                residues[r2]['Charge'] = 0      #Neutralise all charges - assume in permanent salt bridge so not physiologically relevant
                                                residues[r2]['hydrophobic'] = hydrophobicity_annotation["GLY"]
                                                sb_list.append((str(r1[0])+str(r1[1][0])+str(r1[1][1])+": "+str(res1.res_type),str(r2[0])+str(r2[1][0])+str(r2[1][1])+": "+str(res2.res_type),str("{0:.4f}".format(d))))
                                                num_saltbridges += 1

                #ADJACENCY MATRIX FOR JUST THE CDRs
                
                adj_matrix_cdr = dict()
                for i in residue_classes:
                        adj_matrix_cdr[i] = dict()
                        for j in residue_classes:
                                adj_matrix_cdr[i][j] = 0

                for r1 in cdr_vicinity:
                        for r2 in cdr_vicinity:
                                if r1 != r2:        # just for residues in the CDR vicinity

                                        res1 = residues[r1]['res']
                                        res2 = residues[r2]['res']
                                        d = res2.distance(res1,res2)

                                        # Adjacency matrix looks only at values within 7.5 A from each other.
                                        
                                        if d > 7.5:
                                                continue        #ignore pairwise combination if >7.5A apart

                                        for c1 in residue_classes:        
                                                for c2 in residue_classes:
                                                    if c1 == c2:        #Only interested in patches of the same type of interaction (eg. hydrophob, +ve charge or -ve charge)
                                                        if c1 in residues[r1] and c2 in residues[r2]:
                                                                #The Coulombic measure.
                                                                coulombic_score = (residues[r1][c1]*residues[r2][c2])/(d*d)
                                                                adj_matrix_cdr[c1][c2] +=coulombic_score    #THIS IS THE CDR PATCH HYDROPHOBICITY SCORE!
                                                                

                ##########################################################################################################

                #Full statistics used for calculations.                

                stats = dict()
                stats[annotation_index] = dict()
                stats[annotation_index]['Patch_Hydrophob_CDR'] = float("{0:.4f}".format(adj_matrix_cdr['hydrophobic']['hydrophobic']))     #CDR Vicinity Patch HASA
                stats[annotation_index]['Patch_Pos_Charge_CDR'] = float("{0:.4f}".format(adj_matrix_cdr['positive']['positive']))     #CDR Vicinity Patch +CASA
                stats[annotation_index]['Patch_Neg_Charge_CDR'] = float("{0:.4f}".format(adj_matrix_cdr['negative']['negative']))     #CDR Vicinity Patch -CASA

                os.system('rm -rf '+temp_dir)   #Finally, forcibly remove all temporary files + folders

        #WRITE THE ANNOTATED PDB FILE!
        
        string_compat_cdrvic = []
        for q in cdr_vicinity:
            new_tup = (q[0],(str(q[1][0]),q[1][1]))
            string_compat_cdrvic.append(new_tup)

        gen_imgt_cdrs = ['27','28','29','30','31','32','33','34','35','36','37','38',
                         '56','57','58','59','60','61','62','63','64','65',
                         '105','106','107','108','109','110','111','112','113','114','115','116','117'] 


        with open(input_file,"r") as f:
            with open(input_file[:-4] + "_Annotated.pdb","w") as fo:
                for line in f.readlines():
                    if line[0:4] == "ATOM":
                        line = line.strip()
                        residue_id = str(line[23:27].replace(" ",""))
                        res,ins = stripInsertionCode(residue_id)
                        res_tuple = (line[21],(res,ins))
                        
                        if res in gen_imgt_cdrs:
                            fo.write("HETATM"+line[6:55])
                        else:
                            fo.write(line[0:55])


                        if res_tuple in string_compat_cdrvic:

                            if line[17:20] == "LYS" or line[17:20] == "ARG":
                                fo.write(" 1.00 ")
                            elif line[17:20] == "HIS":
                                fo.write(" 0.10 ")
                            elif line[17:20] == "GLU" or line[17:20] == "ASP":
                                fo.write("-1.00 ")
                            else:
                                fo.write("    0 ")
                            fo.write("{0:.3f}".format(hydrophobicity_annotation[line[17:20]])+"\n")
                        else:

                            fo.write("    0     0\n")        

        return stats

##############################################################################
#CODE STARTS HERE#

if __name__ == '__main__':

        #HYDROPHOBICITY SCORING SYSTEMS
        #0: Gravy - Kyte & Doolttle (1982)
        #1: Wimley & White (1996)
        #2: Hessa et al. (2005)
        #3: Eisenberg & McLachlan (1986)
        #4: Black & Mould (1990)

        #USAGE:
        #python Annotator.py [hydrophobic scoring system] [input_file] [ab chains] [region_definition] [output_file]
        #python Annotator.py 0 examples/cristian.pdb IG results.txt

        pH = 7.4    #physiological, here predetermined but could be made an argument
        annotation = sys.argv[1]    #scoring system
        input_file = sys.argv[2]    #input file name
        chains = sys.argv[3]        #Usually IG
        region_definition = sys.argv[4]  #Chothia or IMGT

        CreateAnnotation(int(annotation),float(pH),input_file,chains,region_definition) #Start Main Code
