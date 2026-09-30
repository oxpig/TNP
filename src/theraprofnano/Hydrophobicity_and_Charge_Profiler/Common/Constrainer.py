import urllib.request, urllib.error, urllib.parse
import os,sys
import random
from optparse import OptionParser
import tempfile
import anarci

#Used to determine if a residue is an insertion
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Check if file exists and where it is. Return a full path to it
def check_file(filename):
	#See if it is in the current directory
	cwd = os.getcwd()
	try:
		with open(cwd+"/"+filename,'r'):
			return cwd+"/"+filename
	except IOError:
		#Not in the current directory, check if it is an absolute path then
		try:
			with open(filename,'r'):
				return filename
		except IOError:
			print("File ",filename," does not exist.")
			cleanup()

#Return the sequence of the chain together with the map
def get_sequence(filename,t_chain):
	
	f = open(filename,'r')
	prev_res = ""
	seq = ""
	res_map = dict()
	resnames = []
	for line in f.readlines():
		if line[0:4]=="ATOM":
			chain = line[21]
			if chain != t_chain:
				continue
			line = line.strip()
			resname=line[23:29].strip()
			if resname == prev_res:
				continue
			else:
				prev_res = resname
				
			AA_3 = str(line[17:20])
			AA_1 = AAtable[AA_3]
			resnames.append(resname)
			res_map[resname] = AA_1 
			seq+=AA_1
	f.close()
	return res_map,seq,resnames
	
#Change the b-factor on the line
def color_line(line,bf):
	while(len(bf))<6:
		bf=" "+bf
	line = line[0:60]+bf+line[66:len(line)]
	return line

def rename_line(chid,line):
	chain = chid[0]
	line = line[0:21]+chain+line[22:len(line)]
	
	chid = chid[1:len(chid)]
	try:
		chid = str(int(chid))
		while len(chid)<4:
			chid=' '+chid
		line = line[0:22]+chid+line[26:len(line)]
		
	except ValueError:
		#There is an insertion code
		ic = chid[len(chid)-1]
		chid = chid[0:len(chid)-1]
		
		while len(chid)<4:
			chid=' '+chid
		line = line[0:22]+chid+ic+line[27:len(line)]
	
	return line
	

#only save the residues of the file that are in thelist
def save_colored(filename,chains,mapping,folder):
	"Transform mapping into a real thing..."
	f = open(filename,'r')
	f_w = open(folder+"/red_blue.pdb",'w')
	prev_res = ""
	seq = ""
	res_map = dict()
	resnames = []
	for line in f.readlines():
		if line[0:4]=="ATOM":
			chain = line[21]
			
			if chain not in chains:
				continue
			line = line.strip()
			resname=line[23:29].strip()
			try:
				line = rename_line(id_map[resname+' '+chain],line)
			except KeyError:#This means that the structure is much bigger.
				pass
			resname = resname+" "+chain
			if resname in mapping:
				bf = "100.00"
			else:
				bf = "0.00"
			#line = color_line(line,bf)			
			f_w.write(line+"\n")
				
	f_w.close()
	f.close()	
	return res_map,seq,resnames		
			
#Create a folder with a specified name
def create_folder(fname):
	cwd = os.getcwd()
	foldname = cwd+"/"+fname
	if os.path.exists(foldname): 
		return False			
	else:
		os.makedirs(foldname)
		return foldname

#Create a temporary folder in the user directory - return it so that we can do whatever with it
def create_temp_folder():
	dirpath = tempfile.mkdtemp()
	return dirpath
	
#remove the insertion codes from Chothia codes
def stripChothia(res):
	icode = "ZXCVBNMLKJHGFDSAWQERTYUIOP"
	if res[len(res)-1] in icode:
		res = res[0:(len(res)-1)]
	return res

# Changed to use anarci as default
def chothia_number(tfold,seq):
	Numbering,ChainType  = anarci.number(seq,scheme="chothia")
	print(Numbering, ChainType)
	mapping = []
	
	for entry in Numbering:
		insert = ''
		if entry[0][1]!=' ':
			insert = entry[0][1]
		key = ChainType+str(entry[0][0])+insert
		if entry[1] == '-':
			continue
		mapping.append([key,entry[1]])
	
	return mapping


#Cleanup the temporary files and quit.
def cleanup():
	os.system("rm -rf "+tempfold)
	quit()
#Check if a residue is a CDR residue according to the given definition (deff)
def is_CDR(res,deff):
	
	for CDR in definitions[deff]:
		if res in definitions[deff][CDR]:
			return CDR

	return False

def write_data(filename,data):
	
	f_out = open(filename,'w')
	f_out.write(data)
	f_out.close()	

##MAIN##
#AA mapping
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

#CDR Definitions # NOTE GG removed L chain definitions
definitions = {		
		"kabat" : {"H1" : ["H31", "H32", "H33", "H34", "H35"],
				"H2" : ["H50", "H51", "H52", "H53", "H54", "H55", "H56", "H57", "H58", "H59", "H60", "H61", "H62", "H63", "H64", "H65"],
				"H3" : ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]}
		,
		"chothia" : {"H1" : ["H26", "H27", "H28", "H29", "H30", "H31", "H32"], 
					"H2" : ["H52", "H53", "H54", "H55", "H56"],
					"H3" : ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]} 
		,
		"contact" : {"H1" : ["H30", "H31", "H32", "H33", "H34", "H35"], 
					"H2" : ["H47", "H48", "H49", "H50", "H51", "H52", "H53", "H54", "H55", "H56", "H57", "H58"],
					"H3" : ["H93", "H94", "H95", "H96", "H97", "H98", "H99", "H100", "H101"] }
		}

#All the temporary files are stored here
tempfold = create_temp_folder()
cwd = os.getcwd()

##Read In options##
#Parse the options
usage = "USAGE: python Framer.py --f FILE --c CHAINS"
parser = OptionParser(usage=usage)

#Single file mode
parser.add_option("--f",help="Nanobody file location", dest="file")
parser.add_option("--c",help="Nanobody chain to constrain, e.g. --c ABCD", dest="chains")
parser.add_option("--o",help="Output directory - where to write the results", dest="out_dir")
parser.add_option("--d",help="Definition to use - kabat, chothia or contact", dest="deff")

(options, args) = parser.parse_args()

if (options.file and options.chains and options.out_dir and (options.deff=="kabat" or options.deff=="chothia" or options.deff=="contact")):
	
	filename = check_file(options.file)
	#Attempt to create the output directory
	
	chain_maps = dict()
	output_cdr = ""
	cdr_detailed = "Original ID	Original Chain	AA	Chothia ID	CDR(FR=frame,or CDR id)\n"
	output_frame = ""
	cdr_map = []
	id_map = dict()
	for chain in options.chains:
		chain_maps[chain],seq,resnames = get_sequence(filename,chain)
		#Chothia number the sequence
		ch_map = chothia_number(tempfold,seq)
		zipped = list(zip(resnames,ch_map))
		for elem in zipped:
			
			cdr = is_CDR(stripChothia(elem[1][0]),options.deff)
			id_map[elem[0]+" "+chain] = elem[1][0]
			if cdr!=False:
				output_cdr+= elem[0]+" "+chain+"\n"
				cdr_map.append(elem[0]+" "+chain)
			else:
				cdr = "FR"
				output_frame+= elem[0]+" "+chain+"\n"
			cdr_detailed+=	elem[0]+"\t"+chain+"\t"+elem[1][1]+"\t"+elem[1][0]+"\t"+str(cdr)+"\n"	
			#Sanity check
			#print elem, chain_maps[chain][elem[0]]
			if chain_maps[chain][elem[0]] != elem[1][1]:
				
				print("Something wrong in the mapping:")
				cleanup()				

	
	#Save the pdb file with colored CDRs
	write_data(options.out_dir+"/full_info.txt",cdr_detailed)
	write_data(options.out_dir+"/paratope.txt",output_cdr)
	write_data(options.out_dir+"/framework.txt",output_frame)
	save_colored(filename,options.chains,cdr_map,options.out_dir)
	#print "The CDRs were determined using the ",options.deff," definition.\n The results can be found in ",options.out_dir
			
else:
	print("There's something missing from your input!")

cleanup()



