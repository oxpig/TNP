import math
from scipy.spatial.distance import cdist

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

class Atom:
	def __init__(self,typ,coords,line):
		self.typ = typ
		self.coords = coords
		self.pdb_line = line
	"""
	To avoid "PDBUtils.Atom instance" message...	
	def __repr__(self):
		return "%s" % (self.typ)
	"""
		
#Class to hold a single chain from the PDB
class Residue:
	def __init__(self,atoms,res_type):
		self.atoms = atoms
		self.res_type = res_type
	
	def distance(self,res1,res2):
		mindist = 99999999
		for a1 in res1.atoms:
			for a2 in res2.atoms:
				dist = self.euclid(a1.coords[0]-a2.coords[0],a1.coords[1]-a2.coords[1],a1.coords[2]-a2.coords[2])
				if dist<mindist:
					mindist = dist
					
		return mindist
	
	def sb_distance(self,res1,res2):
		mindist = 99999999
		donors = []
		acceptors = []
		if res1.res_type == "LYS":
			for atom in res1.atoms:
				if atom.typ == "NZ":			
					donors.append(atom)
			if res2.res_type == "ASP":
				for atom in res2.atoms:
					if atom.typ == "OD1" or atom.typ == "OD2":
						acceptors.append(atom)
			elif res2.res_type == "GLU":
				for atom in res2.atoms:
					if atom.typ == "OE1" or atom.typ == "OE2":
						acceptors.append(atom)
		elif res1.res_type == "ARG":
			for atom in res1.atoms:
				if atom.typ == "NH1" or atom.typ == "NH2":
					donors.append(atom)
			if res2.res_type == "ASP":
				for atom in res2.atoms:
					if atom.typ == "OD1" or atom.typ == "OD2":
						acceptors.append(atom)
			elif res2.res_type == "GLU":
				for atom in res2.atoms:
					if atom.typ == "OE1" or atom.typ == "OE2":
						acceptors.append(atom)					
		"""
		print "*******************************************"
		
		for a in donors:
			print a.typ
			print a.coords
			print a.pdb_line
			
		for b in acceptors:
			print b.typ
			print b.coords
			print b.pdb_line
			
		print "*******************************************"		
		"""
		interacters = False
		
		for a1 in donors:
			for a2 in acceptors:
				dist = self.euclid((a1.coords[0]-a2.coords[0]),(a1.coords[1]-a2.coords[1]),(a1.coords[2]-a2.coords[2]))
				if dist<mindist:
					mindist = dist
					interacters = (a1,a2)
		
		#if interacters != False:
		#	print "Minimum Distance between " + str(interacters[0].typ) + " on " + str(res1.res_type) + " coord: " + str(interacters[0].coords) + " and " + str(interacters[1].typ) + " on " + str(res2.res_type) + " coords: " + str(interacters[1].coords) + " is: " + str(mindist)		
		
		return mindist				

	def euclid(self,x,y,z):
		result = math.sqrt(x*x+y*y+z*z)
		
		return result

class PDBchain:
	
	#Get the contact mappings between the two structures
	def contact_map(self,str1,str2):
		mapping = dict()
		types = dict()
		for chain_1 in str1:
			for chain_2 in str2:
				for res_id_1 in str1[chain_1].residues:
					r_1 = str1[chain_1].residues[res_id_1]
					index_id_1 = chain_1+str(res_id_1[0])+res_id_1[1]
					if index_id_1 not in mapping:
						mapping[index_id_1] = dict()
					types[index_id_1] = r_1.res_type
					for res_id_2 in str2[chain_2].residues:
						index_id_2 = chain_2+str(res_id_2[0])+res_id_2[1]
						r_2 = str2[chain_2].residues[res_id_2]
						types[index_id_2] = r_2.res_type
						d = r_1.distance(r_1,r_2)
						mapping[index_id_1][index_id_2] = d
	
		return mapping,types
	#Change the b-factor on the pdb line
	def color_line(self,line,bf):
		while(len(bf))<6:
			bf=" "+bf
		line = line[0:60]+bf+line[66:len(line)]
			
		return line
	
	def __init__(self,path_to_pdb,which_chain):
		
		try:
			residues = dict()
			atoms = []
			curr_id = ""
			res_type = ""
			for line in open(path_to_pdb):
				
				if line[0:4]!='ATOM':
					continue
				line = line.strip()
				#print "Chain:",line[21]
				chain = line[21]
				if chain!=which_chain:
					continue
				#print line
				#print "Sid:",line[22:28].replace(" ","")
				sid = line[22:28].replace(" ","")
				
				if  is_number(sid):
					sid = (int(sid),'')
				else:
					sid = (int(sid[0:len(sid)-1]),sid[len(sid)-1])
				if len(curr_id) ==0:
					curr_id = sid
				#print "atom: ",line[12:17].replace(" ","")
				atom_type = 	line[12:17].replace(" ","")			
				
				#print "x:",float(line[30:38].replace(" ",""))
				x = float(line[30:38].replace(" ",""))
				#print "y:",float(line[38:46].replace(" ",""))
				y = float(line[38:46].replace(" ",""))
				#print "z:",float(line[46:54].replace(" ",""))
				z = float(line[46:54].replace(" ",""))
				atom = Atom(atom_type,(x,y,z),line)				
				#Push the previous residue				
				if sid!=curr_id:
					#push the residue
					res = Residue(atoms,res_type)
					residues[curr_id] = res
					curr_id = sid
					atoms = [atom]
				else:
					atoms.append(atom)
				#print "Res: ",line[17:21].replace(" ","")
				res_type = line[17:21].replace(" ","")
				
				
			#Last residue			
			res = Residue(atoms,res_type)
			residues[curr_id] = res
			#print "Initialized with ",len(residues)," residues"
			self.residues = residues		
		except IOError:
			print("IOError: ",path_to_pdb,"does not exist")
			quit()
		

