### GG modified from https://github.com/oxpig/SAbDab/blob/master/lib/python/ABDB/AB_Utils/sequence_liabilities.py, added polyreactivity motifs

import os
import re
from ABDB.AB_Utils.region_definitions import Accept
from anarci import anarci, anarci_output

# annotate_sequences accounts for if only VH or VL
def annotate_sequences(heavy_sequence, light_sequence, scheme="imgt", restrict_species=True):
    """
    Annotate the sequence. Use the IMGT scheme to do the alignments. Report to the log

    Parameters
    ----------
    heavy_sequence : str
            Antibody VH chain
    light_sequence : str
            Antibody VL chain
    scheme : str, default scheme='imgt'
            ANARCI numbering scheme to annotate sequence
    restrict_species : bool, default restrict_species=True
            If True, this will restric germline species for numbering to default ANARCI values (human and mouse)
            Otherwise, this will override default species
        
    Return
    ------
        lists : ANARCI annotated heavy and light chains

    """
    # Use the IMGT numbering scheme for the modelling
    sequences = [ ("heavy_sequence", heavy_sequence), ("light_sequence", light_sequence)]
    if restrict_species:
        numbered, alignment_details, hit_tables = anarci( sequences, scheme=scheme, allowed_species=None)
    else:
        numbered, alignment_details, hit_tables = anarci( sequences, scheme=scheme)

    heavy_numbered, light_numbered = numbered
    heavy_alignment_details, light_alignment_details = alignment_details
    heavy_hit_tables, light_hit_tables = hit_tables
    hchtype, lchtype = "", ""
        
    if (heavy_numbered is None or len(heavy_numbered) > 1 or heavy_alignment_details[0]["chain_type"] != "H"):
        ##self.message( "Heavy chain sequence was not identified as having a single antibody heavy domain." )
        if heavy_alignment_details:
            hchtype = heavy_alignment_details[0]["chain_type"]
    else:
        start = heavy_alignment_details[0]["query_start"]
        end   = heavy_alignment_details[0]["query_end"]
        hchtype = heavy_alignment_details[0]["chain_type"]

    if light_numbered is None or len(light_numbered) > 1 or light_alignment_details[0]["chain_type"] not in "KL":
        if light_alignment_details:
            lchtype = light_alignment_details[0]["chain_type"]
    else:
        start = light_alignment_details[0]["query_start"]
        end   = light_alignment_details[0]["query_end"]
        #self.message( "-"*start + "*"*(end-start) + "-"*(len(light_sequence)-end) )
        #self.message( "Trimming light sequence to the variable domain (*) only." )
        lchtype = light_alignment_details[0]["chain_type"]

    # Check to make sure that there's at least one numbered sequence. JL 28-07-2015
    if heavy_numbered and light_numbered:
        return heavy_numbered[0], light_numbered[0], hchtype, lchtype
    elif heavy_numbered and not light_numbered:
        #self.message("Single domain antibody identified; heavy chain only.")
        return heavy_numbered[0], None, hchtype, lchtype
    elif not heavy_numbered and light_numbered:
        #self.message("Single domain antibody identified; light chain only.")
        return None, light_numbered[0], hchtype, lchtype
    else:
        #self.message( "No antibody chains have been identified in the sequence. Modelling halted." )
        return None, None, hchtype, lchtype


# Macros for the verniers and the nterminus
positions_macros={"verniers":{"H":[ (2," "), (28," "), (29," "), (54," "), (55," "), (78," "), (88," "), (105," "), (106," "), (118," ") ],
                                "L":[(4," "), (27," "),(28," "),(29," "),(30," "),(31," "),(32," "),(33," "),(34," "),(35," "),(36," "),(41," "),(42," "),(52," "),(53," "),(55," "),(84," "),(94," "),(118," ")]},
                    "nterminus":{"H":[(1," ")],"L":[(1," ")]}}


# Get a way of telling if the identified position is in a region of interest
def get_acceptor(regions, chain, scheme="imgt"):
    """
    choth_hverniers = [ (2," "), (27," "), (28," "), (49," "), (50," "), (69," "), (79," "), (93," "), (94," "), (103," ") ]
    imgt_hverniers = [ (2," "), (28," "), (29," "), (54," "), (55," "), (78," "), (88," "), (105," "), (106," "), (118," ") ]
    choth_lverniers = [ (4," "), (27," "), (28," "), (29," "), (30," "), (35," "), (36," "), (46," "), (47," "), (49," "), (68," "), (78," "), (98," ")]
    imgt_lverniers =[(4," "), (27," "),(28," "),(29," "),(30," "),(31," "),(32," "),(33," "),(34," "),(35," "),(36," "),(41," "),(42," "),(52," "),(53," "),(55," "),(84," "),(94," "),(118," ")]
    """
    accept = Accept( numbering_scheme=scheme, definition="north" ) # We use the north still - used to use kabat.
    for region in regions:
        # Add positions to the accept class using the custom macros above or those built in to it.
        if region in positions_macros:
            accept.add_positions(positions_macros[region][chain], chain)
        else:
            accept.add_regions( [region] )
    return accept


# NOTE edit here to deal with if no light_sequence input
def get_liabilities(heavy_sequence, light_sequence, scheme="imgt", save=False, outfile="sequence_liabilities.csv", restrict_species=True):
    """
    Returns list of sequence liabilities provided input heavy and light sequences
    
    Parameters
    ----------
    heavy_sequence : str
            Antibody VH chain
    light_sequence : str
            Antibody VL chain
    scheme : str, default scheme='imgt'
            ANARCI numbering scheme to annotate sequence
    restrict_species : bool, default restrict_species=True
            If True, this will restric germline species for numbering to default ANARCI values (human and mouse)
            Otherwise, this will override default species
        
    Return
    ------
        list of tuples 
            found_liabilities in input VH/VL sequences

        int
            number of defined liabilities in CSV file 
    """

    # Generate CSV file if not present in same dir
    liability_path = os.path.join( os.path.split(__file__)[0], "liabilities.csv" )
    if not os.path.isfile(liability_path):
        make_liabilities_csv()
        liability_path = os.path.join( os.path.split(__file__)[0], "liabilities.csv" )
    
    liabilities = []
    with open( liability_path ) as lfile:
        for line in lfile:
            if stripped := line.strip():
                if stripped[0] != "#":
                    name, region, regex = stripped.split(",")
                    liabilities.append( [name, region.split(";"), re.compile( regex ) ] )

    # annotate_sequences deals with if nanobody
    heavy_numbered, light_numbered, hchtype, lchtype = annotate_sequences( heavy_sequence, light_sequence, scheme=scheme, restrict_species=restrict_species)
                
    if heavy_numbered:
        heavy_numbered = ([ (n,a) for n,a in heavy_numbered[0] if a != "-" ], heavy_numbered[1], heavy_numbered[2])
    if light_numbered:
        light_numbered = ([ (n,a) for n,a in light_numbered[0] if a != "-" ], light_numbered[1], light_numbered[2])

    # NOTE changed here:
    if heavy_numbered and light_numbered:
        target_numbering = {"H": heavy_numbered[0], "L": light_numbered[0]}
        light_numbering=target_numbering["L"]
        heavy_numbering=target_numbering["H"]
    elif heavy_numbered and not light_numbered:
        target_numbering = {"H": heavy_numbered[0], "L": light_numbered}
        light_numbering=target_numbering["L"]
        heavy_numbering=target_numbering["H"]
    elif light_numbered and not heavy_numbered:
        target_numbering = {"H": heavy_numbered, "L": light_numbered[0]}
        light_numbering=target_numbering["L"]
        heavy_numbering=target_numbering["H"]


    # Other variables
    found_liabilities, termHglut = [], False
    s, i2p = {"H":"","L":""}, {"H":{}, "L":{}}
    numbering = {"H":heavy_numbering,"L":light_numbering}

    # Convert the numbering back to a sequence and map the indices onto the positions for lookup
    for chain in "HL":
        i = 0
        if numbering[chain] != None: # NOTE here account for if nanobody
            for p , aa in numbering[chain]:
                if aa != "-":
                    s[chain]+= aa 
                    i2p[chain][i] = p            
                    i+=1

    # Run through each of the liabilities find if they are present
    for li, (liability, regions, motif) in enumerate(liabilities):
        # Do for both chains
        for chain in "HL": 
            if numbering[chain] != None: # NOTE here account for if nanobody
                acceptor_obj = get_acceptor( regions, chain )
                # Find hits using the regex
                for hit in motif.finditer(s[chain]):
                    # Now find whether we care about it using the rules of the scoring regions.   
                    start, end = hit.start(), hit.end()-1         
                    pos_start, pos_end = i2p[chain][start], i2p[chain][end]
                    if liability == "Unpaired Cys (C)" and pos_start in [(23, " "),(104, " ")]: # ignore the known cysteines.
                        continue
                    if acceptor_obj.accept( pos_start, chain ):
                        if liability == "N-terminal glutamate (VH and VL) (E)": # requires both heavy and light chain
                            if chain=="L" and termHglut:
                                found_liabilities.extend([[liability, start, end, list(pos_start), list(pos_end), chain, li], termHglut])

                            else:
                                # fix 04-03-2016; termHglut has to be in the same syntax as everything else in found_liabilities
                                termHglut=[liability, start, end , list(pos_start), list(pos_end), chain, li ]
                        else:
                            found_liabilities.append( [liability, start, end, list(pos_start), list(pos_end) , chain, li ] )

    if save:
        i2pa = {}
        for chain in "HL":
            if numbering[chain] != None: # NOTE here account for if nanobody
                i2pa[chain] = {}
                i = 0
                for p , aa in numbering[chain]:
                    if aa != "-":
                        i2pa[chain][i] = (p,aa)
                        i+=1

        # Write liability file
        with open( os.path.join( outfile ) ,'w' ) as lfile:
            print(",".join(["chain", "position", "aa", "liability"]), file=lfile)
            for liability, start, end, _, _, chain, li in found_liabilities:
                for i in range( start , end+1 ):
                    print(",".join([chain, "%d%s"%(i2pa[chain][i][0]), i2pa[chain][i][1], liability]), file=lfile)


    return found_liabilities, len(liabilities)


# create csv of liabilities and their frequency
def make_liabilities_csv():
    """Function to regenerate liabilities.csv file if unavailable in the same dir"""

    content = ("# Define liabilities, name, region(s), regex motif\n"
            "Unpaired Cys (C),fv,C\n"
            "N-linked glycosylation (NXS/T X not P),fv,N[^P][ST]\n"
            "Met oxidation (M),cdrs;verniers,M\n"
            "Trp oxidation (W),cdrs;verniers,W\n"
            "Asn deamidation (NG NS NT),cdrs;verniers,N[GST]\n"
            "Asp isomerisation (DG DS DT DD DH),cdrs;verniers,D[GSTDH]\n"
            "Lysine Glycation (KE KD EK ED),cdrs;verniers,KE|KD|EK|ED\n"
            "N-terminal glutamate (VH and VL) (E),nterminus,E\n"
            "Integrin binding (RGD RYD LDV),fv,RGD|RYD|LDV\n"
            "CD11c/CD18 binding (GPR),fv,GPR\n"
            "Fragmentation (DP),cdrs;verniers,DP\n"
            "Polyreactivity (RR VG VV VVV WW WWW WXW),fv,RR|VG|VV|VVV|WW|WWW|WXW\n")

    liability_path = os.path.join( os.path.split(__file__)[0], "liabilities.csv" )
    with open(liability_path,'w') as f:
        f.write(content)
    f.close()
