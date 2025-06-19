import pickle as pickle
from ABDB.AB_Utils import Accept    #Required to assign CDR or Not in a given definition
from anarci import run_anarci   #Required to number the chains
from optparse import OptionParser   #Required for arguments

##############################################################################
#USE THE FUNCTION accept() WITHIN ANARCI TO ASSIGN CDRs OR FRAMEWORK REGIONS)#

a = Accept()	#Instantiate Accept()

#SPECIFY YOUR DEFAULT NUMBERING SCHEME HERE!

a.numbering_scheme = "imgt"

#############################################################
#Called on each single sequence after ANARCI has numbered it#

def find_all_regions(output_name,numbering,Chain): # Chain should always be H for VHH Nbs

        definitions = ['imgt','north','chothia','kabat']  #list your desired definitions here
	
        regions = {}
        cdr_length_dict = {}

        for definition in definitions:
         
            a.definition = definition
            regions[definition] = {}
            cdr_length_dict[definition] = {}

            if Chain == "H":
                list_of_regions = ["cdrh3","fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3","fwh4"]
                list_of_cdrs = ["cdrh1","cdrh2","cdrh3"]
            else:
                print('ERROR: Chain was not H')


            for r in list_of_regions:
	
                region = ''
		
                a.set_regions([r])		# Feed each region at a time into accept

                for n in range(len(numbering)):
                    if a.accept(numbering[n][0], Chain) == 1:       # If this position exists in the chain in question
                        if numbering[n][1] == '-':		# nothing to add to the string if a gap present at that number!
                            continue
                        else:		
                            region += numbering[n][1]		# use string concatenation to build up the region primary sequence
                            
                    
                regions[definition][r] = region   # We collect the entire primary sequence for each antibody region
                if r in list_of_cdrs:
                    cdr_length_dict[definition][r] = len(regions[definition][r]) 
                    
        # Returns a dictionary of regions as strings.
        return regions,cdr_length_dict

###################################################################
#The ANARCI function#
def anarci_fun(name,sequence, ncores): # NOTE changed allowed_species to any for VHH
    anarci_tuple = (name,sequence)
    output = run_anarci([anarci_tuple], scheme='imgt', assign_germline = True, allowed_species = None, ncpu= ncores) 
    return output

###################################################################
#This is automatically called after __main__#

def main(output_name, sequence, chain, output_dest, ncpu = 1, verbose=True): 
        #aa_dict = {}    #Used to store amino acid by position data
        region_dict = {}	#Used to store strings and lengths of IMGT/North CDRs
        region_dict['imgt'] = {}
        region_dict['north'] = {}

        if chain == "H":
            chain_id = "Heavy"
            list_of_regions = ["cdrh1","cdrh2","cdrh3","fwh1","fwh2","fwh3","fwh4"]
        else:
            print('ERROR: chain type not H')

	#################################
        output = anarci_fun(output_name,sequence, ncpu) #Submit to ANARCI, saving the output!
		    
        assert output[1] != [None]	    #If ANARCI failed, notify user and raise exception
	    
        numbering = output[1][0][0][0]		    #this is the sequence numbering part of the output

	# use accept() to determine CDR or Framework strings
        primary_seqs, length_dict = find_all_regions(output_name,numbering,chain)

        return primary_seqs, length_dict #amino_dict, length_dict

###########################################################
#CODE STARTS HERE#
if __name__ == '__main__':

        parser = OptionParser()		#instantiate the argument parser
	
        #Define arguments here
	
        parser.add_option("-n", "--ncores", dest="ncores", type="int", default=1, help="number of cores")

        parser.add_option("-i", "--input", action="store", type="string", dest="input_sequence")

        parser.add_option("-s", "--species", action="store", type="string", dest="current_species")#, default= "human")

        (options, args) = parser.parse_args()

        print("--------------------------------------")
	
        #if number of cores specified, set it, else use 1
	
        if options.ncores:
            ncores = options.ncores
            print("The number of cores is set to %d " % (ncores))

        if options.input_sequence:
            input_sequence = options.input_sequence
            print("The sequence being analysed is %s " % (input_sequence))

        #Now call the main function to get ANARCI to do some analysis on this sequence!
        main(input_sequence,options.current_species,ncores)

