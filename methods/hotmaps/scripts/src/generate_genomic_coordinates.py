'''
Created on Feb 2, 2017

@author: fran
'''

import MySQLdb
import argparse
import sys
import os



def arguments():
    info = 'Generate file with the genomic coordinates of each structure in the input file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-s', '--structures',
                        type=str, required=True,
                        help='File with the structure IDS. One line per entry.')
    parser.add_argument('--host',
                        type=str, required=False,default="localhost",
                        help='MySQL host')
    parser.add_argument('--db',
                        type=str, default='mupit_modbase',
                        help='MySQL MuPIT database name (default: mupit_modbase)')
    parser.add_argument('--mysql-user',
                        type=str, required=False,default="fran",
                        help='MySQL user name')
    parser.add_argument('--mysql-passwd',
                        type=str, required=False,default="platano",
                        help='MySQL password')
    parser.add_argument('-o', '--output-file',
                        type=str, required=True,default="melon",
                        help='File to output the result')
    args = parser.parse_args()
    return vars(args)

def calculate_segments(l):
    '''
    Given a list of coordinates return the segments of the nucleotides 
    
    '''
    l_output = []
    #list_coordinates = list(set(l))
    list_coordinates = list(l)
    if len(list_coordinates) % 3 != 0:
        print "ERROR"
        
        return []
        
    list_coordinates.sort(key=long)
    
    previous = list_coordinates[0]
    first = previous
    i = 1
    while(i < len(list_coordinates)):
        
        if(not(previous+1==list_coordinates[i])):
            
            l_output.append((first,previous))
            first = list_coordinates[i]
        
        previous = list_coordinates[i]
        i = i+1
    l_output.append((first,previous))
    return l_output
            
def generate_file_genomic_information(list_pdb_ids,opts,dict_d):
    
    """
    Function to parse command line arguments
    from the user

    Parameters
    
    list_pdbs_ids: List of the pdb ids to generate the genomic coordinates
    file_output: Output file to write the segments 
    dict_d: dictionary with non-redudant chains per structure. 
    
    Returns
    -------
    Format output. Generate only one entry per unique chain. 
    
    chr1    Start    End    Strand    PDB_ID    CHAIN_ID    OTHER

    """
    
    
    # make mysql connection
    db = MySQLdb.connect(host=os.environ['MYSQL_HOST'],
                         port=int(os.environ['MYSQL_PORT']),
                         user=os.environ['MYSQL_USER'],
                         passwd=os.environ['MYSQL_PASSWD'],
                         db=os.environ['MYSQL_DB'])
    cursor = db.cursor()
    # query for getting all the genomic information of the structure
    out_f = opts["output_file"]
    f = open(out_f,'a')
    for id_struct in list_pdb_ids:
        dict_chains ={}
        dict_res = {}
        myquery = (
            "select * from Genome2PDB where PDBID LIKE '{pdb_id}_%';"
       
        ).format(pdb_id=id_struct)
        cursor.execute(myquery)
        chr_orig = ""
        # iterate through all mappings
        for result in cursor.fetchall():
            (chr_id,complete_id,res,resInt,c1,c2,c3) = result
            strand = "+"
            if(chr_orig==""):
                chr_orig = str(chr_id)
            chain = complete_id[len(complete_id)-1]
            if(not(chain in dict_d[id_struct]) or chr_orig!=chr_id): # IF the chain is redudant with other chains in the structure not include it
                
                
                continue
            chr= chr_id.replace('chr','')
            if(c3<c1):
                strand = "-"
               
            if(chain in dict_chains): # The chain had been already included append a new residue to the chain
                if(res in dict_res[chain]): # There are two positions for this residue!, DO NOT INCLUDE IT THEN!
                    continue
                dict_res[chain].append(res)
                dict_chains[chain].append((res,chr,c1,c2,c3,strand))
            else:
                dict_res[chain]= [res]
                dict_chains[chain] = [(res,chr,c1,c2,c3,strand)]
        dict_present = {}    
        # Now generate the continuous segments
        list_coordinates = []
        for key in dict_chains.keys():
            # Sort the tuples
            dict_present[id_struct+"_"+key] = []
            strand = dict_chains[key][0][5]
            chr = dict_chains[key][0][1]   
            list_coordinates = []
            dict_res = {}
            for (res,chr,c1,c2,c3,strand) in dict_chains[key][0:len(dict_chains[key])]:
            
                list_coordinates.append(c1)
                list_coordinates.append(c2)
                list_coordinates.append(c3)
                
            list_segments = calculate_segments(list_coordinates)
            if(len(list_segments)) == 0:
                print key,id_struct
                sys.exit(1)     
                
            for (beg,end) in list_segments:     
                dict_present[id_struct+"_"+key].append((chr,beg,end,strand))
            # Print the values
      
        for key in dict_present.keys():
            
            protein=key.split("_")
            if len(protein) == 2:
                struct = protein[0]
                chain = protein[1]
            else:
                struct = "_".join(protein[0:len(protein)-1])
                chain = protein[len(protein)-1]
            sep = "\t"
            for (chr,first,last,strand) in dict_present[key]: 
                
                f.write( chr+sep+str(first)+sep+str(last)+sep+strand+sep+struct+sep+chain+"\n")
    f.close()
    cursor.close()
    db.close()
               
def generate_unique_chains(list_pdbs,opts):
   
    # make mysql connection
    db = MySQLdb.connect(host=os.environ['MYSQL_HOST'],
                         port=int(os.environ['MYSQL_PORT']),
                         user=os.environ['MYSQL_USER'],
                         passwd=os.environ['MYSQL_PASSWD'],
                         db=os.environ['MYSQL_DB'])
    cursor = db.cursor()
    # query for getting all the genomic information of the structure
    
    dict_d = {}
    dict_descript = {}
    for id_struct in list_pdbs:
        dict_chains ={}
        myquery = (
            "select pdbId,hugo,pdbTitle from PDB_Info where pdbId LIKE '{pdb_id}_%';"
       
        ).format(pdb_id=id_struct)
        cursor.execute(myquery)
        dict_d[id_struct] = []
        dict_descript[id_struct] = []
        # iterate through all mappings
        for result in cursor.fetchall():
            (pdb_id,hugo,descript) = result
            pdb = pdb_id[0:4]
            chain = pdb_id[len(pdb_id)-1]
            if(hugo in dict_descript[id_struct] ):
                continue
            else:
                dict_d[id_struct].append(chain)
                dict_descript[id_struct].append(hugo)

    cursor.close()
    db.close()

    return dict_d 
   
def check_proteins_present_outputfile(list_pdbs,fileoutput):
    # Check
    set_pdbs = set()
    if(not(os.path.exists(fileoutput))):
        return list_pdbs
    f = open(fileoutput)
    for line in f:
        line = line.rstrip()
        data = line.split("\t")
        pdb = data[4]
        set_pdbs.add(pdb)
    
    f.close()
    list_new_pdbs = []
    for pdb in list_pdbs:
        if pdb in set_pdbs:
            continue
        list_new_pdbs.append(pdb)
    return list_new_pdbs
def main(opts):
    
    file_data = opts["structures"]
    f = open(file_data)
    list_pdbs = []
    for line in f: 
        line = line.rstrip()
        list_pdbs.append(line)
          
    f.close()
    list_pdbs_n = check_proteins_present_outputfile(list_pdbs,opts["output_file"])
    
    dict_d = generate_unique_chains(list_pdbs_n,opts)
    
    generate_file_genomic_information(list_pdbs_n,opts,dict_d)     
 



if __name__ == '__main__':
    opts = arguments()
    
    main(opts)