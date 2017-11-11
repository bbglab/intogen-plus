'''
Created on Apr 3, 2017

@author: fran

Module that given the ouput of generating random mutations in a structure retunrs the mapping of the aminoacids and chains affected by those mutations'''
import MySQLdb
import os

import time


def read_generated_mutations(structure_id,file_mutations,cancer_type,number_mutations):
    """Read the generated mutations.

    Parameters
    ----------
    structure_id : string
        pdb of the structure
    cancer_Type: string
        cancer type of the mutations
    file_mutations: strnig
        path to the file containing the generated mutations
    opts: 
        Parameters for the query to the database
    Returns
    -------
    list_mutations:
        list of tuples containing the mutaitons (RES,CHAIN). Duplicates are allowed. 
    """
    # make mysql connection
    retries = 5
    while retries > 0:
        try:
            db = MySQLdb.connect(host=os.environ['MYSQL_HOST'],
                                 port=int(os.environ['MYSQL_PORT']),
                                 user=os.environ['MYSQL_USER'],
                                 passwd=os.environ['MYSQL_PASSWD'],
                                 db=os.environ['MYSQL_DB'])
            break
        except Exception:
            time.sleep(5)
            retries -= 1

    if retries == 0:
        raise RuntimeError("Impossible to connect")

    cursor = db.cursor()
    # query for getting all the genomic information of the structure
    
    list_results = []
    f = open(file_mutations)
    list_output = []
    header = True
    for line in f: 
        if(header):
            header=False
            continue
        line = line.rstrip()
        data = line.split("\t")
        # CHROMOSOME    POSITION    REF    ALT    PROTEIN_ID    CHAIN    STRAND    SIG_REF    SIG_ALT    POSITION_UPSTREAM    POSITION_DOWNSTREAM
        pos = int(data[1])+1
        
        pdb = data[4]
        chain = data[5]
        #cancer = data[11]
        if pdb == structure_id: # and cancer == cancer_type: 
            list_output.append((pdb,chain,pos))
    f.close()
    
    for (pdbid,chain,pos) in list_output:
        myquery = (
            "select pdbPos from Genome2PDB B, Sprot2PDB A where A.pdbId=B.PDBId and B.seqRes=(A.sprotPos-1)  and A.pdbId = '{pdb_id}_{chain}' and (pos1={pos} or pos2={pos} or pos3={pos});" # THINK ABOUT a better query 
	#	"select seqRes,pos1,pos2,pos3,pdbId from Genome2PDB where pdbId LIKE '{pdb_id}_%';" # THINK ABOUT a better query       
        ).format(pdb_id=pdbid,chain=chain,pos=pos)
        cursor.execute(myquery)
        
        # iterate through all mappings
        for result in cursor.fetchall():# Now get the PDB position 
            pos = int(result[0])
            list_results.append((chain,pos))

    cursor.close()
    db.close()

    return list_results

def map_generated_mutations(structure_id,list_output,d_correspondence):
    """Read the generated mutations.

    Parameters
    ----------
    structure_id : string
    list_results: list of dictionaries output of the randomizer_aa
    opts: 
        Parameters for the query to the database
    Returns
    -------
    list_mutations:
        list of tuples containing the mutaitons (RES,CHAIN). Duplicates are allowed. 
    """

    dict_results = {}

    header = True


    for dict in list_output:
        pdbid = dict["PDB_ID"]
        chain = dict["CHAIN"]
        pos = str(dict["POSITION"])
        sim = dict["SIMULATION_ID"]

        try:
            posn = int(d_correspondence[pdbid+"_"+chain][pos])
        except:
            print "Error, position" + pos + " " + pdbid + chain + " Not mapped"
            continue

        if(sim in dict_results):
            dict_results[sim].append((chain,posn))
        else:
            dict_results[sim] = [(chain,posn)]
                
   
    return dict_results


def generate_correspondence(structure_id, db):

    cursor = db.cursor()
    # query for getting all the genomic information of the structure

    d_correspondence = {}
    myquery = (
            #"select pdbPos,pos1,pos2,pos3,A.pdbId from Genome2PDB B, Sprot2PDB A where A.pdbId=B.PDBId and B.seqRes=(A.sprotPos-1)  and A.pdbId LIKE '{pdb_id}_%';" # THINK ABOUT a better query
         	"select seqRes,pos1,pos2,pos3,pdbId from Genome2PDB where pdbId LIKE '{pdb_id}_%';" # THINK ABOUT a better query
        ).format(pdb_id=structure_id)

    cursor.execute(myquery)
    for result in cursor.fetchall():# Now get the PDB position
        try:
            posPDB = int(result[0])
        except:
            posPDB = 1

        pos1 = str(result[1])
        pos2 = str(result[2])
        pos3 = str(result[3])
        pdb_id =result[4]

        if not(pdb_id in d_correspondence):
            d_correspondence[pdb_id] = {}

        d_correspondence[pdb_id][pos1] = posPDB
        d_correspondence[pdb_id][pos2] = posPDB
        d_correspondence[pdb_id][pos3] = posPDB

    cursor.close()

    return d_correspondence

#print read_generated_mutations("1a08","/home/fran/Documents/clusters3D/coordinates/out_test.results","GBM",{"host":"localhost","mysql_user":"fran","mysql_passwd":"platano","db":"mupit_modbase"})    
