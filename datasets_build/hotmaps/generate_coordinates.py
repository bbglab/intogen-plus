import collections
import functools
import multiprocessing
import os
import sqlite3
import sys
import tqdm


def calculate_segments(l):
    """Given a list of coordinates return the segments of the nucleotides"""
    l_output = []
    # list_coordinates = list(set(l))
    list_coordinates = list(l)
    if len(list_coordinates) % 3 != 0:
        print("ERROR")
        return []

    list_coordinates.sort()

    previous = list_coordinates[0]
    first = previous
    i = 1
    while i < len(list_coordinates):

        if not (previous + 1 == list_coordinates[i]):
            l_output.append((first, previous))
            first = list_coordinates[i]

        previous = list_coordinates[i]
        i = i + 1
    l_output.append((first, previous))
    return l_output


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable
    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk


def generate_unique_chains(pdbs, db_file):
    # make sqlite connection
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    result = []

    for id_struct in pdbs:

        myquery = (
            "select * from Genome2PDB where PDBID LIKE '{pdb_id}_%';"

        ).format(pdb_id=id_struct)
        cursor.execute(myquery)
        mappings = cursor.fetchall()

        if len(mappings) == 0:
            # do not go further
            continue

        # query for getting all the genomic information of the structure
        chains = set()
        seen_hugo = set()
        myquery = (
            "select pdbId,hugo,pdbTitle from PDB_Info where pdbId LIKE '{pdb_id}_%';"

        ).format(pdb_id=id_struct)
        cursor.execute(myquery)
        # iterate through all mappings
        for pdb_id, hugo, descript in cursor.fetchall():
            chain = pdb_id[-1]
            if hugo in seen_hugo:
                continue
            else:
                chains.add(chain)
                seen_hugo.add(hugo)

        if len(chains) == 0:
            # do not go further
            continue

        dict_chains = collections.defaultdict(list)
        seen_res = collections.defaultdict(set)

        # iterate through all mappings
        chr_orig = None
        for chr_id, complete_id, res, resInt, c1, c2, c3 in mappings:
            if chr_orig is None:
                chr_orig = str(chr_id)
            chain = complete_id[- 1]
            if chain not in chains or chr_orig != chr_id:  # IF the chain is redudant with other chains in the structure not include it
                continue
            chr = chr_id.replace('chr', '')
            strand = '-' if c3 < c1 else '+'

            if res in seen_res[chain]:  # There are two positions for this residue!, DO NOT INCLUDE IT THEN!
                continue
            seen_res[chain].add(res)
            dict_chains[chain].append((res, chr, c1, c2, c3, strand))

        dict_present = collections.defaultdict(list)
        # Now generate the continuous segments
        for chain in dict_chains.keys():
            # Sort the tuples
            list_coordinates = []

            for (res, chr, c1, c2, c3, strand) in dict_chains[chain]:
                list_coordinates.append(c1)
                list_coordinates.append(c2)
                list_coordinates.append(c3)

            list_segments = calculate_segments(list_coordinates)
            if (len(list_segments)) == 0:
                print(chain, id_struct)
                sys.exit(1)

            for beg, end in list_segments:
                result.append((chr, str(beg), str(end), strand, id_struct, chain))

    conn.close()
    return result


def main(db_file, structures_file, output_file, cpus):
    pdbs = []
    with open(structures_file) as f:
        for line in f:
            line = line.rstrip()
            pdbs.append(line)

    chainer = functools.partial(generate_unique_chains, db_file=db_file)

    with open(output_file, 'w') as f, multiprocessing.Pool(processes=cpus) as pool:
        for r in tqdm.tqdm(pool.imap_unordered(chainer, chunkizator(pdbs, size=100)), total=len(pdbs)//100+1):
            for val in r:
                f.write("\t".join(val) + "\n")


if __name__ == '__main__':
    args = sys.argv[1:]
    db_file = args[0]
    structures_file = args[1]
    output_file = args[2]
    if len(args) > 3:
        cpus = int(args[3])
    else:
        cpus = 1

    main(db_file, structures_file, output_file, cpus)
