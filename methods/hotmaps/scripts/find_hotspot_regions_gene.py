import Bio.PDB
import src.pdb_structure as pstruct
import src.utils as utils
import src.graph as graph
import argparse
import csv
import itertools as it

# get logger
import logging
import os
logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = 'Uses BFS to connect hotspot residues into connected regions.'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-m', '--multiple-testing',
                        type=str, required=True,
                        help='File that corrects for multiple testing')
    parser.add_argument('-a', '--annotation-dir',
                        type=str, required=True,
                        help='Annotation directory from CRAVAT')
    parser.add_argument('-p', '--pdb-info',
                        type=str, required=True,
                        help='PDB information file (contains paths to PDBs)')
    parser.add_argument('-r', '--radius',
                        default=10.0,
                        type=float,
                        help='Sphere radius in angstroms for connecting link between two residues (Default: 10.0)')
    parser.add_argument('-q', '--q-value',
                        default=.01,
                        type=float,
                        help='Q-value for FDR (Default: .01)')
    parser.add_argument('-o', '--output',
                        default='output.txt',
                        type=str,
                        help='Output result file for hotspot regions')

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging

    opts = vars(args)
    return opts


def read_delim(path):
    """Read in tab delimited file."""
    data = []
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        data = list(myreader)
    return data


def read_mupit_file(path, signif_res):
    """Reads in the mupit annotation file, but only the lines corresponding
    to significant residues. This reduces memory usage by a lot.

    Parameters
    ----------
    path : str
        path to mupit annotation file
    signif_res :
        significant residues
    """
    data = []
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        anot_header = next(myreader)
        pdb_ix = anot_header.index('pdb_id')
        gene_ix = anot_header.index('HUGO symbol')
        tx_ix = anot_header.index('Reference Transcript')
        res_ix = anot_header.index('Reference Codon Position')
        chain_ix = anot_header.index('chain')
        pdb_res_ix = anot_header.index('residue')

        # read in data for significant lines
        skip_ct = 0
        for l in myreader:
            # skip lines that don't have correct annotation
            if len(l) < res_ix:
                skip_ct += 1
                continue

            # add the residue information if it is significant
            res_info = (l[gene_ix], l[tx_ix], int(l[res_ix]))
            if res_info in signif_res:
                data.append(l)
            #gene_info = l[gene_ix]
            #if gene_info in signif_res:
                #data.append(l)
        logger.info('Skipped {0} lines'.format(skip_ct))

    # record the position of the columns in the header
    column_dict = {
        'pdb': pdb_ix,
        'chain': chain_ix,
        'pdb_res': pdb_res_ix,
        'gene' : gene_ix,
        'tx': tx_ix,
        'res': res_ix
    }

    return data, column_dict


def update_graph(gene2graph, cog, signif_struct_info, struct, radius):
    """Updates the residue neighbor graph based on the current structure.

    Residues are linked by edges if they are within the provided radius
    and are on the same gene.

    Parameters
    ----------
    gene2graph : dict
        dictionary with genes as keys pointing to significant hotspot residue
        neighbor graph
    signif_struct_info : dict
        identifies which residues are significant hotspots
    struct : Bio.PDB structure
        structure under consideration when populating the graph
    radius : float
        radius deemed close enough to add an edge between two residues

    Returns
    -------
    gene2graph : dict
        updated graph based on the provided structure
    """
    # get which residues are significant
    signif_pdb_pos = signif_struct_info.keys()
    possible_res = set(signif_pdb_pos)

    # find neighbor residues
    cog = {k: cog[k] for k in cog
           if (k[2], k[3][1]) in signif_pdb_pos}
    neighbors = pstruct.find_neighbors(cog, radius)

    #struct_info = struct_chain[pdb_id]

    # add edge if residues are neighbors
    avail_models = [m.id for m in struct]
    for s in signif_pdb_pos:
        tmp_chain, tmp_res = s
        cur_res = signif_struct_info[s]
        cur_gene = cur_res[0]

        # update gene2graph
        gene2graph.setdefault(cur_gene, {})
        gene2graph[cur_gene].setdefault(cur_res, set())

        for m in avail_models:

            try:
                # get neighbors
                tmp_id = struct[m][tmp_chain][int(tmp_res)].get_full_id()
                tmp_neighbors = set([(n[2], n[3][1]) for n in neighbors[tmp_id]])

                # get only neighbors that are significant and in the
                # same gene
                signif_neighbors = set([signif_struct_info[o]
                                        for o in (tmp_neighbors & possible_res)])
                signif_neighbors_gene = set([s for s in signif_neighbors
                                             if s[0] == cur_gene])

                # add result to the graph
                gene2graph[cur_gene][cur_res] = gene2graph[cur_gene][cur_res] | signif_neighbors_gene
            except KeyError:
                # skip deleted chains, or models without a chain
                # be careful this catches all keyerrors
                pass
    return gene2graph


def retrieve_components(graph_dict, tumor_type):
    """Get the connected components and format the output."""
    ttype_output = []
    for mygene in graph_dict:
        g = graph_dict[mygene]
        components = graph.connected_components(g)
        tmp = [mygene, tumor_type]
        for component in components:
            format_str = ';'.join('{0}:{1}'.format(n[1], n[2]) for n in component)
            tmp.append(format_str)
        ttype_output.append(tmp)
    return ttype_output


def main(opts):
    # read in the PDB info file
    pdb_info = utils.read_pdb_info(opts['pdb_info'])

    # read in multiple testing file
    mtc = read_delim(opts['multiple_testing'])
    header = mtc.pop(0)
    ttype_ix = header.index('Tumor Type')
    qval_ix = header.index('q-value')
    gene_ix = header.index('HUGO Symbol')
    tx_ix = header.index('Sequence Ontology Transcript')
    res_ix = header.index('CRAVAT Res')
    #mtc.sort(key=lambda x: x[0])

    # iterate through each tumor type
    output = []
    gene2graph_all = {}  # graphs for combined tumor types
    uniq_ttypes = set(m[ttype_ix] for m in mtc)
    for ttype in uniq_ttypes:
        logger.info('Working on {0} . . .'.format(ttype))
        # initialize the graph to empty
        gene2graph = {}  # graph for an individual tumor type

        # get the significant residues for the tumor type
        mtc_ttype = [m for m in mtc
                     if (m[ttype_ix] == ttype) and (float(m[qval_ix])<=opts['q_value'])]
        significant_res = set([(m[gene_ix], m[tx_ix], int(m[res_ix]))
                               for m in mtc_ttype])

        # read annotation file
        # HACK annotation_file = os.path.join(opts['annotation_dir'], 'mupit_mutations_' + ttype)
        annotation_file = opts['annotation_dir']

        annotation, col_pos = read_mupit_file(annotation_file, significant_res)
        pdb_ix = col_pos['pdb']
        anot_gene_ix = col_pos['gene']
        anot_tx_ix = col_pos['tx']
        anot_res_ix = col_pos['res']

        # sort by structure
        annotation.sort(key=lambda x: x[pdb_ix])

        for pdb_id, grp in it.groupby(annotation, lambda x: x[pdb_ix]):
            if pdb_id not in pdb_info:
                continue
            struct_info = pdb_info[pdb_id].copy()
            pdb_path = struct_info.pop('path')
            struct_chains = []
            for d in struct_info:
                struct_chains.extend(struct_info[d])
            #pdb_path = pdb2path[pdb_id]

            struct = utils.read_structure(pdb_path, pdb_id)
            if struct is None:
                continue  # skip if pdb file not found

            # calculate the centers of geometry
            cog = pstruct.calc_center_of_geometry(struct, struct_chains)

            # contains relevant mupit annotations for this pdb
            tmp = list(grp)

            # get significant residues
            signif_struct_info = {}
            for s in tmp:
                try:
                   tmp_pos =  (s[col_pos['chain']], int(s[col_pos['pdb_res']]))
                except:
                    print 'int error'
                    continue
                signif_struct_info[tmp_pos] = (s[anot_gene_ix], s[anot_tx_ix], s[anot_res_ix])

            # update the graph to reflect info from the current structure
            gene2graph = update_graph(gene2graph, cog, signif_struct_info,
                                      struct, opts['radius'])
            # update graph for the combined cross-tumor type regions
            banned_ttypes = ['COAD', 'READ', 'PANCAN12', 'CHOL', 'SARC',
                             'TGCT', 'THYM', 'UVM']
            if ttype not in banned_ttypes:
                gene2graph_all = update_graph(gene2graph_all, cog, signif_struct_info,
                                              struct, opts['radius'])

        # format the results into the output list
        tmp_out = retrieve_components(gene2graph, ttype)
        output += tmp_out
        logger.info('Finished {0}'.format(ttype))

    # update output to contain cross-tumor type reference regions
    tmp_out = retrieve_components(gene2graph_all, 'REF')
    output += tmp_out

    # write output
    with open(opts['output'], 'wb') as handle:
        for line in output:
            handle.write('\t'.join(line)+'\n')

    logger.info('Finished Successfully!!!')



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
