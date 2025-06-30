import argparse
import csv


def parse_arguments():
    info = 'Extracts the residues that were significant hotspots'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Full output from hotspot script which includes p-values')
    parser.add_argument('-s', '--significance-level',
                        type=float, required=True,
                        help='The significance level to call it a hotspot residue')
    parser.add_argument('-o', '--output',
                        type=str, default=None,
                        help='Output listing only the residues that were hotspots')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    output = [['Structure', 'Tumor Type', 'Model', 'Chain',
               'Residue Position', 'p-value']]
    with open(opts['input']) as handle:
        # get csv reader
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)

        # define column position from the header info
        struct_ix = header.index('Structure')
        ttype_ix = header.index('Tumor Type')
        model_ix = header.index('Model')
        chain_ix = header.index('Chain')
        res_ix = header.index('Mutation Residues')
        pval_ix = header.index('Hotspot P-value')

        # iterate over every structure/tumor type pair
        for line in myreader:
            # skip if no p-values
            if not line[pval_ix]:
                continue
            chains = line[chain_ix].split(',')
            models = line[model_ix].split(',')
            res_pos = line[res_ix].split(',')
            res_pval = map(float, line[pval_ix].split(','))

            # iterate over p-values
            for i, p in enumerate(res_pval):
                if p <= opts['significance_level']:
                    output.append([line[struct_ix], line[ttype_ix], models[i],
                                   chains[i], res_pos[i], p])

    # write output to file
    if opts['output'] is not None:
        with open(opts['output'], 'w') as handle:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(output)
    else:
        return output



if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
