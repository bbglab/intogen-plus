import csv
import gzip
from os import path


def main(input, output):
    with gzip.open(input) as fin, open(output, 'w') as fout:
        reader_in = csv.DictReader(fin, delimiter='\t')
        reader_out = csv.DictWriter(fout, delimiter='\t')

        for info in reader_in:
            id_ = info['PDBId']
            if id_.startswith('ENSP'):
                info['NonBioPDBPath'] = ''
                info['PDBPath'] = 'ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model/{}.pdb.gz'.format(id_)
            elif id_.startswith('NP'):
                info['NonBioPDBPath'] = ''
                info['PDBPath'] = 'ModBase_H_sapiens_2013_refseq/models/model/{}.pdb.gz'.format(id_)
            else:
                pdb_path = 'biounit/coordinates/all/{}.pdb1.gz'.format(id_)
                if path.exists(path.join(path.dirname(output), pdb_path)):
                    info['NonBioPDBPath'] = 'structures/all/pdb/pdb{}.ent.gz'.format(id_)
                    info['PDBPath'] = 'ModBase_H_sapiens_2013_refseq/models/model/{}.pdb.gz'.format(id_)
                else:
                    info['NonBioPDBPath'] = ''
                    info['PDBPath'] = ''

            reader_out.writerow(info)