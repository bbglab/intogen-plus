# Import modules

import logging

import pandas as pd
import numpy as np

from intogen_combination.config import CONF

logger = logging.getLogger(__name__)


class Parser:

    def __init__(self, method, gene_coordinates):
        """Initialize an instance of the Parser class
        :param method: str, method to parse
        :param gene_coordinates: coordinates of CDS
        :return: None
        """
        self.name = method
        self.gene_coordinates = gene_coordinates
        self.gene_id = CONF[method]['GENE_ID']
        self.pvalue = CONF[method]['PVALUE']
        self.qvalue = CONF[method]['QVALUE']

    def read(self, input_file):
        """Read an output file. 
        :param input_file: path, file to read with results
        :return: pandas.DataFrame
        """
        try:
            df = pd.read_csv(input_file, header=0, sep="\t")
        except OSError as e:
            logger.warning('File {} not found'.format(input_file))
            return None
        # P-value
        try:
            df = df[np.isfinite(df[self.pvalue])]
        except Exception:
            logger.warning('No finite p-value'.format(input_file))
        df = df[[self.gene_id, self.pvalue, self.qvalue]]
        df.rename(columns={
            self.gene_id: 'GENE_ID', self.pvalue: 'PVALUE', self.qvalue: 'QVALUE'
        }, inplace=True)

        return df
