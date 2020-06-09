import csv
import gzip
import pickle

import click
import pandas as pd

from intogen_combination import parser, grid_optimizer, schulze, \
    stouffer_script, create_tiers_drivers
from intogen_combination.config import METHODS


def get_value(dict_values, key, position=0):
    if key in dict_values:
        return dict_values[key][position]  # Position = 0 for pvalue, Position = 1 for Qvalue
    else:
        return None


def main(outprefix='f', **files):

    # STEP 1
    ranking, pvalues = parser.parse(**files)

    with gzip.open("{}.step1".format(outprefix), "wb") as fd:
        pickle.dump(ranking, fd)

    with gzip.open("{}.step1b".format(outprefix), "wt") as fd:
        writer = csv.writer(fd, delimiter='\t')
        writer.writerow(
            ['SYMBOL'] + ["PVALUE_{}".format(h) for h in METHODS] + ["QVALUE_{}".format(h) for h in METHODS])
        for gene, values in pvalues.items():
            writer.writerow(
                [gene] + [get_value(values, m, 0) for m in METHODS] + [get_value(values, m, 1) for m in METHODS])

    # STEP 2
    optimized = grid_optimizer.run(ranking.copy(), **files)

    optimized.to_csv("{}.step2".format(outprefix), sep="\t", index=False, compression="gzip")

    # STEP 3
    df, ranking1 = schulze.optimal_ranking(optimized, ranking.copy(), borda=True)

    df.sort_values("RANKING").to_csv("{}.step3".format(outprefix), sep="\t", index=False, compression="gzip")

    with gzip.open("{}.step3b".format(outprefix), "wb") as fd:
        pickle.dump(ranking1, fd)

    # STEP 4
    data = pd.read_csv("{}.step1b".format(outprefix), sep='\t', compression="gzip")
    df = stouffer_script.run(data, df, optimized,
                        files['oncodrivefml'],
                        files['dndscv'],
                        brown=True, fisher=True)

    df.to_csv("{}.stouffer.out.gz".format(outprefix), sep='\t', index=False, compression="gzip")

    # STEP 5
    df_tiers = create_tiers_drivers.run(df.copy(), threshold=0.01, column_filter='QVALUE_stouffer_w')

    df_tiers.to_csv("{}.01.out.gz".format(outprefix), sep="\t", index=False, compression="gzip")

    # STEP 6
    df_tiers2 = create_tiers_drivers.run(df.copy(), threshold=0.05, column_filter='QVALUE_stouffer_w')

    df_tiers2.to_csv("{}.05.out.gz".format(outprefix), sep="\t", index=False, compression="gzip")


@click.command()
@click.option('--oncodriveclustl')
@click.option('--dndscv')
@click.option('--oncodrivefml')
@click.option('--hotmaps')
@click.option('--smregions')
@click.option('--cbase')
@click.option('--mutpanning')
@click.option('-o', '--output', required=True)
def cli(output, **kwargs):
    main(output, **{k: v for k, v in kwargs.items() if v is not None})


if __name__ == '__main__':
    cli()
