import click
import json
import itertools

CB = dict(zip('ACGT','TGCA'))

def get_context():
    ## get the 96 trinucleotides list
    subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0] != z[1]]
    flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
    contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
    contexts_formatted = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]

    return contexts_formatted


def generate_complementary_triplet(sequence_sub):
    # Convert the sequence to uppercase for consistency
    sequence, sub = sequence_sub.split(">")
    sequence = sequence.upper()

    # Generate the complementary triplet
    complementary_triplet = ''.join(CB[nucleotide] for nucleotide in sequence)[::-1]
    subs = CB[sub]

    return complementary_triplet+'>'+subs


def run(infile, outfile):
    """
    Takes as input the output of bgsignature and prepares it for boostDM run
    """

    with open(infile, 'r') as bgsig_file:
        bgsig = json.load(bgsig_file)

    bgsig_dict=dict()

    for context in get_context():
        complementary_triplet = generate_complementary_triplet(context)
        
        bgsig_context = bgsig[context] if bgsig.get(context) != None else 0
        bgsig_complementary = bgsig[complementary_triplet] if bgsig.get(complementary_triplet) != None else 0

        bgsig_dict[context] = bgsig_context + bgsig_complementary


    with open(outfile + '.mutrate.json', 'w') as o:
        json.dump(bgsig_dict, o, indent=4)



@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(input, output):
    run(input, output)
    

if __name__ == "__main__":
    cli()