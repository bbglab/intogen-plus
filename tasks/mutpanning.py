import os
import gzip

from pyliftover import LiftOver

from .base import Task, valid_consequence, run_command


class MutPanningTask(Task):

    KEY = 'mutpanning'
    INPUT_HEADER = [
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
    ]
    SAMPLES_HEADER = [
        'ID', 'Sample', 'Cohort',  # 'Subtype', 'ConfidenceLevel', 'Study'
    ]
    OUTPUT_HEADER = [
        'Name', 'TargetSize', 'TargetSizeSyn', 'Count',
        'CountSyn', 'Significance', 'FDR'
    ]
    LIFTOVER = None

    def input_start(self):
        super().input_start()

        genome = os.environ['INTOGEN_GENOME'].lower()
        if genome != "hg19":
            self.LIFTOVER = LiftOver(os.path.join(os.environ['INTOGEN_DATASETS'], 'preprocess', '{}ToHg19.over.chain.gz'.format(genome)))

    def input_write(self, _, value):
        if value['Feature'] in self.transcripts:
            consequence = value['Consequence'].split(',')[0].replace('missense_variant', 'Missense_Mutation')
            chromosome, position = value['Location'].split(':')
            _, sample, ref, alt = value['#Uploaded_variation'].split('__')
            variant_type = None
            if ref == '-':
                variant_type == "INS"
            elif alt == '-':
                variant_type == "DEL"
            elif len(ref) == len(alt) and len(ref) == 1:
                variant_type == "SNP"
            elif len(ref) == len(alt) and len(ref) == 2:
                variant_type == "DNP"
            elif len(ref) == len(alt) and len(ref) == 2:
                variant_type == "TNP"
            else:
                return

        if self.LIFTOVER:
            strand = '-' if value['STRAND'] == '-1' else '+'
            hg19_position = self.LIFTOVER.convert_coordinate("chr{}".format(chromosome), int(position) - 1, strand)
            if hg19_position is None or len(hg19_position) != 1:
                return
            position = str(hg19_position[0][1] + 1)

        self.write_row([
            value['SYMBOL'],      # Hugo_Symbol
            chromosome,           # Chromosome
            position,             # Start_Position in Hg19
            position,             # End_Position
            value['STRAND'],      # Strand
            consequence,          # Variant_Classification,
            variant_type,         # Missense_Mutation
            ref,                  # Reference_Allele
            ref,                  # Tumor_Seq_Allele1
            alt,                  # Tumor_Seq_Allele2
            sample,               # Tumor_Sample_Barcode
        ])

    def create_input_files(self, mutations_file, sample_file):
        samples = set([])
        with gzip.open(self.in_file, 'rb') as fd, open(mutations_file, 'w') as out:
            for line in fd:
                out.write(line.decode())
                samples.add(line.decode().strip().split('\t')[-1])

        with open(sample_file, 'w') as fd:
            fd.write('\t'.join(self.SAMPLES_HEADER) + '\n')
            for _id, sample in enumerate(samples):
                fd.write(f"{_id}\t{sample}\t{self.name}\n")

    def create_output_files(self, result_file):
        """
        Create an empty output file, with only the header
        :param result_file:
        :return:
        """
        os.makedirs(os.path.dirname(result_file), exist_ok=True)
        with open(result_file, 'w') as fd:
            fd.write('\t'.join(self.OUTPUT_HEADER) + '\n')

    def run(self):
        # Input parameters
        mutations_file = os.path.join(self.output_work, f'{self.name}.mutations')
        sample_file = os.path.join(self.output_work, f'{self.name}.samples')
        self.create_input_files(mutations_file, sample_file)
        database = os.path.join(os.environ.get("INTOGEN_DATASETS"), 'mutpanning', 'Hg19/')

        # Output parameters
        result_file = os.path.join(
            self.output_folder,
            'SignificanceFiltered',
            f'Significance{self.name}.txt'
        )

        # Create an output file with only the header. If the method fails,
        # this is the file that will be compressed as self.out_file
        self.create_output_files(result_file)

        # Run MutPanning
        run_command(f"""
        {self.cmdline} {self.output_folder} {mutations_file} {sample_file} {database} && \
        cat {result_file} | gzip > {self.out_file}
        """)
