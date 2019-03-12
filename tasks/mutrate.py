import os

from .base import Task, run_command


class MutrateTask(Task):

    KEY = 'mutrate'
    
    def run(self):

        weights = os.path.join(
            self.output_project,
            'deconstructsig', '{}.out.gz'.format(self.name)
        )

        annotmuts = self.in_file.replace(".out.gz", ".annotmuts.gz")
        genemuts = self.in_file.replace(".out.gz", ".genemuts.gz")
        # weights = self.weight_file
        cores = os.environ.get('INTOGEN_CPUS', 4)
        output = os.path.splitext(os.path.splitext(self.out_file)[0])[0]

        run_command(f"{self.cmdline} {annotmuts} {genemuts} {weights} {cores} {output}")