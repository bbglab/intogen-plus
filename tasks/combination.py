import os
import sys
import signal
import subprocess

from os import path
from .base import Task


class CombinationTask(Task):

    KEY = 'combination'

    def __init__(self, output_folder):

        super().__init__(output_folder)

        self.name = None
        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, self.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):
        pass

    def input_write(self, _, value):
        pass

    def input_end(self):
        pass

    def run(self):

        # Run Schulze
        if not path.exists(self.out_file):
            cmd = "bash -c 'export SCHULZE_DATA={0} && {1}/schulze.sh {2} {3}'".format(
                os.path.join(os.environ['INTOGEN_DATASETS'], 'combination'),
                os.path.join(os.environ['INTOGEN_METHODS'], 'combination'),
                self.output_folder, self.name)

            with subprocess.Popen(cmd, shell=True, stdin=sys.stdin, stderr=sys.stderr) as p:
                with open(self.out_file + ".pid", "wt") as fd:
                    fd.write("{}\n".format(p.pid))
                try:
                    errcode = p.wait()
                except:
                    p.kill()
                    p.wait()
                    raise
            os.unlink(self.out_file + ".pid")

            if errcode != 0:
                raise RuntimeError("{} [error] - code {}".format(self, errcode))

        return self.out_file

    def clean(self):

        if path.exists(self.out_file + ".pid"):
            with open(self.out_file + ".pid", "rt") as fd:
                pid = int(fd.readline().strip())
                os.kill(pid, signal.SIGTERM)
            os.unlink(self.out_file + ".pid")

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)
