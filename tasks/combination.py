import os
import sys
import signal
import subprocess

from os import path
from .base import Task, run_command


class CombinationTask(Task):

    KEY = 'combination_mutpanning'

    def run(self):
        run_command(f"{self.cmdline} {self.output_folder} {self.name}")