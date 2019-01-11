
from os import path

from setuptools import setup, find_packages
from smdeg import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

setup(
    name='smregions',
    version=__version__,
    packages=find_packages(),
    package_data={'smregions': ['*.conf.template', '*.conf.template.spec']},
    # url="https://bitbucket.org/bbglab/simReg",
    # download_url="https://bitbucket.org/bbglab/simReg/get/"+__version__+".tar.gz",
    license='UPF Free Source Code',
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description='Identify signals of positive selection in somatic mutations',
    install_requires=required,
    entry_points={
        'console_scripts': [
            'smregions = smregions.main:cmdline'
        ]
    }
)
