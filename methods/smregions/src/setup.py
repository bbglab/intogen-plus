
from os import path

from setuptools import setup, find_packages
from smregions import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

setup(
    name='smregions',
    version=__version__,
    packages=find_packages(),
    package_data={'smregions': ['*.conf.template', '*.conf.template.spec']},
    license='',
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description='',
    install_requires=required,
    entry_points={
        'console_scripts': [
            'smregions = smregions.main:cmdline'
        ]
    }
)
