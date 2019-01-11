
from os import path

from setuptools import Extension
from setuptools import setup, find_packages
from simReg import __version__


directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

setup(
    name='simReg',
    version=__version__,
    packages=find_packages(),
    package_data={'simReg': ['*.txt.gz', '*.conf.template', '*.conf.template.spec', '*.pyx']},
    url="https://bitbucket.org/bbglab/simReg",
    download_url="https://bitbucket.org/bbglab/simReg/get/"+__version__+".tar.gz",
    license='UPF Free Source Code',
    author='BBGLab (Barcelona Biomedical Genomics Lab)',
    author_email='bbglab@irbbarcelona.org',
    description='Identify signals of positive selection in somatic mutations',
    setup_requires=[
        'cython',
    ],
    install_requires=required,
#    ext_modules=[Extension('simReg.walker_cython', ['simReg/walker_cython.pyx'])],
    entry_points={
        'console_scripts': [
            'simreg = simReg.main:cmdline'
        ]
    }
)
