from os import path
from setuptools import setup, find_packages, Extension


VERSION = "0.1"
DESCRIPTION = "IntOGen Combination functionality"

directory = path.dirname(path.abspath(__file__))

# Get requirements from the requirements.txt file
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()


# Get the long description from the README file
with open(path.join(directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='intogen_combination',
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="",
    author="Barcelona Biomedical Genomics Lab",
    author_email="bbglab@irbbarcelona.org",
    license="Apache Software License 2.0",
    packages=find_packages(),
    package_data={'intogen_combination': ['*.cfg']},
    install_requires=required,
    setup_requires=[
        'cython==0.29.36',
    ],
    ext_modules=[Extension('intogen_combination.schulze_strongest_path_cython', ['intogen_combination/schulze_strongest_path_cython.pyx'])],
    entry_points={
            'console_scripts': [
                'intogen-combine = intogen_combination.main:cli'
            ]
        },
)
