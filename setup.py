from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='Martinoid',
    version='0.0.5',
    description='This module was inspired by martinize (http://cgmartini.nl/index.php/tools2/proteins-and-bilayers/204-martinize) and has been created to perform automatic topology building of peptoids within the MARTINI forcefield (v2.1) in the GROMACS program.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Hamish-cmyk/MartinoidPeptoidCG',
    author='Hamish W. A. Swanson, Alexander van Teijlingen',
    author_email='a.vant@linuxmail.org',
    license='BSD 2-clause',
    packages=['Martinoid'],
    install_requires=['ase',
                      'pandas',
                      'numpy',
                      'mdtraj',
                      'matplotlib',
                      'argparse'
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
    ],
    include_package_data=True
)
