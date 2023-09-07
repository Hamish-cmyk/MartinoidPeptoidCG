from setuptools import setup

setup(
    name='Martinoid',
    version='0.0.1',
    description='This module was inspired by martinize (http://cgmartini.nl/index.php/tools2/proteins-and-bilayers/204-martinize) and has been created to perform automatic topology building of peptoids within the MARTINI forcefield (v2.1) in the GROMACS program.',
    url='https://github.com/avanteijlingen/ORCA-Parser',
    author='Alexander van Teijlingen',
    author_email='a.vant@linuxmail.org',
    license='BSD 2-clause',
    packages=['orca_parser'],
    install_requires=['ase',
                      'pandas',
                      'numpy',
                      'mdtraj',
                      'matplotlib'
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
    ],
)
