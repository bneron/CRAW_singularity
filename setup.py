# -*- coding: utf-8 -*-
import sys
if sys.version_info[0] == 2:
    sys.exit("Sorry, Python 2 is not supported")

from distutils.core import setup
import time

setup(
    name="craw",
    version='{}-dev'.format(time.strftime('%Y%m%d')),
    author='Bertrand Néron',
    author_email='bneron@pasteur.fr',
    url="https://gitlab.pasteur.fr/bneron/craw",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    description="Compute the coverage at specified regions",
    long_description="""From a bam file and a file containing regions,
coverage compute the coverage for each positions of these regions in sense and antisense.""",
    platforms=["Unix"],
    install_requires=['pysam >= 0.9.1.4'],
    packages=['craw'],
    scripts=['bin/coverage'],
    data_files=[
        ("share/craw/doc/html", ['doc/build/html/index.html'])
    ]
)