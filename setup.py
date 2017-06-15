# -*- coding: utf-8 -*-

###########################################################################
#                                                                         #
# This file is part of Counter RNAseq Window (craw) package.              #
#                                                                         #
#    Authors: Bertrand Néron                                              #
#    Copyright © 2017  Institut Pasteur (Paris).                          #
#    see COPYRIGHT file for details.                                      #
#                                                                         #
#    craw is free software: you can redistribute it and/or modify         #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         #
#    craw is distributed in the hope that it will be useful,              #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                 #
#    See the GNU General Public License for more details.                 #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with craw (see COPYING file).                                  #
#    If not, see <http://www.gnu.org/licenses/>.                          #
#                                                                         #
###########################################################################


import sys
if sys.version_info[0] == 2:
    sys.exit("Sorry, Python 2 is not supported")

import os
import time
import sysconfig

from distutils.errors import DistutilsFileError
from distutils.util import subst_vars as distutils_subst_vars

from setuptools import setup
from setuptools.dist import Distribution
from setuptools.command.install_scripts import install_scripts as _install_scripts
from setuptools.command.install_lib import install_lib as _install_lib



class install_lib(_install_lib):

    def finalize_options(self):
        _install_lib.finalize_options(self)

    def run(self):
        def subst_file(_file, vars_2_subst):
            input_file = os.path.join(self.build_dir, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)

        inst = self.distribution.command_options.get('install')
        if inst:
            if self.distribution.fix_lib is not None:
                vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                                'VERSION': self.distribution.get_version()
                                }
                for _file in self.distribution.fix_lib:
                    subst_file(_file, vars_2_subst)
        _install_lib.run(self)


class install_scripts(_install_scripts):

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        _install_scripts.finalize_options(self)

    def run(self):
        def subst_file(_file, vars_2_subst):
            input_file = os.path.join(self.build_dir, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)

        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        if self.distribution.fix_scripts is not None:
            vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                            'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'craw'),
                            }
            for _file in self.distribution.fix_scripts:
                subst_file(_file, vars_2_subst)
        _install_scripts.run(self)


class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        # It's important to define options before to call __init__
        # otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        self.fix_lib = None
        self.fix_scripts = None
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py install    will install the package
"""


def get_install_data_dir(inst):
    """
    :param inst: installation option
    :type inst: dict
    :return: the prefix where to install data
    :rtype: string
    """

    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])
    elif 'user' in inst:
        import site
        inst['prefix'] = ('command line', site.USER_BASE)
    elif 'root' in inst:
        inst['prefix'] = ('command line',
                          os.path.join(inst['root'][1],
                                       sysconfig.get_path('data').strip(os.path.sep)
                                       )
                          )

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join(sysconfig.get_path('data'), 'share')
    return install_dir


def subst_vars(src, dst, vars):
    """
    substitute variables (string starting with $) in file
    :param src: the file containing variable to substitute
    :type src: string
    :param dst: the destination file
    :type dst: string
    :param vars: the variables to substitute in dict key are variable name
    :type vars: dict
    """
    try:
        src_file = open(src, "r")
    except os.error as err:
        raise DistutilsFileError("could not open '{0}': {1}".format(src, err))
    try:
        dest_file = open(dst, "w")
    except os.error as err:
        raise DistutilsFileError("could not create '{0}': {1}".format(dst, err))
    with src_file, dest_file:
        for line in src_file:
            new_line = distutils_subst_vars(line, vars)
            dest_file.write(new_line)


def expand_data(data_to_expand):
    """
    From data structure like setup.py data_files (see http://)
     [(directory/where/to/copy/the/file, [path/to/file/in/archive/file1, ...]), ...]
    but instead of the original struct this one accept to specify a directory in elements to copy.

    This function will generate one entry for all *content* of the directory and subdirectory
    recursively, to in fine copy the tree in archive in dest on the host

    the first level of directory itself is not include (which allow to rename it)
    :param data_to_expand:
    :type  data_to_expand: list of tuple
    :return: list of tuple
    """
    def remove_prefix(prefix, path):
        prefix = os.path.normpath(prefix)
        path = os.path.normpath(path)
        to_remove = len([i for i in prefix.split(os.path.sep) if i])
        truncated = [i for i in path.split(os.path.sep) if i][to_remove:]
        truncated = os.path.sep.join(truncated)
        return truncated

    data_struct = []
    for base_dest_dir, src in data_to_expand:
        base_dest_dir = os.path.normpath(base_dest_dir)
        for one_src in src:
            if os.path.isdir(one_src):
                for path, _, files in os.walk(one_src):
                    if not files:
                        continue
                    path_2_create = remove_prefix(one_src, path)
                    data_struct.append(
                        (os.path.join(base_dest_dir, path_2_create), [os.path.join(path, f) for f in files]))
            if os.path.isfile(one_src):
                data_struct.append((base_dest_dir, [one_src]))
    return data_struct

try:
    from pypandoc import convert
    def read_md(f):
        return convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, "
          "could not convert Markdown to RST")

    def read_md(f): return open(f, 'r').read()

setup(name="craw",
      version='0.9.0rc4',
      author='Bertrand Néron',
      author_email='bneron@pasteur.fr',
      url="https://gitlab.pasteur.fr/bneron/craw",
      keywords=['bioinformatics', 'RNAseq', 'coverage'],
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
      description="Counter RNA seq Window is a package which aim to compute and visualize the coverage of RNA seq experiment.",
      # gitlab and github use md format for README whereas pypi use restructuredtext
      # so to display correctly the readme on the pypi page we need
      # to convert md -> rst using pandoc/pypandoc
      long_description=read_md('README.md'),
      platforms=["Unix"],
      install_requires=['pysam>=0.9.1.4', 'matplotlib>=1.5.3', 'pandas>=0.17.1', 'numpy>=1.11.2', 'pillow>=3.4.2',
                        'psutil>=4.0.0'],
      packages=['craw'],
      scripts=['bin/craw_coverage', 'bin/craw_htmp'],

      data_files=expand_data([('share/craw/doc/html', ['doc/build/html/']),
                              ('share/craw/doc/pdf/', ['doc/build/latex/CounterRNAseqWindow.pdf'])]),

      # library file where some variable must be fix by install_lib
      fix_lib=['craw/__init__.py'],
      # scripts file where some variable must be fix by install_scripts
      fix_scripts=['craw_coverage'],

      cmdclass={'install_lib':  install_lib,
                'install_scripts':  install_scripts,
               },
      distclass=UsageDistribution
      )
