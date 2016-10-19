# -*- coding: utf-8 -*-

import sys
if sys.version_info[0] == 2:
    sys.exit("Sorry, Python 2 is not supported")

import os
import time

from distutils import log
from distutils.core import setup
from distutils.core import Command
from distutils.dist import Distribution
from distutils.command.install import install as _install
from distutils.command.install_lib import install_lib as _install_lib
from distutils.command.install_data import install_data as _install_data
from distutils.command.install_scripts import install_scripts as _install_scripts
from distutils.errors import DistutilsFileError, DistutilsOptionError, DistutilsPlatformError
from distutils.util import subst_vars as distutils_subst_vars
from distutils.util import change_root, convert_path, get_platform


class test(Command):

    description = "run the unit tests against the build library"

    user_options = [('verbosity', 'v', 'verbosity of outputs (cumulative option ex -vv)', 1),
                    ('build-base=', 'b', "base build directory (default: 'build.build-base')"),
                    ('build-lib=', None, "build directory for all modules (default: 'build.build-lib')"),
                    ('plat-name=', 'p', "platform name to build for, if supported (default: {0})".format(get_platform())),
                    ]

    help_options = []

    def initialize_options(self):
        self.verbosity = None
        self.build_base = 'build'
        self.build_lib = None
        self.build_purelib = None
        self.build_platlib = None
        self.plat_name = None
        self.skip_build = 0
        self.warn_dir = 1

    def finalize_options(self):
        if self.verbosity is None:
            self.verbosity = 0
        else:
            self.verbosity = int(self.verbosity)

        if self.plat_name is None:
            self.plat_name = get_platform()
        else:
            # plat-name only supported for windows (other platforms are
            # supported via ./configure flags, if at all).  Avoid misleading
            # other platforms.
            if os.name != 'nt':
                raise DistutilsOptionError(
                            "--plat-name only supported on Windows (try "
                            "using './configure --help' on your platform)")

        plat_specifier = ".{0}-{1}".format(self.plat_name, sys.version[0:3])

        # Make it so Python 2.x and Python 2.x with --with-pydebug don't
        # share the same build directories. Doing so confuses the build
        # process for C modules
        if hasattr(sys, 'gettotalrefcount'):
            plat_specifier += '-pydebug'

        # 'build_purelib' and 'build_platlib' just default to 'lib' and
        # 'lib.<plat>' under the base build directory.  We only use one of
        # them for a given distribution, though --
        if self.build_purelib is None:
            self.build_purelib = os.path.join(self.build_base, 'lib')
        if self.build_platlib is None:
            self.build_platlib = os.path.join(self.build_base,
                                              'lib' + plat_specifier)

        # 'build_lib' is the actual directory that we will use for this
        # particular module distribution -- if user didn't supply it, pick
        # one of 'build_purelib' or 'build_platlib'.
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib


    def run(self):
        """
        """
        if not self.skip_build:
            self.run_command('build')
            # If we built for any other platform, we can't install.
            build_plat = self.distribution.get_command_obj('build').plat_name
            # check warn_dir - it is a clue that the 'install' is happening
            # internally, and not to sys.path, so we don't check the platform
            # matches what we are running.
            if self.warn_dir and build_plat != get_platform():
                raise DistutilsPlatformError("Can't test when "
                                             "cross-compiling")
        from tests import run_tests
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib

        log.info("running test")
        os.environ['CRAW_HOME'] = os.path.dirname(os.path.abspath(__file__))
        kind_of_skipped = {}
        test_res = run_tests.run_functional_tests(False, verbosity=self.verbosity)
        for test in test_res.skipped:
            kind_of_skipped[test[1]] = True
        for skip_reason in kind_of_skipped.keys():
            msg = skip_reason
            log.warn(msg)
        if not test_res.wasSuccessful():
            sys.exit("some tests fails. Run python setup.py test -vv to have more details")



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
                            'PREFIXCONF': os.path.join(get_install_conf_dir(inst), 'craw'),
                            'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'craw'),
                            'PREFIXDOC': os.path.join(get_install_doc_dir(inst), 'craw')
                            }
            for _file in self.distribution.fix_scripts:
                subst_file(_file, vars_2_subst)
        _install_scripts.run(self)


class install_data(_install_data):

    user_options = [
        ('install-dir=', 'd',
         "base directory for installing data files "
         "(default: installation base dir)"),
        ('root=', None,
         "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ]

    boolean_options = ['force']

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        self.install_dir = get_install_data_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.files_2_install = self.distribution.data_files
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_data = {}\n'.format(self.install_dir))

    def run(self):
        self.mkpath(self.install_dir)
        for f in self.files_2_install:
            if isinstance(f, str):
                if not os.path.exists(f):
                    log.warn("WARNING the document {} cannot be found, installation skipped".format(f))
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'{0}' -- installing right in '{1}'".format((f, self.install_dir)))
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of path
                dir_ = convert_path(f[0])
                if not os.path.isabs(dir_):
                    dir_ = os.path.join(self.install_dir, dir_)
                elif self.root:
                    dir_ = change_root(self.root, dir_)
                self.mkpath(dir_)
                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(dir_)
                else:
                    # Copy files, adding them to the list of output files.
                    for data in f[1]:
                        data = convert_path(data)  # return name that will work on the native filesystem
                        if not os.path.exists(data):
                            log.warn("WARNING the document {} cannot be found, installation skipped".format(data))
                            continue
                        if os.path.isdir(data):
                            out = self.copy_tree(data, dir_)
                            self.outfiles.extend(out)
                        else:
                            (out, _) = self.copy_file(data, dir_)
                            self.outfiles.append(out)


class install_conf(_install_data):

    _install.sub_commands += [('install_conf', lambda self:True)]
    description = "installation directory for configuration files"
    setattr(_install, 'install_conf', None)
    _install.user_options.append(('install-conf=', None, description))

    user_options = [
        ('install-conf=', 'd',
         "base directory for installing configuration files "
         "(default: installation base dir etc)"),
        ('root=', None,
         "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ]

    boolean_options = ['force']

    def initialize_options(self):
        _install_data.initialize_options(self)
        self.conf_files = self.distribution.conf_files


    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        self.install_dir = get_install_conf_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )

    def run(self):
        if not self.conf_files:
            return
        self.mkpath(self.install_dir)
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        vars_2_subst = {
                        'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'craw'),
                        }
        for f in self.conf_files:
            if isinstance(f, str):
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'{0}' -- installing right in '{1}'".format(f, self.install_dir))
                dest = os.path.join(self.install_dir, f + ".new")
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of files
                _dir = convert_path(f[0])
                if not os.path.isabs(_dir):
                    _dir = os.path.join(self.install_dir, _dir)
                elif self.root:
                    _dir = change_root(self.root, _dir)
                self.mkpath(_dir)

                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(_dir)
                else:
                    # Copy files, adding them to the list of output files.
                    for conf in f[1]:
                        conf = convert_path(conf)
                        dest = os.path.join(_dir, os.path.basename(conf) + ".new")
                        (out, _) = self.copy_file(conf, dest)
                        if conf in self.distribution.fix_conf:
                            input_file = out
                            output_file = input_file + '.tmp'
                            subst_vars(input_file, output_file, vars_2_subst)
                            if os.path.exists(input_file):
                                os.unlink(input_file)
                            self.move_file(output_file, input_file)
                            self.outfiles.append(input_file)


class install_doc(install_data):

    _install.sub_commands += [('install_doc', lambda self: not self.no_doc)]

    description = "installation directory for documentation files"

    setattr(_install, 'install_doc', None)
    setattr(_install, 'no_doc', None)

    _install.user_options.append(('install-doc=', None, description))
    _install.user_options.append(('no-doc', None, 'do not install documentation'))

    user_options = [
        ('install-doc=', 'd', "base directory for installing documentation files "
                              "(default: installation base dir share/doc)"),
        ('root=', None, "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ('no-doc', None, 'do not install documentation')
        ]

    boolean_options = ['force']

    def initialize_options(self):
        install_data.initialize_options(self)
        self.install_doc = None
        self.no_doc = None
        self.files_2_install = self.distribution.doc_files

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        self.install_dir = get_install_doc_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.no_doc = inst.get('no_doc', ('command line', False))[1]

    def run(self):
        install_data.run(self)




class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        # It's important to define options before to call __init__
        # otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        self.conf_files = None
        self.doc_files = None
        self.fix_lib = None
        self.fix_conf = None
        self.fix_scripts = None
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py test       will run the tests on the newly build library
  setup.py install    will install the package
"""


def get_install_data_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join('/', 'usr', 'share')
    return install_dir


def get_install_conf_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_conf' in inst:
        install_dir = inst['install_conf'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'etc')
    else:
        install_dir = '/etc'
    return install_dir


def get_install_doc_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_doc' in inst:
        install_dir = inst['install_doc'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share', 'craw', 'doc')
    else:
        install_dir = os.path.join('/', 'usr', 'share', 'craw', 'doc')
    return install_dir


def subst_vars(src, dst, vars):
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


require_python = ['python (>3.1)']
require_packages = ['pysam (>=0.9.1)']

setup(name="craw",
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

      #data_files=[],
      #conf_files=[],
      doc_files=[('html', ['doc/build/html/']),
                 ('pdf', ['doc/build/latex/craw.pdf'])],
      # configuration files where some variables must be fix by install_conf
      #fix_conf=['conf/config.ini'],
      # library file where some variable must be fix by install_lib
      fix_lib=['craw/__init__.py'],
      # scripts file where some variable must be fix by install_scripts
      fix_scripts=['coverage'],
      cmdclass={'test': test,
               'install_lib':  install_lib,
               'install_scripts':  install_scripts,
               'install_data': install_data,
               'install_conf': install_conf,
               'install_doc': install_doc,
               },
      distclass=UsageDistribution
      )
