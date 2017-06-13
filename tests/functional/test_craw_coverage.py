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

import shutil
import tempfile
import os
from subprocess import Popen, PIPE
from itertools import zip_longest

from tests import CRAWTest, which


class Test(CRAWTest):

    def setUp(self):
        if 'CRAW_HOME' in os.environ:
            self.craw_home = os.environ['CRAW_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.craw_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))
        self.tmp_dir = tempfile.gettempdir()
        self.bin = os.path.join(self.craw_home, 'bin', 'craw_coverage') if self.local_install else which('craw_coverage')


    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
            pass
        except:
            pass


    def test_bam_with_fixed_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'small.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --bam={bam_file} --annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 bam_file=os.path.join(self._data_dir, 'small.bam'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_wo_start.txt'),
                                                 ref_col='Position',
                                                 before=5,
                                                 after=3,
                                                 qual=0,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_bam_with_chr_strand_col(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'coverage_fix_window_chr_strand_col.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --bam={bam_file} --annot={annot_file} " \
                  "--chr-col={chr_col} " \
                  "--strand-col={strand_col} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 bam_file=os.path.join(self._data_dir, 'small.bam'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_wo_start_chr_strand_col.txt'),
                                                 ref_col='Position',
                                                 chr_col='chr',
                                                 strand_col='brin',
                                                 before=5,
                                                 after=3,
                                                 qual=0,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_bam_with_var_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'coverage_var_window.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --bam={bam_file} --annot={annot_file} " \
                  "--ref-col={ref_col} " \
                  "--start-col={start_col} " \
                  "--stop-col={stop_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 bam_file=os.path.join(self._data_dir, 'small.bam'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_w_start.txt'),
                                                 ref_col='Position',
                                                 start_col='beg',
                                                 stop_col='end',
                                                 qual=15,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_wig_with_fixed_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_fixed_window.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig={wig_file} --annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 wig_file=os.path.join(self._data_dir, 'small_fixed.wig'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
                                                 ref_col='Position',
                                                 before=5,
                                                 after=3,
                                                 qual=0,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_wig_with_var_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_var_window.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig={wig_file} --annot={annot_file} " \
                  "--ref-col={ref_col} " \
                  "--start-col={start_col} " \
                  "--stop-col={stop_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 wig_file=os.path.join(self._data_dir, 'small_variable.wig'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_4_wig_var_win.txt'),
                                                 ref_col='Position',
                                                 start_col='beg',
                                                 stop_col='end',
                                                 qual=15,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_2wig_with_fixed_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_splited_fixed_window.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig-for={wig_for} --wig-rev={wig_rev} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 wig_for=os.path.join(self._data_dir, 'small_fixed.wig'),
                                                 wig_rev=os.path.join(self._data_dir, 'small_fixed_reverse.wig'),
                                                 annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
                                                 ref_col='Position',
                                                 before=5,
                                                 after=3,
                                                 qual=0,
                                                 out_file=test_result_path
                                                )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_only_forward_wig(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_only_forward.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig-for={wig_for} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
            bin=self.bin,
            wig_for=os.path.join(self._data_dir, 'small_fixed.wig'),
            annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
            ref_col='Position',
            before=5,
            after=3,
            qual=0,
            out_file=test_result_path
        )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                             cov_process.returncode,
                             command,
                             ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                         ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_only_reverse_wig(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_only_reverse.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig-rev={wig_rev} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
            bin=self.bin,
            wig_rev=os.path.join(self._data_dir, 'small_fixed_reverse.wig'),
            annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
            ref_col='Position',
            before=5,
            after=3,
            qual=0,
            out_file=test_result_path
        )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                             cov_process.returncode,
                             command,
                             ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                         ))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()

        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self._check_coverage_file(expected_result, test_result)


    def test_forward_n_mixed_wig(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_only_forward.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --wig={wig_file} " \
                  "--wig-for={wig_for} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
            bin=self.bin,
            wig_file=os.path.join(self._data_dir, 'small_variable.wig'),
            wig_for=os.path.join(self._data_dir, 'small_fixed.wig'),
            annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
            ref_col='Position',
            before=5,
            after=3,
            qual=0,
            out_file=test_result_path
        )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 1)
        last_line = cov_process.stderr.readlines()[-1].decode('utf-8')
        self.assertEqual(last_line,
                        "argparse.ArgumentError: '--wig' option cannot be specified in the same time as '--wig-for' or '--wig-rev' options.\n")


    def test_bam_n_wig(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_only_forward.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--bam={bam_file} " \
                  "--wig={wig_file} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
            bin=self.bin,
            bam_file=os.path.join(self._data_dir, 'small.bam'),
            wig_file=os.path.join(self._data_dir, 'small_variable.wig'),
            annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
            ref_col='Position',
            before=5,
            after=3,
            qual=0,
            out_file=test_result_path
        )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 1)
        last_line = cov_process.stderr.readlines()[-1].decode('utf-8')
        self.assertEqual(last_line,
                         "argparse.ArgumentError: '--bam' option cannot be specified in the same time as '--wig', '--wig-for' or '--wig-rev' options.\n")

    def test_no_input_file(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'wig_only_forward.cov'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--annot={annot_file} " \
                  "--before={before} " \
                  "--after={after} " \
                  "--ref-col={ref_col} " \
                  "--qual-thr={qual} " \
                  "--quiet " \
                  "--output={out_file} ".format(
            bin=self.bin,
            annot_file=os.path.join(self._data_dir, 'annotation_4_wig_fixed_win.txt'),
            ref_col='Position',
            before=5,
            after=3,
            qual=0,
            out_file=test_result_path
        )
        # print("\n@@@", command)
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=PIPE,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 1)
        last_line = cov_process.stderr.readlines()[-1].decode('utf-8')
        self.assertEqual(last_line,
                         "argparse.ArgumentError: At least one of these options must be specified '--bam', '--wig' , '--wig-for', '--wig-rev'.\n")

    def _check_coverage_file(self, expected_result, test_result):
        for expected, result in zip_longest(expected_result, test_result, fillvalue=''):
            if expected.startswith("# Version:"):
                continue
            if expected.startswith("# --annot="):
                continue
            elif expected.startswith("# --bam="):
                continue
            elif expected.startswith("# --wig="):
                continue
            elif expected.startswith("# --wig-for="):
                continue
            elif expected.startswith("# --wig-rev="):
                continue
            elif expected.startswith("# --output="):
                continue
            else:
                self.assertEqual(expected, result)

