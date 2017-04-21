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
import hashlib
from PIL import Image

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
        self.bin = os.path.join(self.craw_home, 'bin', 'craw_htmp') if self.local_install else which('craw_htmp')


    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass


    def test_raw(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'htmp_raw_lin.png'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--size {size} " \
                  "--out={out_file} " \
                  "--quiet " \
                  "{cov_file}".format(
                              bin=self.bin,
                              size='raw',
                              out_file=test_result_path,
                              cov_file=os.path.join(self._data_dir, '4_htmp.cov')
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
            msg = "craw_htmp execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "craw_htmp finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))
        for sense in ('sense', 'antisense'):
            filename, suffix = os.path.splitext(output_filename)
            filename = "{}.{}{}".format(filename, sense, suffix)
            expected_result_path = os.path.join(self._data_dir, filename)
            expected_md5 = hashlib.md5()
            expected_md5.update(open(expected_result_path,'rb').read())

            result_path, suffix = os.path.splitext(test_result_path)
            result_path = "{}.{}{}".format(result_path, sense, suffix)
            result_md5 = hashlib.md5()
            result_md5.update(open(result_path,'rb').read())
            self.assertEqual(expected_md5.hexdigest(), result_md5.hexdigest())


    def test_raw_log(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'htmp_raw_log.png'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--size {size} " \
                  "--norm {norm} " \
                  "--out={out_file} " \
                  "--quiet " \
                  "{cov_file}".format(
                              bin=self.bin,
                              size='raw',
                              norm='log',
                              out_file=test_result_path,
                              cov_file=os.path.join(self._data_dir, '4_htmp.cov')
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
            msg = "craw_htmp execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "craw_htmp finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))
        for sense in ('sense', 'antisense'):
            filename, suffix = os.path.splitext(output_filename)
            filename = "{}.{}{}".format(filename, sense, suffix)
            expected_result_path = os.path.join(self._data_dir, filename)
            expected_md5 = hashlib.md5()
            expected_md5.update(open(expected_result_path,'rb').read())

            result_path, suffix = os.path.splitext(test_result_path)
            result_path = "{}.{}{}".format(result_path, sense, suffix)
            result_md5 = hashlib.md5()
            result_md5.update(open(result_path,'rb').read())
            self.assertEqual(expected_md5.hexdigest(), result_md5.hexdigest())


    def test_raw_log(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'htmp_raw_log.png'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--size {size} " \
                  "--norm {norm} " \
                  "--out={out_file} " \
                  "--quiet " \
                  "{cov_file}".format(
                              bin=self.bin,
                              size='raw',
                              norm='log',
                              out_file=test_result_path,
                              cov_file=os.path.join(self._data_dir, '4_htmp.cov')
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
            msg = "craw_htmp execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "craw_htmp finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))
        for sense in ('sense', 'antisense'):
            filename, suffix = os.path.splitext(output_filename)
            filename = "{}.{}{}".format(filename, sense, suffix)
            expected_result_path = os.path.join(self._data_dir, filename)
            expected_md5 = hashlib.md5()
            expected_md5.update(open(expected_result_path,'rb').read())
            expected_im = Image.open(expected_result_path)

            result_path, suffix = os.path.splitext(test_result_path)
            result_path = "{}.{}{}".format(result_path, sense, suffix)
            result_md5 = hashlib.md5()
            result_md5.update(open(result_path,'rb').read())
            result_im = Image.open(result_path)
            self.assertEqual(expected_md5.hexdigest(), result_md5.hexdigest())


    def test_raw_log_row(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw_test')
        os.makedirs(self.out_dir)
        output_filename = 'htmp_raw_log+row.png'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} " \
                  "--size {size} " \
                  "--norm {norm} " \
                  "--out={out_file} " \
                  "--quiet " \
                  "{cov_file}".format(
                              bin=self.bin,
                              size='raw',
                              norm='log+row',
                              out_file=test_result_path,
                              cov_file=os.path.join(self._data_dir, '4_htmp.cov')
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
            msg = "craw_htmp execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "craw_htmp finished with non zero exit code: {0} command launched=\n{1}\n{2}".format(
                            cov_process.returncode,
                            command,
                            ''.join([l.decode('utf-8') for l in cov_process.stderr.readlines()]),
                            ))
        for sense in ('sense', 'antisense'):
            filename, suffix = os.path.splitext(output_filename)
            filename = "{}.{}{}".format(filename, sense, suffix)
            expected_result_path = os.path.join(self._data_dir, filename)
            expected_md5 = hashlib.md5()
            expected_md5.update(open(expected_result_path,'rb').read())

            result_path, suffix = os.path.splitext(test_result_path)
            result_path = "{}.{}{}".format(result_path, sense, suffix)
            result_md5 = hashlib.md5()
            result_md5.update(open(result_path,'rb').read())
            self.assertEqual(expected_md5.hexdigest(), result_md5.hexdigest())


