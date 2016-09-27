# -*- coding: utf-8 -*-

import shutil
import tempfile
import os
from subprocess import Popen, PIPE ,DEVNULL


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
        self.bin = os.path.join(self.craw_home, 'bin', 'coverage') if self.local_install else which('coverage')

    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass


    def test_with_fixed_window(self):
        """
        | test if returncode of coverage is 0 and
        | then test if the generated file is the same as a reference file
        """
        self.out_dir = os.path.join(self.tmp_dir, 'craw')
        os.makedirs(self.out_dir)
        output_filename = 'coverage_with_fixed_window.cov.txt'
        test_result_path = os.path.join(self.out_dir, output_filename)
        command = "{bin} --bam={bam_file} --annot={annot_file} " \
                  "--before={before} --after={after} " \
                  "--ref-col={ref_col} " \
                  "--output={out_file} ".format(
                                                 bin=self.bin,
                                                 bam_file=os.path.join(self._data_dir, 'crac_tac4_20160624_plus_uv_sorted_unspliced.bam'),
                                                 annot_file=os.path.join(self._data_dir, 'annotations_cerevisiae_10.txt'),
                                                 ref_col='Position',
                                                 before=5,
                                                 after=3,
                                                 out_file=test_result_path
                                                )
        #print "\n", command
        if not self.bin:
            raise RuntimeError('coverage not found, CRAW_HOME must be either in your path or CRAW_HOME must be defined '
                               'command launch: \n{}'.format(command))

        try:
            cov_process = Popen(command,
                                shell=True,
                                stdin=None,
                                stderr=DEVNULL,
                                close_fds=False
                                )
        except Exception as err:
            msg = "coverage execution failed: command = {0} : {1}".format(command, err)
            print()
            print(msg)
            raise err from None

        cov_process.wait()
        self.assertEqual(cov_process.returncode, 0,
                         "coverage finished with non zero exit code: {0} command launched=\n{1}".format(
                         cov_process.returncode,
                         command))

        expected_result_path = os.path.join(self._data_dir, output_filename)
        with open(expected_result_path) as expected_result_file:
            expected_result = expected_result_file.readlines()


        with open(test_result_path) as test_result_file:
            test_result = test_result_file.readlines()

        self.assertListEqual(expected_result, test_result)


