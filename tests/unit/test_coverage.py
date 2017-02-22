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


import os
import pysam

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.coverage import get_coverage
from craw.annotation import new_entry_type


class TestCoverage(CRAWTest):

    def test_get_coverage_fix_window(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        sam_file = pysam.AlignmentFile(sam_path, "rb")
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col)
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 14415],
                       ['YEL071W', 'DLD3', 'chrV', '+', 17848],
                       ['YEL071W', 'DLD3', 'chrV', '+', 4]
                       ]

        expected = [{'for': [0, 0, 0, 0, 0, 0, 0, 0],
                    'rev': [0, 0, 0, 0, 0, 0, 0, 0]
                     },
                    {'for': [227, 227, 227, 227, 226, 225, 224, 224],
                     'rev': [0, 0, 0, 0, 0, 0, 0, 0]
                     },
                    {'for': [None, 0, 0, 0, 0, 0, 0, 0],
                     'rev': [None, 0, 0, 0, 0, 0, 0, 0]
                     }
                    ]
        for values, exp_val in zip(value_lines, expected):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_coverage(sam_file,
                                                    annot_entry,
                                                    start=values[-1] - 5,
                                                    stop=values[-1] + 3,
                                                    qual_thr=0,
                                                    max_left=0,
                                                    max_right=0)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])


    def test_get_coverage_var_window(self):
        sam_path = os.path.join(self._data_dir, 'small.bam')
        sam_file = pysam.AlignmentFile(sam_path, "rb")
        annot_fields = ['name', 'gene', 'chromosome', 'strand', 'Position','beg', 'end']
        entry_cls_name = 'foo'
        ref_col = 'Position'
        ne_class = new_entry_type(entry_cls_name, annot_fields, ref_col, start_col='beg', stop_col='end')
        value_lines = [['YEL072W', 'RMD6', 'chrV', '+', 14415, 14412, 14419],
                       ['YEL071W', 'DLD3', 'chrV', '+', 17848, 17840, 17850]]

        expected = [{'for': [None, None, None, None, None, 0, 0, 0, 0, 0, 0, 0, None],
                     'rev': [None, None, None, None, None, 0, 0, 0, 0, 0, 0, 0, None]
                    },
                    {'for': [227, 227, 227, 227, 227, 227, 227, 226, 225, 224, None, None, None],
                     'rev': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, None, None, None]
                    }
                   ]

        for values, exp_val in zip(value_lines, expected):
            annot_entry = ne_class([str(v) for v in values])
            forward_cov, reverse_cov = get_coverage(sam_file,
                                                    annot_entry,
                                                    start=annot_entry.start,
                                                    stop=annot_entry.stop,
                                                    qual_thr=0,
                                                    max_left=8,
                                                    max_right=5)
            self.assertListEqual(forward_cov, exp_val['for'])
            self.assertListEqual(reverse_cov, exp_val['rev'])
