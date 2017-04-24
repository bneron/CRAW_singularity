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

import logging

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

import craw

class MyTestCase(CRAWTest):


    def test_get_version_message(self):
        self.assertTrue(craw.get_version_message().startswith(
            "craw NOT packaged, it should be a development version | Python 3."))
        craw.__version__ = 1.0
        self.assertTrue(craw.get_version_message().startswith("craw 1.0 | Python 3."))


    def test_init_logger(self):
        log_level = 30
        craw.init_logger(log_level)

        craw_log = logging.getLogger('craw')
        self.assertEqual(len(craw_log.handlers), 1)
        self.assertEqual(craw_log.handlers[0].__class__.__name__, 'StreamHandler')
        self.assertEqual(craw_log.getEffectiveLevel(), log_level)