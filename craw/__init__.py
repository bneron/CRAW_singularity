import sys

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

__version__ = '$VERSION'

def get_version_message():
    """

    :return: A human readable version of the craw version
    :rtype: string
    """
    # if I keep '$ VERSION' (without space) as is
    # the setup.py will replace it by the value set in setup
    # so the test become True even if craw is installed using setup.py
    if __version__ == '$' + 'VERSION':
        version = "NOT packaged, it should be a development version"
    else:
        version = __version__
    version_text = "craw {0} | Python {1}.{2}".format(version, sys.version_info.major, sys.version_info.minor)
    return version_text