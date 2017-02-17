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


import argparse


class VersionAction(argparse._VersionAction):
    """Class to allow argparse to handel more complex version output"""

    def __call__(self, parser, namespace, values, option_string=None):
        """Override the :meth:`argparse._VersionAction.__call__` to use
           a RawTextHelpFormatter only for version action whatever the class_formatter
           specified for the :class:`argparse.ArgumentParser` object.
        """
        version = self.version
        if version is None:
            version = parser.version
        formatter = argparse.RawTextHelpFormatter(parser.prog)
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), argparse._sys.stdout)
        parser.exit()


