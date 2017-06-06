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


from collections import namedtuple


class Entry:
    """Handle one entry (One line) of annotation file """

    def __init__(self, values):
        """

        :param values: the values parsed from one line of the annotation file
        :type values: list of string
        """
        self._reverse_values = set(('-', '-1', 'rev'))
        self._forward_values = set(('+', '1', 'for'))

        if len(values) != len(self._fields):
            raise RuntimeError("the number of values ({}) does not match with number of fields ({}): {}".format(len(values),
                                                                                                            len(self._fields),
                                                                                                                values))
        self._values = [self._convert(f, v) for f, v in zip(self._fields, values)]
        if self.start is not None:
            if not (self.start <= self.ref <= self.stop or self.start >= self.ref >= self.stop):
                raise RuntimeError("error in line '{line}': {ref_col} {ref} is not "
                                   "between {start_col}: {start} and {stop_col}: {stop}".format(
                                    line=self,
                                    start=self.start,
                                    start_col=self._fields_idx['start'].col_name,
                                    ref=self.ref,
                                    ref_col=self._fields_idx['ref'].col_name,
                                    stop=self.stop,
                                    stop_col=self._fields_idx['stop'].col_name)
                                    )
            if self.start > self.stop and self.strand == '-':
                self._switch_start_stop()
            elif self.start > self.stop:
                raise RuntimeError("error in line '{line}': {start_col}:{start} > {stop_col}: {stop}"
                                   "on forward strand".format(
                                                              line=self,
                                                              start=self.start,
                                                              start_col=self._fields_idx['start'].col_name,
                                                              stop=self.stop,
                                                              stop_col=self._fields_idx['stop'].col_name
                                                             )
                                   )


    def _convert(self, field, value):
        """
        Convert field parsed from annotation file in Entry internal value

        :param field: the field name associated to the value.
        :type field: string
        :param value: the value to convert
        :type value: string
        :return: the converted value
        :rtype: any
        :raise: RuntimeError or value Error  if a value cannot be converted
        """
        if field == self._fields_idx['ref'].col_name:
            value = int(value)
        elif field == self._fields_idx['strand'].col_name:
            v = value.lower()
            if v in self._forward_values:
                value = '+'
            elif v in self._reverse_values:
                value = '-'
            else:
                raise RuntimeError("strand must be '+/-', '1/-1' or 'for/rev' got '{}'".format(value))
        elif 'start' in self._fields_idx and field == self._fields_idx['start'].col_name:
            value = int(value)
        elif 'stop' in self._fields_idx and field == self._fields_idx['stop'].col_name:
            value = int(value)
        return value


    def _switch_start_stop(self):
        """
        Switch start and stop value if self.start > self.stop
        This situation can occur if annotation regards the reverse strand
        """
        start, stop = self.stop, self.start
        self._values[self._fields_idx['start'].idx] = start
        self._values[self._fields_idx['stop'].idx] = stop

    @property
    def chromosome(self):
        """The name of the Chromosome"""
        return self._values[self._fields_idx['chr'].idx]

    @property
    def ref(self):
        """The position of reference"""
        return self._values[self._fields_idx['ref'].idx]

    @property
    def strand(self):
        """the strand +/-"""
        return self._values[self._fields_idx['strand'].idx]

    @property
    def start(self):
        """The Position to start the coverage computation"""
        if 'start' in self._fields_idx:
            return self._values[self._fields_idx['start'].idx]
        else:
            return None

    @property
    def stop(self):
        """The position to end the coverage computation (included)"""
        if 'stop' in self._fields_idx:
            return self._values[self._fields_idx['stop'].idx]
        else:
            return None

    @property
    def header(self):
        """The header of the annotation file"""
        return '\t'.join(self._fields)

    def __str__(self):
        return '\t'.join([str(v) for v in self._values])

    def __eq__(self, other):
        for v, vo in zip(self._values, other._values):
            if v != vo:
                return False
        return True

Idx = namedtuple('Idx', ('col_name', 'idx'))


def new_entry_type(name, fields, ref_col,
                   strand_col='strand',
                   chr_col='chromosome',
                   start_col=None,
                   stop_col=None):
    """
    From the header of the annotation line create a new Entry Class inherited from Entry Class

    :param name: The name of the new class of entry.
    :type name: str
    :param fields: The fields constituting the new type of entry.
    :type fields: list of string
    :param ref_col: The name of the column representing the position of reference (default is 'position').
    :type ref_col: string
    :param strand_col: The name of the column representing the strand (default is 'strand').
    :type strand_col: string
    :param chr_col: The name of the column representing the name of chromosome (default is 'chromosome').
    :type chr_col: string
    :param start_col: The name of the column representing the position of the first base to compute the coverage (inclusive).
    :type start_col: string
    :param stop_col: The name of the column representing the position of the last base to compute the coverage (inclusive).
    :type stop_col: string
    :return: a new class child of :class:`Entry` which is able to store information corresponding to the header.
    """
    fields_idx = {}
    if any((start_col, stop_col)) and not all((start_col, stop_col)):
        raise RuntimeError("if start_col is specified stop_col must be specified too and vice versa")
    try:
        fields_idx['ref'] = Idx(ref_col, fields.index(ref_col))
    except ValueError:
        raise RuntimeError("The ref_col '{}' does not match any fields: '{}'\n"
                           "You must specify the '--ref-col' option".format(ref_col, ', '.join(fields))) from None
    try:
        fields_idx['strand'] = Idx(strand_col, fields.index(strand_col))
    except ValueError:
        raise RuntimeError("The strand_col '{}' does not match any fields: '{}'".format(strand_col, ', '.join(fields))) from None
    try:
        fields_idx['chr'] = Idx(chr_col, fields.index(chr_col))
    except ValueError:
        raise RuntimeError("The chr_col '{}' does not match any fields: '{}'\n"
                           "You must specify the '--chr-col' option".format(chr_col, ', '.join(fields))) from None

    if start_col:
        try:
            fields_idx['start'] = Idx(start_col, fields.index(start_col))
        except ValueError:
            raise RuntimeError("The start_col '{}' does not match any fields: '{}'".format(start_col, ', '.join(fields))) from None
        try:
            fields_idx['stop'] = Idx(stop_col, fields.index(stop_col))
        except ValueError:
            raise RuntimeError("The stop_col '{}' does not match any fields: '{}'".format(stop_col, ', '.join(fields))) from None
    return type(name, (Entry,), {'_fields_idx': fields_idx, '_fields': fields})


class AnnotationParser:
    """
    Parse the annotation file
         - create new type of Entry according to the header
         - create one Entry object for each line of the file
    """

    def __init__(self, path, ref_col, strand_col='strand', chr_col='chromosome',
                 start_col=None, stop_col=None, sep='\t'):
        """

        :param path: the path to the annotation file to parse.
        :type path: string
        :param ref_col: the name of the column for the reference position
        :type ref_col: string
        :param chr_col: the name of the column for the chromosome
        :type chr_col: string
        :param strand_col: the name of the column for the strand
        :type strand_col: string
        :param start_col: the name of the column for start position
        :type start_col: string
        :param stop_col: the name of the column for the stop position
        :type stop_col: string
        :param sep: The separator tu use to split fields
        :type sep: string
        """
        self.path = path
        self.ref_col = ref_col
        self.chr_col = chr_col
        self.strand_col = strand_col
        self.start_col = start_col
        self.stop_col = stop_col
        self._sep = sep
        with open(self.path, 'r') as annot_file:
            self.header = annot_file.readline().rstrip('\n').split(self._sep)

    def get_annotations(self):
        """
        Parse an annotation file and yield a :class:`Entry` for each line of the file.

        :return: a generator on a annotation file.
        """
        with open(self.path, 'r') as annot_file:
            _ = annot_file.readline()
            MyEntryClass = new_entry_type('MyEntry', self.header, self.ref_col,
                                          strand_col=self.strand_col,
                                          chr_col=self.chr_col,
                                          start_col=self.start_col,
                                          stop_col=self.stop_col)
            for line in annot_file:
                yield MyEntryClass(line.rstrip('\n').split(self._sep))

    def max(self):
        """
        :return: the maximum of bases to take in count before and after the reference position.
        :rtype: tuple of 2 int
        """
        if self.start_col is not None:
            max_left = max_right = 0
            for entry in self.get_annotations():
                left = entry.ref - entry.start
                right = entry.stop - entry.ref
                if entry.strand == '-':
                    # for coverages on reverse strand the sense of scores
                    # must be reversed, so the computation of padding must be too
                    left, right = right, left
                if left > max_left:
                    max_left = left
                if right > max_right:
                    max_right = right
            return max_left, max_right
        else:
            return 0, 0
