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

import re
from abc import ABCMeta, abstractmethod


class WigException(Exception):
    pass


class ExpandedStrand:

    def __init__(self, start, stop, span):
        """

        :param start: start position of the chunk (1 based position)
        :param stop: stop position of the chunk (1 based position)
        :param span: length of span
        """
        self._start = start
        # be careful the last position count in the span
        self.coverage = [0] * (stop + span - start)

    def __getattr__(self, item):
        return getattr(self.coverage, item)

    def __len__(self):
        return len(self.coverage)

    @property
    def stop(self):
        return len(self.coverage) + self._start - 1

    @property
    def start(self):
        return self._start

    def __eq__(self, other):
        return self._start == other._start and self.coverage == other.coverage

    def __setitem__(self, idx, value):
        self.coverage[idx - self.start] = value


class Chunk(metaclass=ABCMeta):

    def __init__(self, **kwargs):
        self.span = 1
        for k, v in kwargs.items():
            if k in ('span', 'start', 'step'):
                v = int(v)
            setattr(self, k, v)
        try:
            if not self.chrom:
                raise WigException("'chrom' field  is not present.")
        except AttributeError:
            raise WigException("'chrom' field  is not present.")

        if self.span <= 0:
            raise WigException("'{}' is not allowed as span value.".format(self.span))

        self.forward = []
        self.reverse = []

    @abstractmethod
    def is_fixed_step(self):
        pass

    @abstractmethod
    def parse_data_line(self):
        pass

    def expand(self):
        expanded_strands = {}
        for sense in ('forward', 'reverse'):
            chunk_strand = getattr(self, sense)
            if chunk_strand:
                expanded_strands[sense] = ExpandedStrand(chunk_strand[0][0], chunk_strand[-1][0], self.span)
                for position, cov in chunk_strand:
                    for i in range(position, position + self.span):
                        expanded_strands[sense][i] = cov
            else:
                # there is no data for this strand
                expanded_strands[sense] = ExpandedStrand(0, 0, 0)
        return expanded_strands['forward'], expanded_strands['reverse']


class FixedChunk(Chunk):

    def __init__(self, **kwargs):
        self.step = 1
        super(FixedChunk, self).__init__(**kwargs)
        if self.step == 0:
            raise WigException
        try:
            self.start
        except AttributeError:
            raise WigException("'start' must be defined for 'fixedStep'.")
        if self.span > self.step:
            raise WigException("'span' cannot be greater than 'step'.")
        self._current_pos = self.start

    def is_fixed_step(self):
        return True

    def parse_data_line(self, line):
        cov = float(line)
        if cov >= 0:
            self.forward.append((self._current_pos, cov))
        else:
            self.reverse.append((self._current_pos, abs(cov)))
        self._current_pos += self.step


class VariableChunk(Chunk):

    def is_fixed_step(self):
        return False

    def parse_data_line(self, line):
        pos, cov = line.split()
        pos = int(pos)
        cov = float(cov)
        if cov >= 0:
            self.forward.append((pos, cov))
        else:
            self.reverse.append((pos, abs(cov)))


class Strand:

    def __init__(self):
        self.coverage = []

    def __len__(self):
        return len(self.coverage)

    def __eq__(self, other):
        return self.coverage == other.coverage

    def extend(self, strand):
        self.coverage.extend(strand.coverage)

    def get_coverage(self, start, stop):
        cov_len = len(self.coverage)
        if 0 <= start <= cov_len:
            if stop < cov_len:
                return self.coverage[start: stop + 1]
            else:
                cov = self.coverage[start:]
                right_fill = [0] * (stop - start - len(cov))
                return cov + right_fill
        elif start < 0:
            left_fill = [0] * abs(start)
            if stop < cov_len:
                return left_fill + self.coverage[: stop + 1]
            else:
                cov = self.coverage
                right_fill = [0] * (stop - cov_len)
                return left_fill + cov + right_fill
        elif start > cov_len:
            return [0] * (stop - start)


class Chromosome:

    def __init__(self, name):
        self.name = name
        self._forward = Strand()
        self._reverse = Strand()

    @property
    def forward(self):
        return self._forward

    @property
    def reverse(self):
        return self._reverse

    def __iadd__(self, chunk):
        if not isinstance(chunk, Chunk):
            raise TypeError("can add only Chunk objects, provide: {}".chunk.__class__.__name__)
        forward, reverse = chunk.expand()
        for sense in ('forward', 'reverse'):
            my_strand = getattr(self, sense)
            strand_2_add = locals()[sense]
            # we switch from real position to zero based position
            # and we stop the fill at the position before the fragt start so => -2
            fill = ExpandedStrand(len(my_strand), strand_2_add.start -2, 1)
            my_strand.extend(fill)
            my_strand.extend(strand_2_add)
        return self


class Genome:

    def __init__(self):
        self._chromosomes = {}
        self.infos = {}

    def __getitem__(self, name):
        return self._chromosomes[name]

    def __contains__(self, chrom):
        if isinstance(chrom, str):
            return chrom in self._chromosomes
        elif not isinstance(chrom, Chromosome):
            raise TypeError("'in <Genome>' requires string or Chromosome as left operand, not {}".format(type(chrom)))
        else:
            return chrom.name in self._chromosomes

    def add(self, chrom):
        if not isinstance(chrom, Chromosome):
            raise TypeError("Genome can contains only Chromosome objects")
        self._chromosomes[chrom.name] = chrom


class WigParser:

    def __init__(self, path):
        self.declaration_type_pattern = re.compile('fixedStep|variableStep')
        self._path = path
        self._genome = None
        self._current_chunk = None
        self._current_chrom = None


    def parse(self):
        with open(self.path, 'r') as wig_file:
            self._genome = Genome()
            for line in wig_file:
                line = line.strip()
                if not line or self.is_comment_line(line):
                    continue
                elif self.is_track_line(line):
                    self.parse_track_line(line)
                elif self.is_declaration_line():
                    self.parse_track_line(line)
                else:
                    self._current_chunk['parse_data_line'](line)
            # we are at the end of the file
            # so add the last chunk to the others
            self._current_chrom += self._current_chunk
        return self._genome


    def parse_track_line(self, line):
        fields = line.re.findall('(\w+=\S+|".*"|\'.*\')', line)
        attrs = {}
        for attr, val in fields.split('='):
            attrs[attr] = val
        if 'type' not in attrs:
            raise WigException('wiggle type is not present.')
        else:
            self._genome.infos = attrs


    def parse_declaration_line(self, line):
        if self._current_chunk:
            self._current_chrom += self._current_chunk
            self._current_chunk = None

        fields = line.split()
        type = fields[0]

        kwargs = {attr: val for attr, val in [f.split('=') for f in fields[1:]]}

        if type == 'fixedStep':
            self._current_chunk = FixedChunk(**kwargs)
        else:
            self._current_chunk = VariableChunk(**kwargs)

        chrom_name = self._current_chunk.chrom
        if chrom_name in self._genome:
            chrom = self._genome[chrom_name]
        else:
            chrom = Chromosome(chrom_name)
            self._genome.add(chrom)
        self._current_chrom = chrom


    def parse_data_line(self, line):
        if self._current_chunk is None:
            raise WigException("this data line '{}' is not preceded by declaration".format(line))
        return self._current_chunk.parse_data_line(line)


    @staticmethod
    def is_track_line(line):
        return line.startswith('track')


    def is_declaration_line(self, line):
        return bool(re.match(self.declaration_type_pattern, line))


    def is_comment_line(self, line):
        return line.startswith('#')



