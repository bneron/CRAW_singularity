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
import logging


root_logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(fmt='{levelname} : {name} : {message}', style='{')
handler.setFormatter(formatter)
root_logger.addHandler(handler)
root_logger.setLevel(logging.DEBUG)

_log = logging.getLogger(__name__)
#level = logging.WARNING - (10 * args.verbosity)
#log.setLevel(logging.NOTSET)


class WigException(Exception):
    pass


class Coverage:

    def __init__(self, start, stop, span):
        """

        :param start: start position of the chunk (1 based position)
        :param stop: stop position of the chunk (1 based position)
        :param span: length of span
        """
        self._start = start
        # be careful the last position count in the span
        self._coverage = [0] * (stop + span - start)
        self._stop = len(self.coverage) + self._start - 1

    def __getattr__(self, item):
        return getattr(self.coverage, item)

    def __len__(self):
        return len(self.coverage)

    @property
    def stop(self):
        return self._stop

    @property
    def start(self):
        return self._start

    @property
    def coverage(self):
        return self._coverage[:]

    def __eq__(self, other):
        return self._start == other._start and self.coverage == other.coverage

    def _translate_pos(self, pos):
        if isinstance(pos, int):
            idx = pos - self._start
        elif isinstance(pos, slice):
            idx = slice(pos.start - self._start, pos.stop - self._start, None)
        return pos

    def __setitem__(self, pos, value):
        self.coverage[self._translate_pos(pos)] = value

    def __getitem__(self, pos):
        return self.coverage[self._translate_pos(pos)]


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
        coverages = {}
        for sense in ('forward', 'reverse'):
            chunk_strand = getattr(self, sense)
            if chunk_strand:
                coverages[sense] = Coverage(chunk_strand[0][0], chunk_strand[-1][0], self.span)
                for position, cov in chunk_strand:
                    for i in range(position, position + self.span):
                        coverages[sense][i] = cov
            else:
                # there is no data for this strand
                coverages[sense] = Coverage(0, 0, 0)
        return coverages['forward'], coverages['reverse']


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


    def get_coverage(self, start, stop):
        cov_len = len(self.coverage)
        right_fill = []
        left_fille = []
        first_chunk_start = self.coverage[0][0]
        last_chunk_stop = self.coverage[-1][1]

        #start < first_chunk_start
        #   right_fill = [0] * (first_chunk_start - start)
        #   start = first_chunk_start
        #start > last_chunk_stop

        #stop < fisrt_chunk_start
        #stop > last_chunk_stop
        #   left_fill = [0] * (stop - last_chunk_stop)
        #   stop = last_chunk_stop

        # recherche de l'inteval contenant le start
        # rechecher de l'interval contenant le stop
        # generer la list des interval
        # expandre la liste
        # ajouter r_fill + left_fill
        # retourner le resultat
        first_chunk = None
        previous_idx = None
        coverage_walker = enumerate(self.coverage)
        for idx, chunk in coverage_walker:
            if start < chunk.start:
                continue
            elif chunk.start <= start <= chunk.stop:
                first_chunk = idx
            elif start < chunk.start:
                first_chunk = previous_idx
            previous_idx = idx

        # le stop est forcement apres le start
        # pas la peine de reparcourir la liste depuis le debut
        # par contre il peut etre dans le meme chunk que le start
        # il faut donc le tester

        for dx, chunk in coverage_walker:
            if stop > chunk.stop:
                continue
            elif chunk.start <= stop <= chunk.stop:
                last_chunk = idx
            elif stop > chunk.stop:
                last_chunk = previous_idx
            previous_idx = idx

class Chromosome:
    """
    Handle chromosomes. a chromosome as a name and contains chunks (forward and reverse)
    """

    def __init__(self, name):
        self.name = name
        self._forward = []
        self._reverse = []

    @property
    def forward(self):
        return self._forward

    @property
    def reverse(self):
        return self._reverse

    def add_chunk(self, chunk):
        if not isinstance(chunk, Chunk):
            raise TypeError("can add only Chunk objects, provide: {}".chunk.__class__.__name__)
        forward, reverse = chunk.expand()
        for sense in ('forward', 'reverse'):
            my_strand = getattr(self, sense)
            one_strand_cov = locals()[sense]
            my_strand.append(one_strand_cov)
        return self


    @staticmethod
    def _find_start_chunk(start, strand):
        first_chunk_start = strand[0][0]
        last_chunk_stop = strand[-1][1]
        if start < first_chunk_start:
            return 0
        if start > last_chunk_stop:
            return None
        else:
            previous_idx = None
            for idx, cov in enumerate(strand):
                if start < cov.start and previous_idx is None:
                    continue
                elif cov.start <= start <= cov.stop:
                    return idx
                elif start < cov.start:
                    return previous_idx
                else:
                    # start > cov.start => test next chunk
                    previous_idx = idx
            # we are sure to find a chunk
            # because we first test that start <= last_chunk_stop


    @staticmethod
    def _find_stop_chunk(stop, start_idx, strand):
        first_chunk_start = strand[0][0]
        last_chunk_stop = strand[-1][1]
        if stop > last_chunk_stop:
            return 0
        if stop < first_chunk_start:
            return None
        else:
            previous_idx = None
            for idx in range(start_idx, len(strand) - start_idx):
                # stop is necessarily greater than start
                # so don't search stop search from the beginning
                cov = strand[idx]
                if stop > cov.stop and previous_idx is None:
                    continue
                elif cov.start <= stop <= cov.stop:
                    return idx
                elif stop < cov.stop:
                    return previous_idx
                else:
                    # stop < cov.stop
                    previous_idx = idx


    @staticmethod
    def _expand_cov(start, stop, strand):
        chunks = strand[start, stop + 1]
        cov = [0] * (stop - start)
        for chunk in chunks:
            cov[chunk.start - start: chunk.stop - start] = chunk.covergae
        return cov

    def get_coverage(self, start, stop, sense=None):
        covs = []
        for sense in (self._forward, self._reverse):
            first_chunk_idx = self._find_start_chunk(start, sense)
            last_chunk_idx = self._find_stop_chunk(stop, start, sense)
            covs.append(self._expand_cov(first_chunk_idx, last_chunk_idx, sense))
        return covs


class Genome:
    """
    A genome is made of chromosomes and some metadata, called infos
    """

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
    """
    class to parse file in wig format.
    at the end of parsing it returns a :class:`Genome` object.
    """

    def __init__(self, path):
        """
        :param path: the path of the wig file to parse.
        :type path: string
        """
        self.declaration_type_pattern = re.compile('fixedStep|variableStep')
        self._path = path
        self._genome = None
        self._current_chunk = None
        self._current_chrom = None


    def parse(self):
        """
        Open a wig file and parse it.
        read wig file line by line check the type of line
        and call the corresponding method accordingly the type of the line:
        - comment
        - track
        - declaration
        - data
        see
        - https://wiki.nci.nih.gov/display/tcga/wiggle+format+specification
        - http://genome.ucsc.edu/goldenPath/help/wiggle.html
        for wig specifications.
        This parser does not fully follow these specification. When a score is negative,
        it means that the coverage is on the reverse strand. So some positions can apear twice
        in one block of declaration (what I call a chunk).

        :return: a Genome coverage corresponding to te wig file
        :rtype: :class:`Genome` object
        """
        with open(self._path, 'r') as wig_file:
            self._genome = Genome()
            for line in wig_file:
                line = line.strip()
                if not line or self.is_comment_line(line):
                    continue
                elif self.is_track_line(line):
                    self.parse_track_line(line)
                elif self.is_declaration_line(line):
                    self.parse_declaration_line(line)
                else:
                    self._current_chunk.parse_data_line(line)
            # we are at the end of the file
            # so add the last chunk to the others
            self._current_chrom += self._current_chunk
        return self._genome


    def parse_track_line(self, line):
        """
        fill the genome infos with the information found on the track.

        :param line: line to parse. The method :meth:`is_track_line` must return True with this line.
        """
        _log.info('parsing : {}'.format(line))
        fields = re.findall("""(\w+)=(".+?"|'.+?'|\S+)""", line)
        attrs = {}
        for attr, val in fields:
            attrs[attr] = val
        if 'type' not in attrs:
            raise WigException('wiggle type is not present.')
        else:
            self._genome.infos = attrs


    def parse_declaration_line(self, line):
        """
        Create a new chunk, and set the current_chunk and current_chromosome (create a new one if necessary).

        :param line: line to parse. The method :meth:`is_declaration_line` must return True with this line.
        """
        _log.info("parsing : {}".format(line))
        if self._current_chunk:
            self._current_chrom.add_chunk(self._current_chunk)
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
        """
        :param line: line to parse. It must not a comment_line, neither a track line nor a declaration line.
        :type line: string
        :return:
        :rtype: 
        """
        if self._current_chunk is None:
            raise WigException("this data line '{}' is not preceded by declaration".format(line))
        return self._current_chunk.parse_data_line(line)


    @staticmethod
    def is_track_line(line):
        """
        :param line: line to parse.
        :type line: string
        :return: True if line is a track line. False otherwise.
        :rtype: boolean
        """
        return line.startswith('track')


    def is_declaration_line(self, line):
        """
        :param line: line to parse.
        :type line: string
        :return: True if line is a decalration line. False otherwise.
        :rtype: boolean
        """
        return bool(re.match(self.declaration_type_pattern, line))


    def is_comment_line(self, line):
        """
        :param line: line to parse.
        :type line: string
        :return: True if line is a comment line. False otherwise.
        :rtype: boolean
        """
        return line.startswith('#')



