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
import collections
from abc import ABCMeta, abstractmethod
import logging

import numpy as np


root_logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(fmt='{levelname} : {name} : {message}', style='{')
handler.setFormatter(formatter)
root_logger.addHandler(handler)
root_logger.setLevel(logging.DEBUG)

_log = logging.getLogger(__name__)
#level = logging.WARNING - (10 * args.verbosity)
#log.setLevel(logging.NOTSET)


class WigError(Exception):
    pass


class Coverage:
    """
    Represent the coverage for a piece of genome (only one strand).
    There is one coverage value for each nucleotide of this genome part.
    """

    def __init__(self, start, stop, span, default_cov=0.0):
        """

        :param start: start position of the chunk (1 based position)
        :param stop: stop position of the chunk (1 based position)
        :param span: length of span
        """

        print("Coverage creating coverage of {} len".format(stop + span - start))
        self._start = start
        # be careful the last position count in the span
        self._coverage = np.array([default_cov] * (stop + span - start))
        self._stop = stop + span - 1

    def __len__(self):
        return len(self._coverage)

    def __eq__(self, other):
        return self._start == other.start and \
               self._stop == other.stop and \
               self.coverage == other.coverage

    def __str__(self):
        s = "start= {}; stop= {}\n{}".format(self.start, self.stop, self._coverage)
        return s

    def _translate_pos(self, pos):
        """
        translate the position on chromosome on index in the list of float
        :param pos: the postion to set value
        :type pos: int or :class:`slice` object
        :return: the index or slice in the list of float
        :rtype: int or slice
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        if isinstance(pos, int):
            if not (self._start <= pos <= self._stop):
                raise IndexError('Coverage {{}} position {} out of range[{}:{}]'.format(pos, self._start, self._stop))
            idx = pos - self._start
        elif isinstance(pos, slice):
            if not (self._start <= pos.start <= self._stop and self._start <= pos.stop <= self._stop + 1):
                raise IndexError('Coverage {{}} position [{}:{}] out of range [{}:{}]'.format(pos.start, pos.stop,
                                                                                              self._start, self._stop))
            idx = slice(pos.start - self._start, pos.stop - self._start, pos.step)
        return idx


    def __setitem__(self, pos, value):
        """

        :param pos: the postion to set value
        :type pos: int or :class:`slice` object
        :param value: value to assign
        :type value: float or iterable of float
        :raise ValueError: when pos is a slice and value have not the same length of the slice
        :raise TypeError: when pos is a slice and value is not iterable
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        if isinstance(pos, slice):
            if isinstance(value, collections.Iterable):
                if (pos.stop - pos.start) != len(value):
                    raise ValueError("can assign only iterable of same length of the slice")
            else:
                raise TypeError('can only assign an iterable')
        try:
            self._coverage[self._translate_pos(pos)] = value
        except IndexError as err:
            raise IndexError(str(err).format('assignment')) from None


    def __getitem__(self, pos):
        """

        :param pos: a position or a slice
        :return: the coverage at this position or corresponding to this slice.
        :rtype: a float o list of float
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        try:
            return self._coverage[self._translate_pos(pos)]
        except IndexError as err:
            raise IndexError(str(err). format('')) from None


    @property
    def stop(self):
        """
        :return: the stop position on the chromosome.
        :rtype: float
        """
        return self._stop

    @property
    def start(self):
        """
        :return: the start position on the chromosome.
        :rtype: float
        """
        return self._start

    @property
    def coverage(self):
        """
        :return: a copy of the coverage values in a list
        :rtype: list of float
        """
        return self._coverage.copy()

    @staticmethod
    def join(coverages, glue=0.0):
        """
        
        :param coverages:
        :type coverages: list or tuple of :class:`Coverage`
        :param glue:
        :type glue: float
        :return: the coverage corresponding to all coverage in coverages, 
                 the gaps between coverages are filled with glue.
        :rtype: :class:`Coverage` object.
        :raise ValueError: if coverages is empty.
        """
        print("Coverage join", coverages)
        if not coverages:
            raise ValueError("No coverage object to join")
        elif len(coverages) == 1:
            return coverages[0]
        start = min([c.start for c in coverages])
        stop = max([c.stop for c in coverages])
        joined_cov = Coverage(start, stop, 1, default_cov=glue)
        for cov in coverages:
            joined_cov[cov.start:cov.stop + 1] = cov.coverage
        return joined_cov


class Chunk(metaclass=ABCMeta):
    """
    Represent the data following a declaration line.
    The a Chunk contains sparse data on coverage 
    on a region of one chromosomes on both strand plus data contains on the declaration line.
    """

    def __init__(self, **kwargs):
        """
        :param kwargs: the key,values pairs found on a Declaration line
        :type kwargs: dictionary
        """
        self.span = 1
        self.start = None
        for k, v in kwargs.items():
            if k in ('span', 'start', 'step'):
                v = int(v)
            setattr(self, k, v)
        try:
            if not self.chrom:
                raise WigError("'chrom' field  is not present.")
        except AttributeError:
            raise WigError("'chrom' field  is not present.")

        if self.span <= 0:
            raise WigError("'{}' is not allowed as span value.".format(self.span))

        self.forward = []
        self.reverse = []

    @property
    @abstractmethod
    def stop(self):
        return NotImplemented

    @abstractmethod
    def is_fixed_step(self):
        """
        This is an abstract methods, must be implemented in inherited class
        :return: True if i's a fixed chunk of data, False otheweise
        :rtype: boolean
        """
        return NotImplemented

    @abstractmethod
    def parse_data_line(self):
        """
        parse a line of data and append the results in the corresponding strand
        This is an abstract methods, must be implemented in inherited class.
        """
        return NotImplemented

    def to_coverages(self):
        """
        :return: The coverage for the forward and reverse (in this order) for this chunk.
        :rtype: tuple (:class:`Coverage`, :class:`Coverage`)
        """
        print("Chunk.to_coverages")
        coverages = {}
        for sense in ('forward', 'reverse'):
            # the chunk stop take already the span in account
            coverages[sense] = Coverage(self.start, self.stop - self.span + 1, self.span)
            chunk_strand = getattr(self, sense)
            if chunk_strand:
                for position, cov in chunk_strand:
                    coverages[sense][position:position + self.span] = [cov] * self.span
        return coverages['forward'], coverages['reverse']


class FixedChunk(Chunk):
    """
    The FixedChunk objects handle data of 'fixedStep' declaration line and it's coverage data 
    """

    def __init__(self, **kwargs):
        self.step = 1
        super().__init__(**kwargs)
        if self.step <= 0:
            raise WigError("'step' must be strictly positive.")
        if self.start is None:
            raise WigError("'start' must be defined for 'fixedStep'.")
        if self.span > self.step:
            raise WigError("'span' cannot be greater than 'step'.")
        self._current_pos = self.start

    @property
    def stop(self):
        return ((len(self.forward) + len(self.reverse) - 1) * self.step) + (self.span - 1) + self.start


    def is_fixed_step(self):
        """
        :return: True
        :rtype: boolean
        """
        return True

    def parse_data_line(self, line):
        """
        parse line of data following a fixedStep Declaration.
        add the result on the corresponding strand (forward if coverage value is positive, reverse otherwise)
        :param line: line of data to parse (the white spaces at the end must be strip)
        :type line: string
        """
        cov = float(line)
        if cov >= 0:
            self.forward.append((self._current_pos, cov))
        else:
            self.reverse.append((self._current_pos, abs(cov)))
        self._current_pos += self.step


class VariableChunk(Chunk):
    """
    The Variable Chunk objects handle data of 'variableStep' declaration line and it's coverage data 
    
    If in data there is negative values this indicate that the coverage match on the reverse strand.
    the chunk start with the smallest position and end to the higest position whatever on wich strand are
    these position. This mean that when the chunk will be convert in Coverage,
    the lacking positions will be filled with 0.0.
    
    for instance: 
      
      variableStep chrom=chr3 span=2
      10 11
      20 22
      20 -30
      25 -50
      
    will give coverages starting at position 10 and ending at 26 for both strands and with 
    the following coverages values
    for = [11.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    rev = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 30.0, 0.0, 0.0, 0.0, 50.0, 50.0]
    """

    @property
    def stop(self):
        return self.forward[-1][0] + self.span - 1

    def is_fixed_step(self):
        """
        :return: False
        :rtype: boolean
        """
        return False


    def parse_data_line(self, line):
        """
        Parse line of data following a variableStep Declaration.
        Add the result on the corresponding strand (forward if coverage value is positive, reverse otherwise)
        :param line: line of data to parse (the white spaces at the end must be strip)
        :type line: string
        """
        pos, cov = line.split()
        pos = int(pos)
        cov = float(cov)
        if self.start is None:
            self.start = pos
        if cov >= 0:
            self.forward.append((pos, cov))
        else:
            self.reverse.append((pos, abs(cov)))


class Chromosome:
    """
    Handle chromosomes. A chromosome as a name and contains :class:`Chunk` objects 
    (forward and reverse)
    """

    def __init__(self, name):
        self.name = name
        self._coverage = np.full((2, 1000000), np.nan)


    def get_coverage(self, start_0_based, stop_0_based):
        """
        To mimic the api of :meth:`pysam.calignmentfile.AlignmentFile.count_coverage`
        The position are 0-based and the stop position is excluded
        Which is different from other wig object which are 1-based position.
        
        :param start: the start of the region (0 based position, included)
        :type start: int
        :param stop: the end of the region (0 based position, excluded)
        :type stop: int
        :return: the coverage corresponding to this region on the both strands.
        :rtype: tuple of 2 lists of floats)
        """
        return


    def __setitem__(self, pos, value):
        """

        :param pos: the postion (0-based) to set value
        :type pos: int or :class:`slice` object
        :param value: value to assign
        :type value: float or iterable of float
        :raise ValueError: when pos is a slice and value have not the same length of the slice
        :raise TypeError: when pos is a slice and value is not iterable
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        strand = 0 if value <= 0 else 1

        if isinstance(pos, slice):
            if isinstance(value, collections.Iterable):
                if (pos.stop - pos.start) != len(value):
                    raise ValueError("can assign only iterable of same length of the slice")
            else:
                raise TypeError('can only assign an iterable')

        while pos > self.coverage.shape[1]:
            self._extend()
        self._coverage[strand, pos] = value


    def __getitem__(self, pos):
        """

        :param pos: a position (0-based) or a slice
        if pos is a slice the left indice is excluded
        :return: the coverage at this position or corresponding to this slice.
        :rtype: a float o list of float
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        return self._coverage[0:pos], self._coverage[1:pos]


    def _extend(self, size=1000000, fill=np.nan):
        chunk = np.full((2, size), fill_value=fill)
        np.hstack(self._coverage, chunk)




class Genome:
    """
    A genome is made of chromosomes and some metadata, called infos
    """

    def __init__(self):
        self._chromosomes = {}
        self.infos = {}

    def __getitem__(self, name):
        """
        :param name: the name of the chromosome to retrieve
        :type name: string
        :return: the chromosome corresponding to the name.
        :rtype: :class:`Chromosome` object.
        """
        return self._chromosomes[name]

    def __contains__(self, chrom):
        if isinstance(chrom, str):
            return chrom in self._chromosomes
        elif not isinstance(chrom, Chromosome):
            raise TypeError("'in <Genome>' requires string or Chromosome as left operand, not '{}'".format(
                chrom.__class__.__name__))
        else:
            return chrom.name in self._chromosomes


    def __delitem__(self, name):
        """
        remove a chromosome from this genome
        
        :param name: the name of the chromosome to remove 
        :type name: string
        :return: None
        
        """
        if name in self._chromosomes:
            del self._chromosomes[name]
        else:
            raise KeyError("The chromosome '{}' is not in this genome.".format(name))


    @property
    def chromosomes(self):
        return list(self._chromosomes.values())


    def add(self, chrom):
        """
        add a chromosome in to a genome.
        if a chromosome with the same name already exist the previous one is replaced silently by this one.

        :param chrom: a chromosome to ad to this genome
        :type chrom: :class:`Chromosome` object.
        :raise: TypeError if chrom is not a :class:`Chromosome` object.
        """
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
            self._current_chrom.add_chunk(self._current_chunk)
        print("start compute coverage")
        self._genome.compute_coverage()
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
            attrs[attr] = val.strip("'").strip('"')
        if 'type' not in attrs:
            raise WigError('wiggle type is not present: {}.'.format(line))
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
            raise WigError("this data line '{}' is not preceded by declaration".format(line))
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



