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
import re
import collections
from abc import ABCMeta, abstractmethod
import logging

import psutil
import numpy as np


_log = logging.getLogger(__name__)


class WigError(Exception):
    """
    Handle error related to wig parsing 
    """
    pass


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


    @abstractmethod
    def is_fixed_step(self):
        """
        This is an abstract methods, must be implemented in inherited class
        :return: True if i's a fixed chunk of data, False otheweise
        :rtype: boolean
        """
        return NotImplemented

    @abstractmethod
    def parse_data_line(self, line, chrom, strand_type):
        """
        parse a line of data and append the results in the corresponding strand
        This is an abstract methods, must be implemented in inherited class.
        
        :param line: line of data to parse (the white spaces at the end must be strip)
        :type line: string
        :param chrom: the chromosome to add coverage data
        :type chrom: :class:`Chromosome` object.
        :param strand_type: which kind of wig is parsing: forward, reverse, or mixed strand
        :type strand_type: string '+' , '-', 'mixed' 
        """
        return NotImplemented

    @staticmethod
    def _convert_cov(strand_type, cov):
        if strand_type == 'mixed':
            cov = float(cov)
        elif strand_type == '-':
            cov = - abs(float(cov))
        elif strand_type == '+':
            cov = abs(float(cov))
        else:
            raise ValueError("value: '{}' is not allowed for strand parameter.".format(strand_type))
        return cov


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
        # we switch from 1-based positions in wig into 0-based position in chromosome
        # to have the same behavior as in bam
        self._current_pos = self.start - 1


    def is_fixed_step(self):
        """
        :return: True
        :rtype: boolean
        """
        return True


    def parse_data_line(self, line, chrom, strand_type):
        """
        parse line of data following a fixedStep Declaration.
        add the result on the corresponding strand (forward if coverage value is positive, reverse otherwise)
        :param line: line of data to parse (the white spaces at the end must be strip)
        :type line: string
        :param chrom: the chromosome to add coverage data
        :type chrom: :class:`Chromosome` object.
        :param strand_type: which kind of wig is parsing: forward, reverse, or mixed strand
        :type strand_type: string '+' , '-', 'mixed' 
        """
        # the line is already striped
        cov = [self._convert_cov(strand_type, line)] * self.span
        # in FixedChunk we translate the origin to a 0-based position at the __init__
        pos = self._current_pos
        chrom[pos:pos + self.span] = cov
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
    
    
    | for = [11.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    | rev = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 30.0, 0.0, 0.0, 0.0, 50.0, 50.0]
    """


    def is_fixed_step(self):
        """
        :return: False
        :rtype: boolean
        """
        return False


    def parse_data_line(self, line, chrom, strand_type):
        """
        Parse line of data following a variableStep Declaration.
        Add the result on the corresponding strand (forward if coverage value is positive, reverse otherwise)
       
        :param line: line of data to parse (the white spaces at the end must be strip)
        :type line: string
        :param chrom: the chromosome to add coverage data
        :type chrom: :class:`Chromosome` object.
        :param strand_type: which kind of wig is parsing: forward, reverse, or mixed strand
        :type strand_type: string '+' , '-', 'mixed'
        :raise ValueError: if strand_type is different than 'mixed', '-', '+'
        """
        pos, cov = line.split()
        # we switch from 1-based positions in wig into 0-based position in chromosome
        # to have the same behavior as in bam
        pos = int(pos) - 1
        cov = [self._convert_cov(strand_type, cov)] * self.span
        chrom[pos:pos + self.span] = cov


class Chromosome:
    """
    Handle chromosomes. A chromosome as a name and contains :class:`Chunk` objects 
    (forward and reverse)
    """

    def __init__(self, name, size=1000000):
        """
        
        :param name: 
        :type name: str
        :param size: 
        :type size: the default size of the chromosome. 
            Each time we try to set a value greater than the chromosome the chromosome size is doubled.
            This is to protect the machine against memory swapping if the user 
            provide a wig file with very big chromosomes.
        """
        self.name = name
        self._pid = os.getpid()
        # 30 is the memory used to allocated new array of shape (2,1)
        # it was empirically determined
        est_avail = self._estimate_memory(size, 30)
        if est_avail <= 0:
            raise MemoryError("Not enough memory to create new chromosome {}".format(self.name))
        self._coverage = np.full((2, size), 0.)

    def __len__(self):
        """
        :return: the actual length of the chromosome
        :rtype: int
        """
        return self._coverage.shape[1]


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
        if isinstance(pos, slice):
            if isinstance(value, collections.Iterable):
                if (pos.stop - pos.start) != len(value):
                    raise ValueError("can assign only iterable of same length of the slice")
            else:
                raise TypeError('can only assign an iterable')

            if value[0] < 0:
                strand = 1
                value = [abs(v) for v in value]
            else:
                strand = 0
            last_pos = pos.stop
        else:
            if value < 0:
                strand = 1
                value = abs(value)
            else:
                strand = 0
            last_pos = pos
        while last_pos >= self._coverage.shape[1]:
            self._extend(size=self._coverage.shape[1])
        self._coverage[strand, pos] = value


    def __getitem__(self, pos):
        """
        :param pos: a position or a slice (0 based)
                   if pos is a slice the left indice is excluded
        :return: the coverage at this position or corresponding to this slice.
        :rtype: a list of 2 list of float [[float,...],[float, ...]]
        :raise IndexError: if pos is not in coverage or one bound of slice is out the coverage
        """
        return self._coverage[:, pos].tolist()


    def _extend(self, size=1000000, fill=0.):
        """
        Extend this chromosome of the size size and fill with fill.
        :param size: the size (in bp) we want to increase the chromosome.
        :type size: int
        :param fill: the default value to fill the chromosome.
        :type fill: float or nan
        :raise MemoryError: if the chromosome extension could overcome the free memory.
        """
        # 10 is the memory used to horizontally extend an array with one col and 2 rows fill with 0.
        # it was empirically determined on linux gentoo plateform with python 3.4.5 and numpy 1.11.2
        est_avail = self._estimate_memory(size, 10)
        tot_k_size = self.__len__() + size
        if est_avail <= 0:
            for unit in ('', 'K', 'M', 'G'):
                if tot_k_size < 1000.0:
                    break
                tot_k_size /= 1000.0
            h_size = "{:.1f}{}bp".format(tot_k_size, unit)
            raise MemoryError("Not enough memory to extend chromosome"
                              " {} to {})".format(self.name, h_size))
        chunk = np.full((2, size), fill_value=fill)
        self._coverage = np.hstack((self._coverage, chunk))


    def _estimate_memory(self, col_nb, mem_per_col):
        """
        :param col_nb: the number of column of the new array or the extension 
        :type col_nb: int
        :param mem_per_col: the memory needed to create or extend an array with one col and 2 rows fill with 0.0  
        :type mem_per_col: int
        :return: the estimation of free memory available after creating or extending chromosome
        :rtype: int 
        """
        my_process = psutil.Process(self._pid)
        est_mem = (col_nb * mem_per_col) + my_process.memory_info().rss
        est_avail = psutil.virtual_memory().available - est_mem
        return est_avail


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
        elif isinstance(chrom, Chromosome):
            return chrom.name in self._chromosomes
        else:
            raise TypeError("'in <Genome>' requires string or Chromosome as left operand, not '{}'".format(
                            chrom.__class__.__name__))


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

    def __init__(self, mixed_wig='', for_wig='', rev_wig=''):
        """

        :param mixed_wig: The path of the wig file to parse.
                          The wig file code for the 2 strands:
                          
                             - The positive coverage values for the forward strand
                             - The negative coverage values for the reverse strand

                          This parameter is incompatible with for_wig and rev_wig parameter.
        :type mixed_wig: string
        :param for_wig: The path of the wig file to parse. 
                        The wig file code for forward strand only.
                        This parameter is incompatible with mixed_wig parameter.
        :type for_wig: string
        :param rev_wig: The path of the wig file to parse. 
                        The wig file code for reverse strand only.
                        This parameter is incompatible with mixed_wig parameter.
        :type rev_wig: string
        """
        if not any((mixed_wig, for_wig, rev_wig)):
            raise WigError("The path for one or two wig files must be specify")
        elif mixed_wig and any((for_wig, rev_wig)):
            raise WigError("Cannot specify the path for mixed wig and forward or reverse wig in same time")

        self.declaration_type_pattern = re.compile('fixedStep|variableStep')
        self.trackline_pattern = re.compile("""(\w+)=(".+?"|'.+?'|\S+)""")
        self.data_line_pattern = re.compile('^-?\d+(\s+-?\d+(\.\d+)?)?$')
        self._path_mixed = mixed_wig
        self._path_for = for_wig
        self._path_rev = rev_wig
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
        it means that the coverage is on the reverse strand. So some positions can appear twice
        in one block of declaration (what I call a chunk).

        :return: a Genome coverage corresponding to the wig files (mixed strand on one wig or two separate wig)
        :rtype: :class:`Genome` object
        """
        self._genome = Genome()
        if self._path_mixed:
            wig_paths = [(self._path_mixed, 'mixed')]
        else:
            wig_paths = []
            if self._path_for:
                wig_paths.append((self._path_for, '+'))
            if self._path_rev:
                wig_paths.append((self._path_rev, '-'))
        for path, strand_type in wig_paths:
            with open(path, 'r') as wig_file:
                for line in wig_file:
                    line = line.strip()
                    if self.is_data_line(line):
                        self.parse_data_line(line, strand_type)
                    elif self.is_declaration_line(line):
                        self.parse_declaration_line(line)
                    elif self.is_track_line(line):
                        self.parse_track_line(line, strand_type=strand_type)
                    elif not line or self.is_comment_line(line):
                        continue
                    else:
                        raise WigError("the line is malformed: {}".format(line))
        return self._genome


    def is_data_line(self, line):
        """

        :param line: line to parse.
        :return: True if it's a data line, False otherwise
        """
        return bool(re.match(self.data_line_pattern, line))


    def parse_data_line(self, line, strand_type):
        """
        :param line: line to parse. It must not a comment_line, neither a track line nor a declaration line.
        :type line: string
        :type strand_type: string '+' , '-', 'mixed'
        :raise ValueError: if strand_type is different than 'mixed', '-', '+'
        """
        if self._current_chunk is None:
            raise WigError("this data line '{}' is not preceded by declaration".format(line))
        self._current_chunk.parse_data_line(line, self._current_chrom, strand_type)


    def is_declaration_line(self, line):
        """
        A single line, beginning with one of the identifiers variableStep or fixedStep, followed by attribute/value pairs
        for instance: ::

          fixedStep chrom=chrI start=1 step=10 span=5
        
        :param line: line to parse.
        :type line: string
        :return: True if line is a declaration line. False otherwise.
        :rtype: boolean
        """
        return bool(re.match(self.declaration_type_pattern, line))


    def parse_declaration_line(self, line):
        """
        Get the corresponding chromosome create one if necessary, 
        and set the current_chunk and current_chromosome.

        :param line: line to parse. The method :meth:`is_declaration_line` must return True with this line.
        """
        _log.info("parsing : {}".format(line))

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
            chrom = Chromosome(chrom_name, size=750000)
            self._genome.add(chrom)
        self._current_chrom = chrom


    @staticmethod
    def is_track_line(line):
        """
        A track line begins with the identifier track and followed by attribute/value pairs 
        for instance: ::

          track type=wiggle_0 name="fixedStep" description="fixedStep format" visibility=full autoScale=off
        
        :param line: line to parse.
        :type line: string
        :return: True if line is a track line. False otherwise.
        :rtype: boolean
        """
        return line.startswith('track')


    def parse_track_line(self, line, strand_type=''):
        """
        fill the genome infos with the information found on the track.

        :param line: line to parse. The method :meth:`is_track_line` must return True with this line.
        """
        _log.info('parsing : {}'.format(line))
        fields = re.findall(self.trackline_pattern, line)
        attrs = {}
        for attr, val in fields:
            attrs[attr] = val.strip("'").strip('"')
        if 'type' not in attrs:
            raise WigError('wiggle type is not present: {}.'.format(line))
        else:
            if strand_type == '+':
                self._genome.infos['forward'] = attrs
            elif strand_type == '-':
                self._genome.infos['reverse'] = attrs
            else:
                self._genome.infos = attrs

    @staticmethod
    def is_comment_line(line):
        """
        :param line: line to parse.
        :type line: string
        :return: True if line is a comment line. False otherwise.
        :rtype: boolean
        """
        return line.startswith('#')








