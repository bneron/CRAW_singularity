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



try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.wig import WigException, ExpandedStrand, FixedChunk, VariableChunk, Strand, Chromosome, Genome, WigParser


class TestExpandedStrand(CRAWTest):

    def test_coverage(self):
        es = ExpandedStrand(1, 5, 3)
        self.assertListEqual(es.coverage, [0, 0, 0, 0, 0, 0, 0])
        es = ExpandedStrand(1, 5, 1)
        self.assertListEqual(es.coverage, [0, 0, 0, 0, 0])

    def test_len(self):
        es = ExpandedStrand(1, 5, 3)
        self.assertEqual(len(es), 7)

    def test_stop(self):
        es = ExpandedStrand(1, 5, 3)
        self.assertEqual(es.stop, 7)

    def test_start(self):
        es = ExpandedStrand(1, 5, 3)
        self.assertEqual(es.start, 1)

    def test_eq(self):
        es1 = ExpandedStrand(1, 5, 3)
        es2 = ExpandedStrand(1, 5, 3)
        self.assertEqual(es1, es2)


class TestFixedChunk(CRAWTest):

    def test_is_fixed_step(self):
        kwargs = {"chrom": "chr3", "start": "400601"}
        fx_ch = FixedChunk(**kwargs)
        self.assertEqual(fx_ch.chrom, kwargs["chrom"])
        self.assertEqual(fx_ch.start, int(kwargs["start"]))
        self.assertEqual(fx_ch.step, 1)
        self.assertEqual(fx_ch.span, 1)
        kwargs = {"chrom": "chr3", "start": "400601", "step": "100"}
        fx_ch = FixedChunk(**kwargs)
        self.assertEqual(fx_ch.chrom, kwargs["chrom"])
        self.assertEqual(fx_ch.start, int(kwargs["start"]))
        self.assertEqual(fx_ch.step, int(kwargs["step"]))
        self.assertEqual(fx_ch.span, 1)
        kwargs = {"chrom": "chr3", "start": "400601", "step": "100", "span": "5"}
        fx_ch = FixedChunk(**kwargs)
        self.assertEqual(fx_ch.chrom, kwargs["chrom"])
        self.assertEqual(fx_ch.start, int(kwargs["start"]))
        self.assertEqual(fx_ch.step, int(kwargs["step"]))
        self.assertEqual(fx_ch.span, int(kwargs["span"]))

        kwargs = {"start": "400601", "step": "100", "span": "5"}
        with self.assertRaises(WigException) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'chrom' field  is not present.")

        kwargs = {"chrom": "chr3", "start": "400601", "step": "100", "span": "0"}
        with self.assertRaises(WigException) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'0' is not allowed as span value.")

        kwargs = {"chrom": "chr3", "step": "100", "span": "1"}
        with self.assertRaises(WigException) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'start' must be defined for 'fixedStep'.")


    def test_is_fixed_step(self):
        kwargs = {"chrom": "chr3", "start": "400601", "step": "100", "span": "5"}
        fx_ch = FixedChunk(**kwargs)
        self.assertTrue(fx_ch.is_fixed_step())

    def test_parse_data_line(self):
        lines = ("11", "22", "30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l)
        self.assertListEqual(fx_ch.forward, [(10, float(11)), (20, float(22)), (30, float(30)), (40, float(50))])
        self.assertListEqual(fx_ch.reverse, [])

    def test_expand(self):
        lines = ("11", "22", "30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l)

        forward, reverse = fx_ch.expand()
        expected_for = ExpandedStrand(10, 40, 2)
        expected_for[10] = 11.0
        expected_for[11] = 11.0
        expected_for[20] = 22.0
        expected_for[21] = 22.0
        expected_for[30] = 30.0
        expected_for[31] = 30.0
        expected_for[40] = 50.0
        expected_for[41] = 50.0
        expected_rev = ExpandedStrand(0, 0, 0)
        self.assertEqual(forward, expected_for)
        self.assertEqual(reverse, expected_rev)


class TestVariableChunk(CRAWTest):

    def test_is_fixed_step(self):
        kwargs = {"chrom": "chr3"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.step, 1)
        self.assertEqual(var_ch.span, 1)
        kwargs = {"chrom": "chr3", "step": "100"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.step, int(kwargs["step"]))
        self.assertEqual(var_ch.span, 1)

        kwargs = {"chrom": "chr3", "span": "5"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.span, int(kwargs["span"]))

        kwargs = {"span": "5"}
        with self.assertRaises(WigException) as ctx:
            VariableChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'chrom' field  is not present.")

        kwargs = {"chrom": "chr3", "span": "0"}
        with self.assertRaises(WigException) as ctx:
            VariableChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'0' is not allowed as span value.")


    def test_is_fixed_step(self):
        kwargs = {"chrom": "chr3", "span": "5"}
        var_ch = VariableChunk(**kwargs)
        self.assertFalse(var_ch.is_fixed_step())

    def test_parse_data_line(self):
        lines = ("10 11", "20 22", "30 -30", "40 -50")
        kwargs = {"chrom": "chr3", "span": "2"}
        var_ch = VariableChunk(**kwargs)
        for l in lines:
            var_ch.parse_data_line(l)
        self.assertListEqual(var_ch.forward, [(10, float(11)), (20, float(22))])
        self.assertListEqual(var_ch.reverse, [(30, float(30)), (40, float(50))])

    def test_expand(self):
        lines = ("10 11", "20 22", "20 -30", "40 -50")
        kwargs = {"chrom": "chr3", "span": "2"}
        var_ch = VariableChunk(**kwargs)
        for l in lines:
            var_ch.parse_data_line(l)

        forward, reverse = var_ch.expand()
        expected_for = ExpandedStrand(10, 20, 2)
        expected_for[10] = 11.0
        expected_for[11] = 11.0
        expected_for[20] = 22.0
        expected_for[21] = 22.0

        expected_rev = ExpandedStrand(20, 40, 2)
        expected_rev[20] = 30.0
        expected_rev[21] = 30.0
        expected_rev[40] = 50.0
        expected_rev[41] = 50.0
        self.assertEqual(forward, expected_for)
        self.assertEqual(reverse, expected_rev)


class TestStrand(CRAWTest):

    def test_get_coverage(self):
        st = Strand()
        st.coverage = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.assertListEqual(st.get_coverage(1, 4), [1, 2, 3, 4])
        self.assertListEqual(st.get_coverage(-1, 4), [0, 0, 1, 2, 3, 4])
        self.assertListEqual(st.get_coverage(5, 12), [5, 6, 7, 8, 9, 0, 0])
        self.assertListEqual(st.get_coverage(-2, 12), [0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0])

    def test_eq(self):
        st1 = Strand()
        st2 = Strand()
        self.assertEqual(st1, st2)
        st1.coverage = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        st2.coverage = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.assertEqual(st1, st2)
        st2.coverage = [0, 0]
        self.assertNotEqual(st1, st2)

    def test_len(self):
        st1 = Strand()
        self.assertEqual(len(st1), 0)
        st1.coverage = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.assertEqual(len(st1), 10)

    def extend(self):
        st1 = Strand()
        ext_strand = ExpandedStrand(0, 4, 1)
        st1.extend(ext_strand)
        self.assertListEqual(st1.coverage, ext_strand.coverage)


class TestChromosome(CRAWTest):

    def test_Chromosome(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        self.assertEqual(ch.name, ch_name)

    def test_iadd(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)

        kwargs = {"chrom": ch_name, "start": "10", "step": "10"}
        fx_ck = FixedChunk(**kwargs)
        lines = ("11", "22")
        for l in lines:
            fx_ck.parse_data_line(l)
        ch += fx_ck
        expected_forward = Strand()
        expected_forward.coverage = ([0] * 9) + [11.] + ([0] * 9) + [22.]
        self.assertEqual(ch.forward, expected_forward)

        ch_name = 'ChrIII'
        span = 2
        ch = Chromosome(ch_name)
        lines = ("10 11", "20 22", "20 -30", "30 -50")
        kwargs = {"chrom": ch_name, "span": str(span)}
        var_ch = VariableChunk(**kwargs)
        for l in lines:
            var_ch.parse_data_line(l)
        ch += var_ch
        expected_forward = Strand()
        expected_forward.coverage = ([0] * 9) + ([11.] * span) + ([0] * 8) + ([22.] * span)
        expected_reverse = Strand()
        expected_reverse.coverage = ([0] * 19) + ([30.] * span) + ([0] * 8) + ([40.] * span)
        self.assertEqual(ch.forward, expected_forward)


class TestGenome(CRAWTest):

    def test_add_and_get(self):
        genome = Genome()
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        genome.add(ch)
        self.assertEqual(ch, genome[ch_name])
        ch2 = Chromosome(ch_name)
        genome.add(ch2)
        self.assertNotEqual(ch, genome[ch_name])
        self.assertEqual(ch2, genome[ch_name])

    def test_membership(self):
        genome = Genome()
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        self.assertFalse(ch in genome)
        genome.add(ch)
        self.assertTrue(ch in genome)


