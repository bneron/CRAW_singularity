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
import logging

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.wig import WigError, Coverage, FixedChunk, VariableChunk, Chromosome, Genome, WigParser, _log


class TestCoverage(CRAWTest):

    def test_coverage_init(self):
        cov = Coverage(1, 5, 3)
        self.assertListEqual(cov.coverage, [0.] * 7)
        cov = Coverage(1, 5, 1)
        self.assertListEqual(cov.coverage, [0.] * 5)

    def test_len(self):
        cov = Coverage(1, 5, 3)
        self.assertEqual(len(cov), 7)

    def test_stop(self):
        cov = Coverage(1, 5, 3)
        self.assertEqual(cov.stop, 7)

    def test_start(self):
        cov = Coverage(1, 5, 3)
        self.assertEqual(cov.start, 1)

    def test_eq(self):
        cov1 = Coverage(1, 5, 3)
        cov2 = Coverage(1, 5, 3)
        self.assertEqual(cov1, cov2)

    def test_setitem(self):
        cov = Coverage(1, 5, 3)
        pos = 2
        value = 3.0
        cov[pos] = value
        expect_cov = [0.] * 7
        expect_cov[pos - 1] = 3.0
        self.assertEqual(cov.coverage, expect_cov)
        cov = Coverage(1, 5, 3)
        cov[pos:pos + 3] = [3.] * 3
        expect_cov = [0.] * 7
        expect_cov[pos - 1: pos + 2] = [3.0, 3.0, 3.0]
        self.assertEqual(cov.coverage, expect_cov)

        with self.assertRaises(TypeError) as ctx:
            cov[cov.start:cov.stop] = 3.0
        self.assertEqual(str(ctx.exception), 'can only assign an iterable')

        with self.assertRaises(ValueError) as ctx:
            cov[cov.start:cov.stop] = (3.0, 2.0)
        self.assertEqual(str(ctx.exception), 'can assign only iterable of same length of the slice')
        with self.assertRaises(IndexError) as ctx:
            cov[cov.stop + 1] = 3.0
        self.assertEqual(str(ctx.exception), 'Coverage assignment position 8 out of range[1:7]')

    def test_getitem(self):
        cov = Coverage(1, 5, 3)
        pos = 2
        value = 3.0
        cov[pos] = value
        self.assertEqual(cov[pos], value)
        self.assertEqual(cov[pos - 1:pos + 2], [0.0, value, 0.0])
        with self.assertRaises(IndexError) as ctx:
            _ = cov[cov.stop + 2]
        self.assertEqual(str(ctx.exception), 'Coverage  position 9 out of range[1:7]')
        with self.assertRaises(IndexError) as ctx:
            _ = cov[cov.start:cov.stop + 2]
        self.assertEqual(str(ctx.exception), 'Coverage  position [1:9] out of range [1:7]')

    def test_str(self):
        cov = Coverage(1, 5, 3)
        expec_str = """start= 1; stop= 7
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]"""
        self.assertEqual(str(cov), expec_str)

    def test_join(self):
        cov1 = Coverage(1, 5, 3)
        cov1[2:4] = [2., 3.]
        cov2 = Coverage(10, 15, 3)
        cov2[12:14] = [12.,13.]
        cov = Coverage.join((cov1, cov2))
        expec_cov = [0.] * 17
        expec_cov[1:3] = [2., 3.]
        expec_cov[11:13] = [12., 13.]
        self.assertEqual(cov.coverage, expec_cov)
        with self.assertRaises(ValueError) as ctx:
            _ = Coverage.join([])
        self.assertEqual(str(ctx.exception), 'No coverage object to join')

class TestFixedChunk(CRAWTest):

    @classmethod
    def setUpClass(cls):
        _log.setLevel(logging.ERROR)

    def test_stop(self):
        lines = ("11", "22", "30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l)
        self.assertEqual(fx_ch.stop, 41)

    def test_fixed_step(self):
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
        with self.assertRaises(WigError) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'chrom' field  is not present.")

        kwargs = {"chrom": "chr3", "start": "400601", "step": "100", "span": "0"}
        with self.assertRaises(WigError) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'0' is not allowed as span value.")

        kwargs = {"chrom": "chr3", "step": "100", "span": "1"}
        with self.assertRaises(WigError) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'start' must be defined for 'fixedStep'.")

        kwargs = {"chrom": "chr3", "start": "400601", "step": "1", "span": "5"}
        with self.assertRaises(WigError) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'span' cannot be greater than 'step'.")

        kwargs = {"chrom": "chr3", "start": "400601", "step": "0", "span": "5"}
        with self.assertRaises(WigError) as ctx:
            FixedChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'step' must be strictly positive.")


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

        lines = ("-11", "22", "-30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l)
        self.assertListEqual(fx_ch.forward, [(20, float(22)), (40, float(50))])
        self.assertListEqual(fx_ch.reverse, [(10, float(11)), (30, float(30))])

    def test_to_coverages(self):
        lines = ("11", "22", "30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l)

        forward, reverse = fx_ch.to_coverages()
        expected_for = Coverage(10, 40, 2)
        expected_for[10] = 11.0
        expected_for[11] = 11.0
        expected_for[20] = 22.0
        expected_for[21] = 22.0
        expected_for[30] = 30.0
        expected_for[31] = 30.0
        expected_for[40] = 50.0
        expected_for[41] = 50.0
        expected_rev = Coverage(10, 40, 2)
        self.assertEqual(forward, expected_for)
        self.assertEqual(reverse, expected_rev)


class TestVariableChunk(CRAWTest):

    def test_stop(self):
        lines = ("10 11", "20 22", "30 -30", "40 -50")
        kwargs = {"chrom": "chr3", "span": "2"}
        var_ch = VariableChunk(**kwargs)
        for l in lines:
            var_ch.parse_data_line(l)
        self.assertEqual(var_ch.stop, 41)

    def test_variable_step(self):
        kwargs = {"chrom": "chr3"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.span, 1)
        kwargs = {"chrom": "chr3", "step": "100"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.span, 1)

        kwargs = {"chrom": "chr3", "span": "5"}
        var_ch = VariableChunk(**kwargs)
        self.assertEqual(var_ch.chrom, kwargs["chrom"])
        self.assertEqual(var_ch.span, int(kwargs["span"]))

        kwargs = {"span": "5"}
        with self.assertRaises(WigError) as ctx:
            VariableChunk(**kwargs)
        self.assertEqual(str(ctx.exception), "'chrom' field  is not present.")

        kwargs = {"chrom": "chr3", "span": "0"}
        with self.assertRaises(WigError) as ctx:
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

    def test_to_coverages(self):
        lines = ("10 11", "20 22", "20 -30", "40 -50")
        kwargs = {"chrom": "chr3", "span": "2"}
        var_ch = VariableChunk(**kwargs)
        for l in lines:
            var_ch.parse_data_line(l)

        forward, reverse = var_ch.to_coverages()
        expected_for = Coverage(10, 40, 2)
        expected_for[10] = 11.0
        expected_for[11] = 11.0
        expected_for[20] = 22.0
        expected_for[21] = 22.0

        expected_rev = Coverage(10, 40, 2)
        expected_rev[20] = 30.0
        expected_rev[21] = 30.0
        expected_rev[40] = 50.0
        expected_rev[41] = 50.0
        self.assertEqual(forward, expected_for)
        self.assertEqual(reverse, expected_rev)


class TestChromosome(CRAWTest):

    def test_Chromosome(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        self.assertEqual(ch.name, ch_name)

    def test_add_chunk(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)

        kwargs = {"chrom": ch_name, "start": "10", "step": "10"}
        fx_ck = FixedChunk(**kwargs)
        lines = ("11", "22")
        for l in lines:
            fx_ck.parse_data_line(l)
        ch.add_chunk(fx_ck)
        self.assertListEqual(ch._chunks, [fx_ck])

        span = 2
        lines = ("10 11", "20 22", "20 -30", "30 -50")
        kwargs = {"chrom": ch_name, "span": str(span)}
        var_ck = VariableChunk(**kwargs)
        for l in lines:
            var_ck.parse_data_line(l)
        ch.add_chunk(var_ck)
        self.assertListEqual(ch._chunks, [fx_ck, var_ck])


    def test_get_coverage(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)

        kwargs = {"chrom": ch_name, "start": "10", "step": "10", "span": "3"}
        fx_ck = FixedChunk(**kwargs)
        lines = ("11", "22", "33")
        for l in lines:
            fx_ck.parse_data_line(l)
        ch.add_chunk(fx_ck)

        span = 2
        lines = ("40 11", "42 22", "40 -30", "42 -50")
        kwargs = {"chrom": ch_name, "span": str(span)}
        var_ck = VariableChunk(**kwargs)
        for l in lines:
            var_ck.parse_data_line(l)
        ch.add_chunk(var_ck)

        fwd = [0.] * 50
        fwd[10:10 + 3] = [11.] * 3
        fwd[20:20 + 3] = [22.] * 3
        fwd[30:30 + 3] = [33.] * 3
        fwd[40:40 + 2] = [11.] * 2
        fwd[42:42 + 2] = [22.] * 2

        rev = [0.] * 50
        rev[40:40 + 2] = [30.] * 2
        rev[42:42 + 2] = [50.] * 2

        # start and stop are in first chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(10, 20)
        self.assertEqual(recv_fwd_cov, fwd[10:21])
        self.assertEqual(recv_rev_cov, rev[10:21])

        # start is before chunk start stop in first chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(1, 20)
        self.assertEqual(recv_fwd_cov, fwd[1:21])
        self.assertEqual(recv_rev_cov, rev[1:21])

        # start and stop are in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(40, 42)
        self.assertEqual(recv_fwd_cov, fwd[40:43])
        self.assertEqual(recv_rev_cov, rev[40:43])

        # stop is beyond in 2nd chunk stop
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(40, 45)
        self.assertEqual(recv_fwd_cov, fwd[40:46])
        self.assertEqual(recv_rev_cov, rev[40:46])

        # start is in first chunk stop is in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(11, 42)
        self.assertEqual(recv_fwd_cov, fwd[11:43])
        self.assertEqual(recv_rev_cov, rev[11:43])

        # start is before first chunk stop is beyond 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(1, 45)
        self.assertEqual(recv_fwd_cov, fwd[1:46])
        self.assertEqual(recv_rev_cov, rev[1:46])

        # start is between first chunk stop is in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(35, 42)
        self.assertEqual(recv_fwd_cov, fwd[35:43])
        self.assertEqual(recv_rev_cov, rev[35:43])

        # start and stop are beyond last chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(50, 52)
        self.assertEqual(recv_fwd_cov, [0.] * 3)
        self.assertEqual(recv_rev_cov, [0.] * 3)

        # start and stop are before first chunk
        recv_fwd_cov, recv_rev_cov = ch.get_coverage(1, 3)
        self.assertEqual(recv_fwd_cov, [0.] * 3)
        self.assertEqual(recv_rev_cov, [0.] * 3)


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

        with self.assertRaises(TypeError) as ctx:
            genome.add(3)
        self.assertEqual(str(ctx.exception), "Genome can contains only Chromosome objects")


    def test_membership(self):
        genome = Genome()
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        self.assertFalse(ch in genome)
        genome.add(ch)
        self.assertTrue(ch in genome)
        self.assertTrue(ch_name in genome)
        with self.assertRaises(TypeError) as ctx:
            3 in genome
        self.assertEqual(str(ctx.exception), "'in <Genome>' requires string or Chromosome as left operand, not 'int'")


class TestWigParser(CRAWTest):

    def test_is_track_line(self):
        wip_p = WigParser('toto')
        line = 'track type=wiggle_0 name=BCMSolidWeCoca48PatientCoverageclean viewLimits=0:1'
        self.assertTrue(wip_p.is_track_line(line))
        line = 'type=wiggle_0 name=BCMSolidWeCoca48PatientCoverageclean viewLimits=0:1'
        self.assertFalse(wip_p.is_track_line(line))


    def test_is_declaration_line(self):
        wip_p = WigParser('toto')
        line = 'fixedStep chrom=chr1 start=58951 step=1'
        self.assertTrue(wip_p.is_declaration_line(line))
        line = 'variableStep chrom=chrI span=1'
        self.assertTrue(wip_p.is_declaration_line(line))
        line = 'undefineStep chrom=chrI span=1'
        self.assertFalse(wip_p.is_declaration_line(line))


    def test_is_comment_line(self):
        wip_p = WigParser('toto')
        line = '#fixedStep chrom=chr1 start=58951 step=1'
        self.assertTrue(wip_p.is_comment_line(line))
        line = 'fixedStep chrom=chr1 start=58951 step=1'
        self.assertFalse(wip_p.is_comment_line(line))


    def test_parse_data_line(self):
        wip_p = WigParser('toto')
        kwargs = {"chrom": "chr3", "start": "10", "step": "100", "span": "5"}
        wip_p._current_chunk = FixedChunk(**kwargs)
        wip_p.parse_data_line("3")
        self.assertEqual(wip_p._current_chunk.forward[-1],  (10, 3.0))
        wip_p.parse_data_line("5")
        self.assertEqual(wip_p._current_chunk.forward[-1],  (110, 5.0))

        kwargs = {"chrom": "chr3", "span": "5"}
        wip_p._current_chunk = VariableChunk(**kwargs)
        wip_p.parse_data_line("10 3")
        self.assertEqual(wip_p._current_chunk.forward[-1],  (10, 3.0))
        wip_p.parse_data_line("10 -5")
        self.assertEqual(wip_p._current_chunk.reverse[-1],  (10, 5.0))

        wip_p = WigParser('toto')
        with self.assertRaises(WigError) as ctx:
            wip_p.parse_data_line("3")
        self.assertEqual(str(ctx.exception), "this data line '3' is not preceded by declaration")


    def test_parse_declaration_line(self):
        wip_p = WigParser('toto')
        wip_p._genome = Genome()
        line = 'fixedStep chrom=chr1 start=10 step=1'
        wip_p.parse_declaration_line(line)
        self.assertTrue(isinstance(wip_p._current_chunk, FixedChunk))
        self.assertTrue('chr1' in wip_p._genome)
        self.assertTrue(wip_p._current_chunk.start, 10)
        self.assertTrue(wip_p._current_chunk.step, 1)
        self.assertTrue(wip_p._current_chunk.span, 1)
        self.assertTrue(wip_p._current_chrom is wip_p._genome['chr1'])

        wip_p = WigParser('toto')
        wip_p._genome = Genome()
        line = 'variableStep chrom=chrI span=2'
        wip_p.parse_declaration_line(line)
        self.assertTrue(isinstance(wip_p._current_chunk, VariableChunk))
        self.assertTrue('chrI' in wip_p._genome)
        self.assertTrue(wip_p._current_chunk.span, 2)
        self.assertTrue(wip_p._current_chrom is wip_p._genome['chrI'])


    def test_parse_fixed_wig(self):
        expected_forward = [0] * 204
        span = 5
        expected_forward[0: 5] = [1.] * span
        for i in range(1, 5):
            pos = (i * 10)
            expected_forward[pos: pos + span] = [float(i + 1)] * span
        for i, pos in enumerate(range(99, 140, 10), 1):
            expected_forward[pos] = float(i)
        for i , pos in enumerate(range(199, 204), 1):
            expected_forward[pos] = float(i)

        wig_parser = WigParser(os.path.join(self._data_dir, 'wig_fixed.wig'))
        genome = wig_parser.parse()

        self.assertTrue('chrI' in genome)
        chrI = genome['chrI']
        recieved_forward = chrI.forward
        self.assertListEqual(recieved_forward.coverage, expected_forward)
        self.assertListEqual(chrI.reverse.get_coverage, [])


    def test_parse_variable_wig(self):
        expec_chrI_forward = [0, 0, 0, 6., 7., 8., 10., 11.]
        expec_chrI_reverse = [1., 2., 3., 4., 5.]

        expec_chrII_forward = [0] * 69
        expec_chrII_forward += [1., 1., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_forward += [2., 2., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_forward += [3., 3.]

        expec_chrII_reverse = [0] * 9
        expec_chrII_reverse += [1., 1., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_reverse += [2., 2., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_reverse += [3., 3., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_reverse += [4., 4., 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_reverse += [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        expec_chrII_reverse += [5., 5.]

        wig_parser = WigParser(os.path.join(self._data_dir, 'wig_variable.wig'))
        genome = wig_parser.parse()

        self.assertTrue('chrI' in genome)
        chrI = genome['chrI']
        self.assertListEqual(chrI.forward.get_coverage, expec_chrI_forward)
        self.assertListEqual(chrI.reverse.get_coverage, expec_chrI_reverse)

        self.assertTrue('chrII' in genome)
        chrII = genome['chrII']
        self.assertListEqual(chrII.forward.get_coverage, expec_chrII_forward)
        self.assertListEqual(chrII.reverse.get_coverage, expec_chrII_reverse)
