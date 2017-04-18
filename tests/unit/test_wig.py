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

import numpy as np

try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.wig import WigError, FixedChunk, VariableChunk, Chromosome, Genome, WigParser, _log




class TestFixedChunk(CRAWTest):

    @classmethod
    def setUpClass(cls):
        _log.setLevel(logging.ERROR)


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
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)

        for l in lines:
            fx_ch.parse_data_line(l, ch)
        exp_cov = np.full((2, 50), 0.)
        exp_cov[0, 9:11] = [11] * 2
        exp_cov[0, 19:21] = [22] * 2
        exp_cov[0, 29:31] = [30] * 2
        exp_cov[0, 39:41] = [50] * 2
        got_forward, got_reverse = ch[:50]
        print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("got_forward", got_forward)
        print("exp_forward",exp_cov[0].tolist())
        print("===================================")
        print("got_reverse", got_reverse)
        print("exp_reverse",exp_cov[1].tolist())
        self.assertListEqual(exp_cov[0].tolist(), got_forward)
        self.assertListEqual(exp_cov[1].tolist(), got_reverse)

        lines = ("-11", "22", "-30", "50")
        kwargs = {"chrom": "chr3", "start": "10", "step": "10", "span": "2"}
        fx_ch = FixedChunk(**kwargs)
        for l in lines:
            fx_ch.parse_data_line(l, ch)
        exp_cov[1, 9:11] = [11] * 2
        exp_cov[0, 19:21] = [22] * 2
        exp_cov[1, 29:31] = [30] * 2
        exp_cov[0, 39:41] = [50] * 2
        self.assertListEqual(exp_cov[0].tolist(), got_forward)
        self.assertListEqual(exp_cov[1].tolist(), got_reverse)



class TestVariableChunk(CRAWTest):


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
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        for l in lines:
            var_ch.parse_data_line(l, ch)
        self.assertListEqual(var_ch.forward, [(10, float(11)), (20, float(22))])
        self.assertListEqual(var_ch.reverse, [(30, float(30)), (40, float(50))])


class TestChromosome(CRAWTest):

    def test_Chromosome(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        self.assertEqual(ch.name, ch_name)


    def test_get_item(self):
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)

        kwargs = {"chrom": ch_name, "start": "10", "step": "10", "span": "3"}
        fx_ck = FixedChunk(**kwargs)
        lines = ("11", "22", "33")
        for l in lines:
            fx_ck.parse_data_line(l, ch)

        span = 2
        lines = ("40 11", "42 22", "40 -30", "42 -50")
        kwargs = {"chrom": ch_name, "span": str(span)}
        var_ck = VariableChunk(**kwargs)
        for l in lines:
            var_ck.parse_data_line(l, ch)

        fwd = [0.] * 50
        fwd[9:9 + 3] = [11.] * 3
        fwd[19:19 + 3] = [22.] * 3
        fwd[29:29 + 3] = [33.] * 3
        fwd[39:39 + 2] = [11.] * 2
        fwd[41:41 + 2] = [22.] * 2

        rev = [0.] * 50
        rev[39:39 + 2] = [30.] * 2
        rev[41:41 + 2] = [50.] * 2

        # start and stop are in first chunk
        recv_fwd_cov, recv_rev_cov = ch[9:20]
        self.assertEqual(recv_fwd_cov, fwd[9:20])
        self.assertEqual(recv_rev_cov, rev[9:20])

        # start is before chunk start stop in first chunk
        recv_fwd_cov, recv_rev_cov = ch[1: 20]
        self.assertEqual(recv_fwd_cov, fwd[1:20])
        self.assertEqual(recv_rev_cov, rev[1:20])

        # start and stop are in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch[40: 42]
        self.assertEqual(recv_fwd_cov, fwd[40:42])
        self.assertEqual(recv_rev_cov, rev[40:42])

        # stop is beyond the 2nd chunk stop
        recv_fwd_cov, recv_rev_cov = ch[40:45]
        self.assertEqual(recv_fwd_cov, fwd[40:45])
        self.assertEqual(recv_rev_cov, rev[40:45])

        # start is in first chunk, stop is in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch[11:42]
        self.assertEqual(recv_fwd_cov, fwd[11:42])
        self.assertEqual(recv_rev_cov, rev[11:42])

        # start is before first chunk stop is beyond 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch[1: 45]
        self.assertEqual(recv_fwd_cov, fwd[1:45])
        self.assertEqual(recv_rev_cov, rev[1:45])

        # start is between first and 2nd chunk, stop is in 2nd chunk
        recv_fwd_cov, recv_rev_cov = ch[34: 42]
        self.assertEqual(recv_fwd_cov, fwd[34:42])
        self.assertEqual(recv_rev_cov, rev[34:42])

        # start and stop are beyond last chunk
        recv_fwd_cov, recv_rev_cov = ch[50:53]
        self.assertEqual(recv_fwd_cov, [0.] * 3)
        self.assertEqual(recv_rev_cov, [0.] * 3)

        # start and stop are before first chunk
        recv_fwd_cov, recv_rev_cov = ch[0:3]
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


    def test_del(self):
        genome = Genome()
        ch_name = 'ChrII'
        ch = Chromosome(ch_name)
        genome.add(ch)
        self.assertEqual(ch, genome[ch_name])
        del genome[ch_name]
        self.assertFalse(ch_name in genome)


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


    def test_chromosomes(self):
        genome = Genome()
        chromosomes = [Chromosome('chrI'), Chromosome('chrII')]
        for ch in chromosomes:
            genome.add(ch)
        self.assertSetEqual(set(chromosomes), set(genome.chromosomes))


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

    def test_parse_track_line(self):
        wip_p = WigParser('toto')
        line = 'track type=wiggle_0 name=BCMSolidWeCoca48PatientCoverageclean viewLimits=0:1'
        kwargs = {"chrom": "chr3", "start": "10", "step": "100", "span": "5"}
        wip_p._genome = Genome()
        wip_p._current_chunk = FixedChunk(**kwargs)
        infos = {'type': 'wiggle_0',
                 'name': 'BCMSolidWeCoca48PatientCoverageclean',
                 'viewLimits': '0:1'}
        wip_p.parse_track_line(line)
        self.assertDictEqual(wip_p._genome.infos, infos)

        wip_p._current_chunk = FixedChunk(**kwargs)
        line = 'track name=BCMSolidWeCoca48PatientCoverageclean viewLimits=0:1'
        with self.assertRaises(WigError) as ctx:
            wip_p.parse_track_line(line)
        self.assertEqual(str(ctx.exception), 'wiggle type is not present: {}.'.format(line))

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
        expected_forward = [0.] * 260
        expected_reverse = [0.] * 260
        span = 5
        for i, pos in enumerate(range(0, 50, 10), 1):
            expected_forward[pos: pos + span] = [float(i)] * span
        for i, pos in enumerate(range(99, 140, 10), 1):
            expected_forward[pos] = float(i)
        for i, pos in enumerate(range(199, 204), 1):
            expected_reverse[pos] = float(i)

        wig_parser = WigParser(os.path.join(self._data_dir, 'wig_fixed.wig'))
        genome = wig_parser.parse()

        self.assertTrue('chrI' in genome)
        chrI = genome['chrI']
        received_forward, received_reverse = chrI.get_coverage(0, 52)
        self.assertListEqual(received_forward, expected_forward[0:52])
        self.assertListEqual(received_reverse, [0.] * 52)

        received_forward, received_reverse = chrI.get_coverage(190, 252)
        self.assertListEqual(received_forward, expected_forward[190:252])
        self.assertListEqual(received_reverse, expected_reverse[190:252])

    def test_parse_variable_wig(self):
        # chrI position 1->5 on rev 4->8 on fwd
        expec_chrI_forward = [0., 0., 0., 6., 7., 8., 10., 11., 0., 0.]
        expec_chrI_reverse = [1., 2., 3., 4., 5., 0., 0., 0., 0., 0.]

        expec_chrII_forward = [0.] * 100
        expec_chrII_reverse = [0.] * 100
        span = 2
        expec_chrII_forward[69:11] = [1.] * span
        expec_chrII_forward[79:21] = [2.] * span
        expec_chrII_forward[89:91] = [3.] * span

        expec_chrII_reverse[9:11] = [1.] * span
        expec_chrII_reverse[19:21] = [2.] * span
        expec_chrII_reverse[29:31] = [3.] * span
        expec_chrII_reverse[39:41] = [4.] * span
        expec_chrII_reverse[59:61] = [5.] * span

        wig_parser = WigParser(os.path.join(self._data_dir, 'wig_variable.wig'))
        genome = wig_parser.parse()

        self.assertTrue('chrI' in genome)
        chrI = genome['chrI']
        recv_chrI_forward, recv_chrI_reverse = chrI.get_coverage(0, 10)
        self.assertListEqual(recv_chrI_forward, expec_chrI_forward)
        self.assertListEqual(recv_chrI_reverse, expec_chrI_reverse)

        self.assertTrue('chrII' in genome)
        chrII = genome['chrII']
        recv_chrII_forward, recv_chrII_reverse = chrII.get_coverage(0, 82)
        # get coverage start and stop are included, and numbered from 1
        self.assertListEqual(recv_chrII_forward, expec_chrII_forward[:82])
        self.assertListEqual(recv_chrII_reverse, expec_chrII_reverse[:82])

    def test_parse(self):
        wig_p = WigParser(os.path.join(self._data_dir, 'wig_fixed_w_comment.wig'))
        genome = wig_p.parse()
        infos = {'type': 'wiggle_0',
                 'name': "wig de test comment line",
                 'color': '96,144,246',
                 'altColor': '96,144,246',
                 'autoScale': 'on',
                 'graphType': 'bar'}
        self.assertDictEqual(genome.infos, infos)
        ch_name = 'chrI'
        ch = Chromosome(ch_name)
        genome.add(ch)
        kwargs = {"chrom": ch_name, "start": "1", "step": "10", "span": "5"}
        fx_ck1 = FixedChunk(**kwargs)
        lines = ("1", "2", "3", "4", "5")
        for l in lines:
            fx_ck1.parse_data_line(l)
        ch.add_chunk(fx_ck1)
        kwargs = {"chrom": ch_name, "start": "100", "step": "10"}
        fx_ck2 = FixedChunk(**kwargs)
        lines = ("1", "2", "3", "4", "5")
        for l in lines:
            fx_ck2.parse_data_line(l)
        ch.add_chunk(fx_ck2)
        kwargs = {"chrom": ch_name, "start": "200"}
        fx_ck3 = FixedChunk(**kwargs)
        lines = ("-1", "-2", "-3", "-4", "-5")
        for l in lines:
            fx_ck3.parse_data_line(l)
        ch.add_chunk(fx_ck3)

        self.assertTrue('chrI' in genome)
        chrI = genome['chrI']
        self.assertEqual(len(chrI._chunks), 3)
        self.assertEqual(chrI._chunks[0], fx_ck1)
        self.assertEqual(chrI._chunks[1], fx_ck2)
        self.assertEqual(chrI._chunks[2], fx_ck3)
