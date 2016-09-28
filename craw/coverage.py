

def get_coverage(sam_file, annot_entry, before=None, after=None, qual_thr=15, max_left=0, max_right=0):
    """

    :param sam_file: the samfile openend with pysam
    :type sam_file: :class:`pysam.AlignmentFile` object.
    :param annot_entry: an entry of the annotation file
    :type annot_entry: :class:`annotation.Entry` object
    :param before: The number of bases to take in account before the reference position
    :type before: int
    :param after: The number of bases to take in account before the reference position
    :type after: int
    :param qual_thr: The quality threshold
    :type qual_thr: int
    :param max_left: The highest number of base before the reference position to take in account.
    :type max_left: int
    :param max_right: The highest number of base after  the reference position to take in account.
    :type max_right: int
    :return: the coverage (all bases)
    :rtype: tuple of 2 list containing int
    """

    def on_forward(al_seg):
        """
        :param al_seg: a pysam aligned segment (the object used by pysam to represent an aligned read)
        :type al_seg: :class:`pysam.AlignedSegment`
        :return: True if read is mapped to forward strand
        :rtype: boolean
        """
        return not al_seg.is_reverse


    def on_reverse(al_seg):
        """
        :param al_seg: a pysam aligned segment (the object used by pysam to represent an aligned read)
        :type al_seg: :class:`pysam.AlignedSegment`
        :return: True if read is mapped to reverse strand.
        :rtype: boolean
        """
        return al_seg.is_reverse


    def coverage_one_strand(sam_file, chromosome, start, end, qual, strand):
        """

        :param sam_file:
        :type sam_file:
        :param chromosome: the name of the chromosome
        :type chromosome: basestring
        :param start:
        :type start: int
        :param end:
        :type end: int
        :param qual:
        :type qual: int
        :param strand:
        :type strand:
        :return: the coverage on forward then on reverse strand.
        The coverage is the sum of all kind bases mapped for each position
        :rtype: tuple of 2 list containing int
        """
        call_back = on_forward if strand == '+' else on_reverse
        coverage = sam_file.count_coverage(reference=chromosome,
                                           start=start,
                                           end=end,
                                           quality_threshold=qual,
                                           read_callback=call_back)
        coverage = [array.tolist() for array in coverage]
        window_cov = []
        for cov_A, cov_T, cov_C, cov_G in zip(*coverage):
            window_cov.append(cov_A + cov_T + cov_C + cov_G)
        return window_cov

    if before is not None:
        # pysam start base numbering at 0
        # so we need to remove 1
        start = annot_entry.ref - before - 1
        # but pysam exclude the end of interval
        # and we want to include it  -1 + 1
        stop = annot_entry.ref + after
        pad_left = []
        pad_right = []
    else:
        start = annot_entry.start -1
        stop = annot_entry.stop
        pad_left = [None] * (max_left - (annot_entry.ref - start))
        pad_right = [None] * (max_right - (stop - annot_entry.ref))
    forward_cov = coverage_one_strand(sam_file,
                                      annot_entry.chromosome,
                                      start,
                                      stop,
                                      qual_thr,
                                      '+'
                                      )
    forward_cov = pad_left + forward_cov + pad_right
    reverse_cov = coverage_one_strand(sam_file,
                                      annot_entry.chromosome,
                                      start,
                                      stop,
                                      qual_thr,
                                      '-'
                                      )
    reverse_cov = pad_left + reverse_cov + pad_right
    return forward_cov, reverse_cov

