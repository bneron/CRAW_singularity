import pysam
import annotation

samfile = pysam.AlignmentFile(args.bam_file, "rb")


with open('toto.cov', 'w') as result_file:
    print('Name', *list(range(0 - args.before, args.after + 1)), sep='\t', file=result_file)
    for annot_entry in annotation.annot_iter(args.annot_file):

        coverage = samfile.count_coverage(reference=annot_entry.chromosome,
                                          start=annot_entry.position - args.before,
                                          end=annot_entry.position + args.after + 1,
                                          quality_threshold=args.qual_thr,
                                          read_callback='all')
        coverage = [array.tolist() for array in coverage]
        window_cov = []
        for cov_A, cov_T, cov_C, cov_G in zip(*coverage):
            window_cov.append(cov_A + cov_T + cov_C + cov_G)

        print(annot_entry.position, *window_cov, sep='\t', file=result_file)