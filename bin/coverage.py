#! /usr/bin/env python3

import argparse
import itertools


def positive_int(string):
    value = int(string)
    if value < 0:
        msg = "must be a positive integer, got: {}".format(string)
        raise argparse.ArgumentTypeError(msg)
    return value

def quality_checker(string):
    value = int(string)
    if not 0 <= value <= 42:
        msg = "must be a integer between 0 and 42, got: {}".format(string)
        raise argparse.ArgumentTypeError(msg)
    return value

parser = argparse.ArgumentParser()
parser.add_argument("--bam",
                    required=True,
                    dest='bam_file',
                    help="the path of the bam file to analyse")
parser.add_argument("--annot",
                    required=True,
                    dest='annot_file',
                    help="the path of the annotation file.")
parser.add_argument("--quality-threshold",
                    dest='qual_thr',
                    type=quality_checker,
                    default=15,
                    help="the minimal quality of read mapping to take it in account")
parser.add_argument("--suffix",
                    default="coverage",
                    help="the name of the suffix to use for the output file.")
region_grp = parser.add_argument_group(title="region of interest",
                                       description="""Parameters which define regions to compute
There is 2 way to define regions:
    - all regions have same length
    - each region have different lengths
in both case  a position of reference must be define (--ref-col default = position)
if all regions have same length:
    --window define the number of nucleotide to take in account before and after the reference position
      (the window will be centered on reference)
    --before define the number of nucleotide to take in account after the reference position.
    --after define the number of nucleotide to take in account after the reference position.
    --before and --after allow to define non centered window.
--after and --before options must be set together and are incompatible with --window option.

If all regions have different lengths. The regions must be specified in the annotation file.
    --start-col define the name of the column in annotation file which define the start position of the region to compute.
    --stop-col define the name of the column in annotation file which define the stop position of the region to compute.
""")
region_grp.add_argument("--ref-col",
                        default="position",
                        help="the name of the column for the reference position.")
region_grp.add_argument("--before",
                    type=positive_int,
                    help="the number of base to compute after the position of reference.")
region_grp.add_argument("--after",
                    type=positive_int,
                    help="the number of base to compute before the position of reference.")
region_grp.add_argument("--window",
                    type=positive_int,
                    help="the number of base to compute around the position of reference.")
region_grp.add_argument("--start-col",
                    help="the name of the column to define the start position.")
region_grp.add_argument("--stop-col",
                    help="the name of the column to define the stop position.")
args = parser.parse_args()

group_one = (args.before, args.after, args.window)
group_two = (args.start_col, args.stop_col)
if all([v is None for v in itertools.chain(group_one, group_two)]):
    raise argparse.ArgumentError("[--window or [--before, --after] or [--start-col, --stop-col] options must be specified")
elif any([v is not None for v in group_one] and any([v is not None for v in group_two])):
    raise argparse.ArgumentError("Options [--before, --after, --window] and [--start-col, --stop-col] are mutually exclusives.")
elif all([v is None for v in group_two]):
    if args.window is None:
        if any([v is None for v in (args.before, args.after)]):
            raise argparse.ArgumentError("The two options --after and --before work together. The both options must be specified in same time")
        else:
            pass
            # window is None, before and after are specify
            # => nothing to do
    else:
        # args.window is not None:
        if any([v is not None for v in (args.before, args.after)]):
            raise argparse.ArgumentError("options [--before, --after] and --window are mutually exclusives.")
        else:
            # --before, --after are None
            args.before = args.after = args.window
elif not all(group_two):
    raise argparse.ArgumentError("The two options --start-col and --stop-col work together. The both options must be specified in same time")
else:
    pass




