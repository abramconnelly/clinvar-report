#!/usr/bin/python
#

from __future__ import absolute_import
from twobit import TwoBitFile

from argparse import ArgumentParser
import sys




def get_reference_allele(chrom, start, hg19_path):
    twobit_file = TwoBitFile(hg19_path)
    end = start + 1
    refallele = twobit_file[chrom][start:end]
    return refallele

def get_reference_seq(chrom, start, end, hg19_path):
    twobit_file = TwoBitFile(hg19_path)
    refallele = twobit_file[chrom][start:end]
    return refallele

def main():
    # Parse options
    parser = ArgumentParser()

    parser.add_argument("-c", "--chrom", dest="chrom",
                      help="chromosome", metavar="CHROM")
    parser.add_argument("-s", "--start", dest="start",
                      help="start", metavar="START")
    parser.add_argument("-e", "--end", dest="end",
                      help="end", metavar="END")
    parser.add_argument("-r", "--ref2bit", dest="twobitref",
                      help="2bit reference genome file",
                      metavar="TWOBITREF")
    options = parser.parse_args()

    if not options.twobitref:
      print "Provide 2bit refernce file\n"
      parser.print_help()
      sys.exit(1)

    print get_reference_seq(options.chrom, int(options.start), int(options.end), options.twobitref)


if __name__ == "__main__":
  main()
