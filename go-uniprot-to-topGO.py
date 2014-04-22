#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#       go-uniprot-to-topGO.py
#==============================================================================
import argparse
import sys
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gene_association", type=str,
                    help="A Uniprot gene association file")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def read_gene_association(source):
    GOdict = defaultdict(list)
    try:
        with open(source, "r") as handle:
            for line in handle:
                if line.startswith("!"):
                    continue
                line = line.rstrip().split()
                GOdict[line[1]].append(line[3])
        return GOdict
    except IOError:
        print("File does not exit!")
#==============================================================================


def main():
    GOdict = read_gene_association(args.gene_association)
    for g in GOdict:
        print g,"\t",",".join(GOdict[g])
if __name__ == "__main__":
    main()
