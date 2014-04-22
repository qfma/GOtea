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
parser.add_argument("gid_to_uniprot", type=str,
                    help="A map from gid to uniprot")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================

def read_map(source):
    MAP = {}
    try:
        with open(source, "r") as handle:
            for line in handle:
                line = line.rstrip().split()
                if len(line) == 2:
                    MAP[line[1]] = line[0]
    except IOError:
        print("File does not exit!")

    return MAP


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
    MAP = read_map(args.gid_to_uniprot)
    for g in GOdict:
        if g in MAP:
            print MAP[g] + "\t" + ",".join(GOdict[g])
if __name__ == "__main__":
    main()
