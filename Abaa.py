#!/usr/bin/env python3

# option

import sys
import os
import re
import numpy
import collections
from Bio import SeqIO
from argparse import ArgumentParser

__version__ = 0.1
__author__ = "Antoine Bodein"

# regex
fasta_reg = re.compile(r"(?i)\.fa$|\.fasta$")


def read_params(args):
    parser = ArgumentParser()
    parser.add_argument("input", metavar="<input dir>", type=str, help="Input directory with joined fasta files")
    parser.add_argument("output", metavar="<output dir>", type=str, help="Output directory to write filterd fasta files")
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    return vars(parser.parse_args())


def fasta_get_length(filepath):
    # return a list of length
    list_length = []
    fasta_seq = SeqIO.parse(filepath, "fasta")
    for seq in fasta_seq:
        list_length.append(len(seq))
    return list_length


def int_list_summary(length_list):
    """
        From an integer list (length of fasta sequences or occurence)
        Return basic stats:
            - Min
            - Max
            - Median
            - Mean
            - Standard Deviation
    """
    MIN = min(length_list)
    MAX = max(length_list)
    MEDIAN = numpy.median(length_list)
    MEAN = numpy.mean(length_list)
    SD = numpy.std(length_list)
    return {"MIN": MIN, "MAX": MAX, "MEDIAN": MEDIAN, "MEAN": MEAN, "SD": SD}


def length_get_frequency(length_list):
    """
        From a list of length,
        Compute frequency of each length,
        and get cutoff
    """
    # int_list_summary(length_list)
    counter = collections.Counter(length_list).most_common()  # return a list of tuple [(value, occ)]
    length_cutoff = numpy.std(counter.values)  # CUTOFF HERE !!!
    good_length = set()

    for value, occ in counter:
        if occ >= length_cutoff:
            good_length.add(value)
    return [a for a in good_length]


def write_fasta_length(filepath_in, filepath_out, length_cutoff):
    # write fasta seq when len(seq) >= cutoff
    with open(filepath_out, "w") as f_out:
        for seq in SeqIO.parse(filepath_in, "fasta"):
            if(len(seq) in length_cutoff):
                SeqIO.write(seq, f_out, "fasta")


def check_input(input_filepath):
    """
        check if filepath is valid,
        and the directory contains any fasta files
    """
    if not os.path.exists(input_filepath):
        raise IOError("<input dir> is not a valid filepath")

    if not os.listdir(input_filepath):
        raise IOError("<input dir> is empty")

    if not [f for f in os.listdir(input_filepath) if fasta_reg.search(f)]:
        # identify .fa or .fasta file
        raise IOError("<input dir> contains any .fa or .fasta files.")


def check_output(output_filepath):
    """
        check if directory exists and is writable
    """

    if not os.path.exists(output_filepath):
        raise IOError("<output dir> is not a valid filepath")

    if not os.access(output_filepath, os.W_OK):
        raise IOError("<output dir> is not a writable directory")


if __name__ == "__main__":
    args = read_params(sys.argv)

    print args

    # verify if both input dir and output dir are valid
    check_input(input_filepath=args["input"])

    # verify output
    check_output(output_filepath=args["output"])

    # foreach fasta files of input dir get length
    length_storing = {}  # dico to store length {filename : [length]}
    good_length_storing = {}  # dico to store good_length {filename : [good_length]}

    # 1) get length
    for f in [f for f in os.listdir(args["input"]) if fasta_reg.search(f)]:
        l = fasta_get_length("/".join([args["input"], f]))
        if l:
            length_storing[f] = l

    # 2) get frequency cutoff + filter length
    for f in length_storing.keys():
        good_length_storing[f] = length_get_frequency(length_storing[f])

    # 3) write good seq
    for f in good_length_storing.keys():
        # get filepath out
        filename, ext = os.path.splitext(f)
        filepath_in = "/".join([args["input"], f])
        filepath_out = "/".join([args["output"], filename + "_filtered" + ext])
        # filepath_stat = "/".join([args["output"], filename + ".stat"])

        write_fasta_length(filepath_in, filepath_out, good_length_storing[f])

exit()
