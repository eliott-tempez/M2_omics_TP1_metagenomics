#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    seq = ""
    is_seq = True
    with gzip.open(amplicon_file, "r") as f:
        for line in f:
            # Restart sequence after passing the '>' line
            if not is_seq:
                seq = ""
            is_seq = True
            # If '>' line
            if line.startswith(b">"):
                is_seq = False
                if len(seq) >= minseqlen:
                    yield seq
            # If seqence line, add it to the sequence
            if is_seq:
                seq += line.decode().strip()
    # Add the last sequence (not followed by '>')
    if len(seq) >= minseqlen:
        yield seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    seq_dict = {}
    # Generate dictionnary with sequence count
    for seq in read_fasta(amplicon_file, minseqlen):
        if seq not in seq_dict:
            seq_dict[seq] = 1
        else:
            seq_dict[seq] += 1
    # Sort dict by count
    seq_list = sorted(seq_dict.items(), key=lambda x: x[1], reverse=True)
    # Yield values
    for item in seq_list:
        if item[1] < mincount:
            break
        yield list(item)


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    match = 0
    seq_1 = alignment_list[0]
    seq_2 = alignment_list[1]
    len_seq = len(seq_1)
    for i_nucl in range(len_seq):
        if seq_1[i_nucl] == seq_2[i_nucl] and seq_1[i_nucl] != "-":
            match += 1
    return ((match / len_seq) * 100)


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    seq_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    for i in range(len(seq_list)):
        for j in range(i, len(seq_list)):
            most_occ_seq = seq_list[i][0]
            most_occ = seq_list[i][1]
            less_occ_seq = seq_list[j][0]
            less_occ = seq_list[j][1]
            alignment = nw.global_align(most_occ_seq, less_occ_seq, gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    pass


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    output_file = args.output_file
    min_seq_len = args.minseqlen
    min_count = args.mincount

    # Cluster
    abundance_greedy_clustering(amplicon_file, min_seq_len, min_count, 0, 0)
    



if __name__ == '__main__':
    main()
