#!/usr/bin/env python3
"""
Functions used in BC_extract.py script



"""

from collections import OrderedDict
import os
import pandas as pd
import numpy as np
import time
from Bio import SeqIO
import re
import gzip

#Global variables
bclength = 8
probe_binding_site = 'CCTGCTAGTCCACGTCCATGTCCACC'

def get_barcodes(barcode_library: str, sheet_name: str ) -> list:
	"""
	Input:
		File path to excel file with barcodes
	Output:
		List of barcodes
	"""
	bc = pd.read_excel(barcode_library, sheet_name=sheet_name)
	bc.drop(0, inplace=True)

	return [barcode for barcode in bc['8nt Barcode Sequence'] if str(barcode) != 'nan']


def get_read_paths(sample_dir: str) -> list:
	"""
	Input:
		File path to directory containing reads
	Output:
		List of file paths to reads
	"""
	read_list = []
	for dirpath, dirnames, filenames in os.walk(sample_dir):

		for name in filenames:

			if name == ".DS_Store" or not name.endswith('.gz'):
				pass
			else:
				read_list.append(f"{dirpath}/{name}")

	return read_list


def return_barcode_counter(barcode_list: list) -> OrderedDict:
	"""
	Return an OrderedDict to store barcode counts for each read
	"""
	od = OrderedDict()
	for bc in barcode_list:
		od[bc] = 0

	return od

def get_pcr_regions(seq, probe_binding_start, probe_len):
    """
    Returns seq of pcr regions as one string
    """
    barcode_start = probe_binding_start + probe_len + 4
    region_1 = seq[probe_binding_start - 4 : probe_binding_start]
    region_2 = seq[barcode_start -4 : barcode_start]
    region_3 = seq[barcode_start + 8 : barcode_start + 11]

    return region_1+region_2+region_3

def calc_bc_count(read_path, barcodes, calc_bias=True):
    """
    Input:
    	File path to read
    	List of barcodes
    Output:
    	OrderedDict with num counts for each barcode in read
        Dictionary with counts of 'N' in barcode and instances when no barcode
            was found in expected location
        Modifies pcr bias df in place
    """
    barcode_counts = return_barcode_counter(barcodes)
    bias_counts = {base : list(np.zeros(11)) for base in ['A', 'C', 'G', 'T']}
    error_log = {key: 0 for key in ['N_count', 'key_errors']}

    with gzip.open(read_path, "rt") as handle:

    	for record in SeqIO.parse(handle, format = "fastq"):

                if probe_binding_site not in record.seq:
                    pass
                else:
                    seq = str(record.seq)
                    probe_binding_start = seq.find(probe_binding_site)
                    probe_len = len(probe_binding_site)

                    try:
                        barcode_start = probe_binding_start + probe_len + 4
                        barcode = seq[barcode_start : barcode_start + bclength]
                        if 'N' in barcode:
                        	error_log['N_count'] += 1
                        else:
                        	barcode_counts[barcode] += 1

                        pcr_regions = get_pcr_regions(
                                        seq=seq,
                                        probe_binding_start=probe_binding_start,
                                        probe_len=probe_len
                        )
                        try:
                            for i in range(11):
                                bias_counts[pcr_regions[i]][i] += 1
                        except IndexError:
                            #print(pcr_regions)
                            pass


                    except KeyError:
                    	error_log['key_errors'] += 1

    return barcode_counts, bias_counts, error_log

if __name__ == "__main__":
    seq = 'TGCTCTCATACGAACTCGTCCCATACCTGCTAGTCCACGTCCATGTCCACCATGATTNGCGTANATGTNGTTAGN'
    probe_binding_start = seq.find(probe_binding_site)
    probe_len = len(probe_binding_site)
    print(get_pcr_regions(seq, probe_binding_start, probe_len))
