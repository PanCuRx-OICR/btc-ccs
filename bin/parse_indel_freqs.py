#! /usr/bin/env python3

import gzip
import sys
import csv
import os
import re
from parse_freqs_functions import print_array_to_csv, add_to_bin, fill_empty_bins, get_info_dict, \
	print_info_to_tab, print_barebones_vcf, print_depths_to_tab, print_mutationTimeR_vcf, \
	get_cellularity_from_celluloid, get_segment_dictionary

# set parameters
HISTOGRAM_BIN_SIZE = 2
MAX_HISTOGRAM_BINS = 100

VCF_REF_FIELD = 3
VCF_ALT_FIELD = 4
VCF_FORMAT_FIELD = 8
TUMOUR_FIELD = 10
GENOTYPE_AF_SUBFIELD = 2

INDEL_DICT = {
	"del_1": "del_1",
	"del_4": "del_4",
	"ins_1": "ins_1",
	"ins_4": "ins_4",
    }

NEW_STRELKA_GENOTYPE_FORMAT = "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50"
STRELKA_GENOTYPE_FORMAT = "GT:AD:AF:DP:F1R2:F2R1:SB"
ALT_GENOTYPE_FORMAT = "GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB"
SVABA_GENOTYPE_FORMAT = "GT:AD:DP:GQ:PL:SR:CR:LR:LO"


CRITICAL_KEYS = ['rmsk','simplerepeat','superdup', 'ANNOVAR_EXONIC', 'ANNOVAR']


def main(indel_file, outfile_root, celluloid_parameters, cnv_file):

	segments = get_segment_dictionary(cnv_file)

	frequency_bins = {}
	counts_at_sites = {}
	these_depths = [0, 0]

	print('filtering indel file')
	with gzip.open(indel_file, 'rt') as file:
		for line in file:
			
			line = line.strip()
			if line.startswith("#"):
				continue
			else:
				vcf_fields = line.split("\t")
				ref_allele = vcf_fields[VCF_REF_FIELD]

				for alt_allele in vcf_fields[VCF_ALT_FIELD].split(","):

					this_change = ""
					this_frequency = None
						
					# excl. in/dels of same length
					if len(ref_allele) == len(alt_allele):
						continue
					elif len(ref_allele) > len(alt_allele):
						if (len(ref_allele) - len(alt_allele) > 3):
							this_change = 'del_4'
						else:
							this_change = 'del_1'
					else:
						if (len(alt_allele) - len(ref_allele) > 3):
							this_change = 'ins_4'
						else:
							this_change = 'ins_1'
					
					if vcf_fields[VCF_FORMAT_FIELD] == NEW_STRELKA_GENOTYPE_FORMAT:  # new strelka indel format
						match = re.match(r'^(.*?):.*?:.*?,.*?:(.*?),.*', vcf_fields[TUMOUR_FIELD])
						if match:
							this_frequency =  int(match.group(2)) / int(match.group(1))
						else:
							raise ValueError(f"Assumed Strelka indel format, couldn't parse frequencies from {vcf_fields[TUMOUR_FIELD]}")
						
						these_depths[1] = int(match.group(2))
						these_depths[0] = int(match.group(1)) - int(match.group(2))

					elif vcf_fields[VCF_FORMAT_FIELD] in [STRELKA_GENOTYPE_FORMAT, ALT_GENOTYPE_FORMAT]:
						match = re.match(r'^(.*?):(.*?):(.*?):.*$', vcf_fields[TUMOUR_FIELD])
						if match:
							this_frequency = float(match.group(3))
					elif vcf_fields[VCF_FORMAT_FIELD] == SVABA_GENOTYPE_FORMAT:
						svaba_tumour_fields = vcf_fields[TUMOUR_FIELD].split(":")
						this_frequency = int(svaba_tumour_fields[1]) / int(svaba_tumour_fields[2])
						
						these_depths[1] = int(svaba_tumour_fields[1])
						these_depths[0] = int(svaba_tumour_fields[2]) - int(svaba_tumour_fields[1])
					else:
						print("Unrecognized Genotype Field {0}, as follows: {1}... skipping".format(vcf_fields[VCF_FORMAT_FIELD], vcf_fields[TUMOUR_FIELD]))
						continue
					
					mutation_id = ":".join((vcf_fields[0], vcf_fields[1], alt_allele))
					if these_depths[0] > 0:
						frequency_bins = add_to_bin(this_frequency, this_change, frequency_bins, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)
						
						key_info_dict = get_info_dict(vcf_fields, CRITICAL_KEYS)

						counts_at_sites[mutation_id] = {
							'id': mutation_id,
							'chr': vcf_fields[0],
							'pos': vcf_fields[1],
							'ref_allele': ref_allele,
							'alt_allele': alt_allele,
							'ref_counts' : these_depths[0],
							'alt_counts' : these_depths[1],
							'major_cn': 'NA',
							'minor_cn': 'NA',
							'info': key_info_dict
						}

						this_chromosome = vcf_fields[0]
						this_position = int(vcf_fields[1])

						if this_chromosome in segments:
							these_segments = segments[this_chromosome]
							for this_segment in these_segments:
								this_start = float(these_segments[this_segment]['start'])
								this_end = float(these_segments[this_segment]['end'])
								if this_start <= this_position and this_end > this_position  :
									counts_at_sites[mutation_id]['major_cn'] = these_segments[this_segment]['major']
									counts_at_sites[mutation_id]['minor_cn'] = these_segments[this_segment]['minor']
					else:
						print("Negative reference reads at {0}, as follows: {1}... skipping".format(mutation_id, vcf_fields[TUMOUR_FIELD]))


	print('printing plotting files')
	frequency_bins = fill_empty_bins(frequency_bins, INDEL_DICT, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)
	print_array_to_csv(".".join((outfile_root,"indel")), frequency_bins, INDEL_DICT, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)

	print('printing summary files')
	print_info_to_tab(".".join((outfile_root,"indel")), counts_at_sites, CRITICAL_KEYS)
	print_barebones_vcf(".".join((outfile_root,"indel")), counts_at_sites)
	cellularity = get_cellularity_from_celluloid(celluloid_parameters)
	print_depths_to_tab(".".join((outfile_root,"indel")), counts_at_sites, cellularity)
	print_mutationTimeR_vcf(".".join((outfile_root,"indel")), counts_at_sites)

if __name__ == "__main__":
	indel_file = sys.argv[1]
	outfile_root = sys.argv[2]
	celluloid_parameters = sys.argv[3]
	cnv_file = sys.argv[4]
	main(indel_file, outfile_root, celluloid_parameters, cnv_file)