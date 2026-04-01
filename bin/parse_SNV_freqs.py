#! /usr/bin/env python3

"""Get VAFs across VCF types"""

import gzip
import sys
import csv
import os
from parse_freqs_functions import print_array_to_csv, add_to_bin, fill_empty_bins, check_path_exists, \
    get_info_dict, print_info_to_tab, print_barebones_vcf, print_depths_to_tab, print_mutationTimeR_vcf, \
    get_cellularity_from_celluloid, get_segment_dictionary

BASE_COMPLEMENTS = {
        "AC": "TG",
        "AG": "TC",
        "AT": "TA",
        "GA": "CT",
        "GC": "CG",
        "GT": "CA"
    }

TUMOUR_COLUMN = 10
HISTOGRAM_BIN_SIZE = 2
MAX_HISTOGRAM_BINS = 100

VCF_REF_FIELD = 3
VCF_ALT_FIELD = 4
VCF_FORMAT_FIELD = 8

MUTECT_GENOTYPE_FORMAT = "GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB" 
NEW_STRELKA_GENOTYPE_FORMAT = "GT:AD:AF:DP:F1R2:F2R1:SB"
STRELKA_GENOTYPE_FORMAT = "DP:FDP:SDP:SUBDP:AU:CU:GU:TU"

GENOTYPE_AF_FIELD = 2
GENOTYPE_STRELKA_DP_FIELD = 0

CRITICAL_KEYS = ['rmsk','simplerepeat','superdup', 'ANNOVAR_EXONIC', 'ANNOVAR']



def get_strelka_field(basepair, tumour_genotype_parts):
    column = 4
    for comparator_basepair in ('A', 'C', 'G', 'T'):
        if basepair == comparator_basepair:
            this_depth = int(tumour_genotype_parts[column].split(",")[0])
        column = column + 1
    return(this_depth)

def get_strelka_frequency(vcf_fields, alt_allele, ref_allele):
    tumour_genotype_parts = vcf_fields[TUMOUR_COLUMN].split(":")
    total_depth = int(tumour_genotype_parts[GENOTYPE_STRELKA_DP_FIELD])

    alt_depth = get_strelka_field(alt_allele, tumour_genotype_parts)
    ref_depth = get_strelka_field(ref_allele, tumour_genotype_parts)

    this_frequency = alt_depth / total_depth
    these_depths = [ref_depth, alt_depth]
    return(this_frequency, these_depths)



def main(snv_file, outfile_root, celluloid_parameters, cnv_file):

    segments = get_segment_dictionary(cnv_file)

    frequency_bins = {}
    counts_at_sites = {}
    counts_at_dbs = {}

    with gzip.open(snv_file, 'rt') as file:
        for line in file:
            
            line = line.strip()
            if line.startswith("#"):
                continue
            vcf_fields = line.split("\t")
            ref_allele = vcf_fields[VCF_REF_FIELD]

            for alt_allele in vcf_fields[VCF_ALT_FIELD].split(","):

                # only incl. single nucleotide variants
                if len(ref_allele) == len(alt_allele) == 1:
                    this_change = ""
                    this_frequency = None

                    this_change = ref_allele + alt_allele
                    if this_change in BASE_COMPLEMENTS:
                        this_change = BASE_COMPLEMENTS[this_change]

                    if vcf_fields[VCF_FORMAT_FIELD] == NEW_STRELKA_GENOTYPE_FORMAT or vcf_fields[VCF_FORMAT_FIELD] ==  MUTECT_GENOTYPE_FORMAT:
                        tumour_genotype_parts = vcf_fields[TUMOUR_COLUMN].split(":")
                        this_frequency = float(tumour_genotype_parts[GENOTYPE_AF_FIELD])
                        these_depths = tumour_genotype_parts[1].split(",")

                    elif vcf_fields[VCF_FORMAT_FIELD] == STRELKA_GENOTYPE_FORMAT:
                        this_frequency, these_depths =  get_strelka_frequency(vcf_fields, alt_allele, ref_allele)

                    frequency_bins = add_to_bin(this_frequency, this_change, frequency_bins, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)

                    this_chromosome = vcf_fields[0]
                    this_position = int(vcf_fields[1])
                    mutation_id = ":".join((this_chromosome, str(this_position), alt_allele))

                    key_info_dict = get_info_dict(vcf_fields, CRITICAL_KEYS)

                    
                    counts_at_sites[mutation_id] = {
                        'id': mutation_id,
                        'chr': this_chromosome,
                        'pos': this_position,
                        'ref_allele': ref_allele,
                        'alt_allele': alt_allele,
                        'ref_counts' : these_depths[0],
                        'alt_counts' : these_depths[1],
                        'major_cn': 'NA',
                        'minor_cn': 'NA',
                        'info': key_info_dict
                    }

                    if this_chromosome in segments:
                        these_segments = segments[this_chromosome]
                        for this_segment in these_segments:
                            this_start = float(these_segments[this_segment]['start'])
                            this_end = float(these_segments[this_segment]['end'])
                            if this_start <= this_position and this_end > this_position  :
                                counts_at_sites[mutation_id]['major_cn'] = these_segments[this_segment]['major']
                                counts_at_sites[mutation_id]['minor_cn'] = these_segments[this_segment]['minor']
                else:

                    if vcf_fields[VCF_FORMAT_FIELD] == NEW_STRELKA_GENOTYPE_FORMAT or vcf_fields[VCF_FORMAT_FIELD] ==  MUTECT_GENOTYPE_FORMAT:
                        tumour_genotype_parts = vcf_fields[TUMOUR_COLUMN].split(":")
                        these_depths = tumour_genotype_parts[1].split(",")

                    this_chromosome = vcf_fields[0]
                    this_position = int(vcf_fields[1])
                    mutation_id = ":".join((this_chromosome, str(this_position), alt_allele))

                    key_info_dict = get_info_dict(vcf_fields, CRITICAL_KEYS)

                    counts_at_dbs[mutation_id] = {
                        'id': mutation_id,
                        'chr': this_chromosome,
                        'pos': this_position,
                        'ref_allele': ref_allele,
                        'alt_allele': alt_allele,
                        'ref_counts' : these_depths[0],
                        'alt_counts' : these_depths[1],
                        'major_cn': 'NA',
                        'minor_cn': 'NA',
                        'info': key_info_dict
                    }

                    if this_chromosome in segments:
                        these_segments = segments[this_chromosome]
                        for this_segment in these_segments:
                            this_start = float(these_segments[this_segment]['start'])
                            this_end = float(these_segments[this_segment]['end'])
                            if this_start <= this_position and this_end > this_position  :
                                counts_at_dbs[mutation_id]['major_cn'] = these_segments[this_segment]['major']
                                counts_at_dbs[mutation_id]['minor_cn'] = these_segments[this_segment]['minor']                    

    frequency_bins = fill_empty_bins(frequency_bins, BASE_COMPLEMENTS, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)
    print_array_to_csv(".".join((outfile_root,"snv")), frequency_bins, BASE_COMPLEMENTS, MAX_HISTOGRAM_BINS, HISTOGRAM_BIN_SIZE)
    cellularity = get_cellularity_from_celluloid(celluloid_parameters)
    print_depths_to_tab(".".join((outfile_root,"snv")), counts_at_sites, cellularity)
    print_info_to_tab(".".join((outfile_root,"snv")), counts_at_sites, CRITICAL_KEYS)
    print_mutationTimeR_vcf(".".join((outfile_root,"snv")), counts_at_sites)
    print_barebones_vcf(".".join((outfile_root,"snv")), counts_at_sites)

    print_depths_to_tab(".".join((outfile_root,"dbs")), counts_at_dbs, cellularity)
    print_info_to_tab(".".join((outfile_root,"dbs")), counts_at_dbs, CRITICAL_KEYS)
    print_mutationTimeR_vcf(".".join((outfile_root,"dbs")), counts_at_dbs)
    print_barebones_vcf(".".join((outfile_root,"dbs")), counts_at_dbs)

if __name__ == "__main__":
    snv_file = sys.argv[1]
    outfile_root = sys.argv[2]
    celluloid_parameters = sys.argv[3]
    cnv_file = sys.argv[4]
    main(snv_file, outfile_root, celluloid_parameters, cnv_file)