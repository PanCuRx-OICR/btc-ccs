#! /usr/bin/env python3

import csv
import os
import gzip

def get_segment_dictionary(cnv_file):
    segments = {}

    with gzip.open(cnv_file, 'rt') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#chrom") | line.startswith("chrom"):
                continue
            seg_fields = line.split("\t")
            this_chromosome = seg_fields[0]
            this_start = seg_fields[1]
            if this_chromosome not in segments:
                segments[this_chromosome] = {}
            if seg_fields[12] != 'NA':
                minor_major = seg_fields[12].split('.')
                segments[this_chromosome]['_'.join((seg_fields[1], seg_fields[2]))] = {
                    'start': seg_fields[1],
                    'end': seg_fields[2],
                    'minor': minor_major[0],
                    'major': minor_major[1],
                }
    
    return(segments)

def print_depths_to_tab(outfile_root, counts_at_sites, cellularity, sample_id = "sample1", germline_cn = 2):

    depth_keys = ['ref_counts', 'alt_counts' ]
    cn_keys = ['major_cn',  'minor_cn' ]

    with open(f"{outfile_root}.depths.txt", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['mutation_id', 'sample_id'] + depth_keys + ['normal_cn'] + cn_keys + ['tumour_content'])
        for outer_key in counts_at_sites.keys():
            if counts_at_sites[outer_key]['major_cn'] != 'NA' and \
                (counts_at_sites[outer_key]['ref_counts'] != 0 and counts_at_sites[outer_key]['alt_counts'] != 0 ):
                row = [outer_key] + [sample_id] + [counts_at_sites[outer_key].get(k, '') for k in depth_keys] + [germline_cn] +  [counts_at_sites[outer_key].get(k, '') for k in cn_keys] + [cellularity]
                writer.writerow(row)

def print_mutationTimeR_vcf(outfile_root, counts_at_sites):

    vcf_keys = ['chr', 'pos', 'id', 'ref_allele', 'alt_allele' ]

    with open(f"{outfile_root}.mutationTimeR.vcf", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', quotechar = "'")
        writer.writerow(['##fileformat=VCFv4.1' ])
        writer.writerow(['##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">' ])
        writer.writerow(['##FORMAT=<ID=AD,Number=2,Type=String,Description=\"Integer Allelic depths (number of reads in each observed allele)\">' ])
        writer.writerow(['##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Integer Total read depth\">' ])
        writer.writerow(['##FORMAT=<ID=FT,Number=1,Type=String,Description=\"String  Variant filters\">' ])
        writer.writerow(['##INFO=<ID=t_alt_count,Number=1,Type=String,Description=\"Integer Tumour alt count\">' ])
        writer.writerow(['##INFO=<ID=t_ref_count,Number=1,Type=String,Description=\"Integer Tumour ref count\">' ])
        writer.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR' ])
        for outer_key in counts_at_sites.keys():
            t_alt_count = "=".join(('t_alt_count',  str(counts_at_sites[outer_key]['alt_counts'])))
            t_ref_count = "=".join(('t_ref_count',  str(counts_at_sites[outer_key]['ref_counts'])))
            info_string = ";".join((t_alt_count, t_ref_count))
            #TODO: confirm that mutant allele is actually ref
            ad_string = ",".join((str(counts_at_sites[outer_key]['ref_counts']), str(counts_at_sites[outer_key]['alt_counts'])))
            site_depth = str(int(counts_at_sites[outer_key]['ref_counts']) + int(counts_at_sites[outer_key]['alt_counts']))
            genotype_string = ":".join(('0/1',ad_string,site_depth,'PASS'))
            row = [counts_at_sites[outer_key].get(k, '') for k in vcf_keys] + ['.', 'PASS'] + [info_string] + ['GT:AD:DP:FT'] + [genotype_string]
            writer.writerow(row)

def get_cellularity_from_celluloid(celluloid_params_file_path):
    celluloid_params_file_path = check_path_exists(celluloid_params_file_path)
    with open(celluloid_params_file_path, 'r') as celluloid_params_file:
        for row in csv.DictReader(celluloid_params_file, delimiter=" "):
            cellularity = "{:.2f}".format(float(row["T1"]))
    return(cellularity)



def add_to_bin(this_frequency, this_change, frequency_bins, max_histogram_bins, histogram_bin_size):
    bin = int((this_frequency * 100) / histogram_bin_size) * histogram_bin_size
    if bin == max_histogram_bins:
        bin = max_histogram_bins - histogram_bin_size
    if this_change not in frequency_bins:
        frequency_bins[this_change] = {}
    if bin not in frequency_bins[this_change]:
        frequency_bins[this_change][bin] = 0
    frequency_bins[this_change][bin] += 1
    return(frequency_bins)

def check_path_exists(path):
    if not os.path.exists(path):
        msg = "Cannot find file: {0}".format(path)
        raise RuntimeError(msg)
    else:
        return(path)

def fill_empty_bins(frequency_bins, variant_dictionary, max_histogram_bins, histogram_bin_size):
    for this_change in variant_dictionary.values():
        for this_bin in range(0, max_histogram_bins, histogram_bin_size):
            if this_change not in frequency_bins:
                frequency_bins[this_change] = {}
            if this_bin not in frequency_bins[this_change]:
                frequency_bins[this_change][this_bin] = 0
    return(frequency_bins)


def get_info_dict(vcf_fields, critical_keys):
    UNSPLITTABLES = ('rmsk','simplerepeat','superdup','SOMATIC', 'STR', 'DBSNP', 'OVERLAP')
    all_info_dict = {}
    for this_info in vcf_fields[7].split(';'):
        if this_info in UNSPLITTABLES:
            all_info_dict[this_info] = True
        else:
            named_info = this_info.split('=')
            all_info_dict[named_info[0]] = named_info[1]

    for this_info in ('rmsk','simplerepeat','superdup', 'ANNOVAR_EXONIC'):
        if this_info not in all_info_dict:
            all_info_dict[this_info] = False

    if 'ANNOVAR' in all_info_dict:
        key_info_dict = dict((k, all_info_dict[k]) for k in critical_keys)
    else:
        key_info_dict = None
    
    return(key_info_dict)
    
    
def print_array_to_csv(outfile_root, frequency_bins, base_complements, max_histogram_bins, histogram_bin_size):
    csvfile = open(f"{outfile_root}.frequency.csv", 'w', newline='')
    writer = csv.writer(csvfile)
    for this_change in base_complements.values():
        for this_bin in range(0, max_histogram_bins, histogram_bin_size):
            writer.writerow([this_change, this_bin, frequency_bins.get(this_change, {}).get(this_bin, 0)])
    csvfile.close()

def print_barebones_vcf(outfile_root, counts_at_sites, sample_id = "sample1"):

    vcf_keys = ['chr', 'pos']
    allele_keys = ['ref_allele', 'alt_allele']

    with open(f"{outfile_root}.barebones.vcf", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', quotechar = "'")
        writer.writerow(["chr","pos1","pos2","ref","alt","patient" ])
        for outer_key in counts_at_sites.keys():
            if len(counts_at_sites[outer_key]['alt_allele']) >= len(counts_at_sites[outer_key]['ref_allele']):
                allele_length = 0
            else:
                allele_length = len(counts_at_sites[outer_key]['ref_allele']) - len(counts_at_sites[outer_key]['alt_allele'])
            pos2 = int(counts_at_sites[outer_key]['pos']) + allele_length
            row = [counts_at_sites[outer_key].get(k, '') for k in vcf_keys] + [pos2] + [counts_at_sites[outer_key].get(k, '') for k in allele_keys]+ [sample_id] 
            writer.writerow(row)


def print_info_to_tab(outfile_root, counts_at_sites, critical_keys):

    with open(f"{outfile_root}.info.txt", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['mutation_id'] + critical_keys)
        for outer_key in counts_at_sites.keys():
            if counts_at_sites[outer_key]['info'] is not None:
                row = [outer_key] + [counts_at_sites[outer_key]['info'].get(k, '') for k in critical_keys]
                writer.writerow(row)

