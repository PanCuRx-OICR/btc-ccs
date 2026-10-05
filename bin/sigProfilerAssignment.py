#! /usr/bin/env python3

import csv
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')

from SigProfilerAssignment import Analyzer as Analyze

def index_of_largest(numbers):
    if not numbers:
        return None  # Return None if the list is empty
    largest_index = 0
    largest_number = 0 
    for i in range(0, len(numbers)):
        if float(numbers[i]) >= float(numbers[largest_index]):
            largest_index = i
            largest_number = numbers[i]
    return largest_index, round(float(largest_number), 2)

def main(sigs_file_path, outfile):
    first_row = True
    with open(sigs_file_path, 'r') as sigs_file, open(outfile, 'w') as outfile:
        for row in csv.reader(sigs_file, delimiter="\t"):
            if first_row:
                signature_indeces = row[4:]
                first_row = False
            if row[0] != 'Sample Names' and row[4] != '':

                chr = ''.join(('chr',row[1]))
                largest_index, largest_number = index_of_largest(row[4:])
                joined_sbs = signature_indeces[largest_index]
                outfile.write('\t'.join((row[0], chr, row[2], row[3], joined_sbs, str(largest_number))) + '\n')

if __name__ == "__main__":
    snv_file = sys.argv[1]
    outfile_root = sys.argv[2]
    sample_id = sys.argv[3]
    signature_set_path = sys.argv[4]
    context_type = sys.argv[5]

    Analyze.cosmic_fit(samples=snv_file, 
                        output=outfile_root, 
                        input_type="vcf", 
                        context_type=context_type,
                        cosmic_version=3.4,                     
                        exome=False,
                        genome_build="GRCh38", 
                        collapse_to_SBS96 = False,
                        signature_database=signature_set_path,
                        export_probabilities=True,
                        export_probabilities_per_mutation=True, 
                        make_plots=False,
                        sample_reconstruction_plots=False, 
                        verbose=True)

    in_file_name = ''.join((('Decomposed_Mutation_Probabilities_',sample_id,'.txt')))
    outfile_file_name = ''.join(((sample_id,'.most_probable_signatures.txt')))
    sigs_file = os.path.join(outfile_root,'Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/', in_file_name)
    this_outfile  = os.path.join(outfile_root,outfile_file_name )
    main(sigs_file, this_outfile)