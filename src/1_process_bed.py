"""
Reorganizing my code.

Sample run: 
python3 1_process_bed.py --outdir /output_path/EnhancerNet --enh_file /path_to_file/Enhancer_K562_5p+60_boundaries_pause_site_b.bed --design full
python3 1_process_bed.py --outdir /output_path/EnhancerNet --enh_file /path_to_file/Enhancer_K562_5p_boundaries_pause_site_b.bed --design 5p
python3 1_process_bed.py --outdir /output_path/new_folder_name --enh_file /path_to_enh_file --design partial

Author: Yutong Zhu
Date: 2024-04-10
"""
import argparse
import os
import sys
import pybedtools
import subprocess
import pandas as pd
import tempfile
import shutil
from multiprocessing import Pool
from helpers import write_params, get_reference, grace_period_bins, align, count_mapped_bins, append_file, combine, select_pairs, generate_partial_elements, merge_pairs


# Create a temporary directory and set pybedtools to use that directory
temp_dir = tempfile.mkdtemp()
pybedtools.helpers.set_tempdir(temp_dir)


def extract_reads(ref_file, file_source, design, dir):
    """ 
    """
    print("Generating files for " + design + " ... ")
    DNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/deep_ATAC_STARR/processing_data_v1/out_DNA_no_UMI/"
    RNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/deep_ATAC_STARR/processing_data_v1/out_RNA_with_UMI/"
    new_RNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/deep_ATAC_STARR/processing_data_v1/out_RNA_with_UMI_2024_04/"
    
    dnas = list(range(1,7))
    rnas = list(range(1,8))

    get_reference_args = [(DNA_path, "DNA", dna_idx) for dna_idx in dnas]
    
    with Pool(10) as p:
        references = p.starmap(get_reference, get_reference_args)
        for dna_idx, (forward, reverse) in zip(dnas, references):
            print("Starting DNA replicate" + str(dna_idx) + " ... ")
            outdir = ("{0}/data/"+design+"/DNA/DNA{1}").format(dir, dna_idx)
            command = "mkdir -p " + outdir
            _ = subprocess.getoutput(command)
            
            reference = grace_period_bins(ref_file, 5, file_source)

            a_forward, a_reverse = p.map(align, [[forward, reference, file_source], [reverse, reference, file_source]])

            count_mapped_bins(False, a_forward, ref_file, forward, "DNA", dna_idx, design, "f", outdir)
            count_mapped_bins(False, a_reverse, ref_file, reverse, "DNA", dna_idx, design, "r", outdir)

        for rna_idx in rnas:
            print("Starting RNA replicate" + str(rna_idx) + " ... ")
            outdir = ("{0}/data/"+design+"/RNA/RNA{1}").format(dir, rna_idx)
            command = "mkdir -p " + outdir
            _ = subprocess.getoutput(command)

            if rna_idx in [1,3]:
                forward, reverse = get_reference(RNA_path, "RNA", rna_idx)
            elif rna_idx in [5,6,7]:
                forward, reverse = get_reference(new_RNA_path, "RNA", rna_idx)
            elif rna_idx == 2:
                forward, reverse = get_reference(RNA_path+"corrected_bam/", "RNA", rna_idx)
            else:
                forward, reverse = get_reference(RNA_path+"corrected_bam_RNA4/", "RNA", rna_idx)
            
            reference = grace_period_bins(ref_file, 5, file_source)
            a_forward, a_reverse = p.map(align, [[forward, reference, file_source], [reverse, reference, file_source]])


            count_mapped_bins(True, a_forward, ref_file, forward, "RNA", rna_idx, design, "f", outdir)
            count_mapped_bins(True, a_reverse, ref_file, reverse, "RNA", rna_idx, design, "r", outdir)
                
    ories = ["f", "r"]
    for orie in ories:
        input_ls = append_file(dir, design, orie, dnas, rnas)
        combine(input_ls, dir, design, orie)


def parse_args():
    parser = argparse.ArgumentParser(description='Generate enhancer read counts bedfile from STARR-seq counts.bed file')
    parser.add_argument("-o", '--outdir', required=True, help="Output directory")
    parser.add_argument('--enh_file', required=True, help="Path to .bed enhancer file, should be in the format of eight or four required fields. (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd).")
    parser.add_argument('--design', required=True, default="full", help="Processing the full or partial element.")

    return parser.parse_args()


def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    ### Write params file
    write_params(args, os.path.join(args.outdir, args.design+"_params.txt"))

    ### Grep enhancer reads - map full enhancers
    ref_file = pybedtools.BedTool(args.enh_file).to_dataframe(disable_auto_names=True, header=None)
    try:
        assert ref_file.shape[1] == 8 or ref_file.shape[1] == 4, "Please input a file with 4 (up to name) or 8 required fields (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd)"
        assert args.design in ["full", "5p", "3p", "partial"], "Please input design with one of (full, 5p, 3p, partial)"
    except AssertionError as e:
        print(e)
        sys.exit()
    
    file_source = "PINTS" if ref_file.shape[1] == 8 else "EnhancerNet"

    if args.design == "full":
        extract_reads(ref_file, file_source, args.design, args.outdir)
        select_pairs(args.outdir, args.design)

    ### Grep the partial enhancer reads
    elif file_source == "PINTS": #partial
        ### func to generate partial deletions - save that to outdir/design_ref
        generate_partial_elements(args.enh_file, args.outdir)
        designs = ["pause_site_b", "pause_site_n", "pause_site_p", "TSS_b", "TSS_n", "TSS_p", "INR_b", "INR_p", "INR_n"]
        
        for design in designs:
            ref_path = args.outdir+"/design_ref/divergent_60bp_without_"+design+".bed"
            ref_file = pybedtools.BedTool(ref_path).to_dataframe(disable_auto_names=True, header=None)
            extract_reads(ref_file, file_source, design, args.outdir)
            merge_pairs(args.outdir, design)

    else: #partial EnhancerNet
        extract_reads(ref_file, file_source, args.design, args.outdir)
        select_pairs(args.outdir, args.design)
        
    pybedtools.cleanup(remove_all=True)
    shutil.rmtree(temp_dir)


    ### Test if these are orientation-independent
    ### Other analyses


if __name__ == '__main__':
    args = parse_args()
    main(args)