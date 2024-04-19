"""
Reorganizing my code.

Author: Yutong Zhu
Date: 2024-04-10
"""
import argparse
import os
import sys
import pybedtools
import subprocess
import pandas as pd
from multiprocessing import Pool
from helpers import write_params, get_reference, grace_period_bins, align, count_mapped_bins, append_file, combine, select_pairs

def extract_reads(ref_file, file_source, design, dir):
    """ 
    """
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
            outdir = ("{0}/"+design+"/DNA/DNA{1}").format(dir, dna_idx)
            command = "mkdir -p " + outdir
            _ = subprocess.getoutput(command)
            
            reference = grace_period_bins(ref_file, 5, file_source)

            a_forward, a_reverse = p.map(align, [[forward, reference], [reverse, reference]])

            count_mapped_bins(False, a_forward, ref_file, forward, "DNA", dna_idx, design, "f", file_source, outdir)
            count_mapped_bins(False, a_reverse, ref_file, reverse, "DNA", dna_idx, design, "r", file_source, outdir)

        for rna_idx in rnas:
            print("Starting RNA replicate" + str(rna_idx) + " ... ")
            outdir = ("{0}/"+design+"/RNA/RNA{1}").format(dir, rna_idx)
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
            a_forward, a_reverse = p.map(align, [[forward, reference], [reverse, reference]])


            count_mapped_bins(True, a_forward, ref_file, forward, "RNA", rna_idx, design, "f", file_source, outdir)
            count_mapped_bins(True, a_reverse, ref_file, reverse, "RNA", rna_idx, design, "r", file_source, outdir)
                
    ories = ["f", "r"]
    for orie in ories:
        input_ls = append_file(dir, design, orie, dnas, rnas)
        combine(input_ls, dir, design, orie)


def parse_args():
    parser = argparse.ArgumentParser(description='Generate enhancer read counts bedfile from STARR-seq counts.bed file')
    parser.add_argument('--outdir', required=True, help="Output directory")
    parser.add_argument('--enh_file', required=True, help="Path to .bed enhancer file, should be in the format of eight or four required fields. (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd).")
    parser.add_argument('--design', required=True, default="full", help="Processing the full or partial element.")

    # parser.add_argument('--juicebox', required=True, default="", help="path to juicebox executable or java command invoking juicer_tools.jar. eg: 'java -jar juicer_tools.jar'")
    # parser.add_argument('--resolution', default=5000, help="Resolution of HiC to download. In units of bp.")
    # parser.add_argument('--include_raw', action="store_true", help="Download raw matrix in addtion to KR")
    # parser.add_argument('--chromosomes', default="all", help="comma delimited list of chromosomes to download")
    # parser.add_argument('--chromosome-prefix', default="", 
    #                     help="Add prefix to chromosome names. For example, ENCODE processed HiC file needs to add chr")
    # parser.add_argument('--skip_gzip', action="store_true", help="dont gzip hic files")

    return parser.parse_args()


def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    ### Write params file
    write_params(args, os.path.join(args.outdir, "params.txt"))

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
        pass
    else: #partial EnhancerNet
        extract_reads(ref_file, file_source, args.design, args.outdir)
        select_pairs(args.outdir, args.design)
        
    pybedtools.cleanup(remove_all=True)
        
    ### Test if these are orientation-independent
    ### Other analyses


if __name__ == '__main__':
    args = parse_args()
    main(args)