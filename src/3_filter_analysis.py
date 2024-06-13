"""
Downstream analysis.

Sample run:
python3 3_filter_analysis.py -o ../new_data --enh_file /path_to_enh_file --design TSS pause_site
Author: Yutong Zhu
Date: 2024-4-24
"""
import os
import pandas as pd
import numpy as np
import argparse
import pybedtools
from helpers import count_reads_from_bigwig


def map_partial_to_full_coordinates(outdir, design, full_ref, full_enh):
    """ 
    design: one of (pause_site, TSS, INR)
    """
    config = {
        ("pause_site", "b"): (-15,-15,"pause_site_b"),
        ("pause_site", "p"): (0,-15,"pause_site_p"),
        ("pause_site", "n"): (-15,0,"pause_site_n"),
        ("TSS", "b"): (32,32,"TSS_b"),
        ("TSS", "p"): (0,32,"TSS_p"),
        ("TSS", "n"): (32,0,"TSS_n"),
        ("INR", "b"): (1,1,"INR_b"),
        ("INR", "p"): (0,1,"INR_p"),
        ("INR", "n"): (1,0,"INR_n"),
    }

    output = pd.DataFrame()
    for (condition, strand) in config.keys():
        if condition == design:
            start, end, name = config.get((condition, strand))
            partial = pybedtools.BedTool(os.path.join(outdir, "data", design+"_"+strand, "srt_"+design+"_"+strand+".bed")).to_dataframe(disable_auto_names=True, header=None)[[0,1,2]]
            partial[[1]] = partial[[1]] - start ## to merge with TSS from full_ref
            partial[[2]] = partial[[2]] + end
            partial.rename(columns={0:"chrom", 1:"thickStart", 2:"thickEnd"}, inplace=True)

            partial["partial_index"] = name+(partial.index+1).astype(str) ## indexes saved to select from DESeq ids

            ## if _b, can merge on chr, start, end; else, e.g. _p, start:>=60bp (ori), end:design
            if strand == "b":
                indexes = pd.merge(full_ref, partial, how ='inner', on = ["chrom", "thickStart", "thickEnd"])
                # print(set(partial["partial_index"])-set(indexes["partial_index"]))
            elif strand == "p":
                partial = partial.loc[:,["chrom", "thickEnd", "partial_index"]]
                indexes = pd.merge(full_ref, partial, how ='inner', on = ["chrom", "thickEnd"])
            else:
                partial = partial.loc[:,["chrom", "thickStart", "partial_index"]]
                indexes = pd.merge(full_ref, partial, how ='inner', on = ["chrom", "thickStart"])
            
            partial = indexes.merge(full_enh, how="inner", on=["chrom", "start", "end"]) ## partial overlapped elements with names in col:3, and coordinates in 6&7, corresponding full element with id index
            
            output = pd.concat([output,partial],ignore_index=True)
    return output


def parse_args():
    parser = argparse.ArgumentParser(description='Downstream filtering of elements')
    parser.add_argument("-o", '--outdir', required=True, help="Output directory")
    parser.add_argument('--enh_file', required=True, help="Path to original .bed enhancer file, should be in the format of eight or four required fields. (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd).")
    parser.add_argument('--design', nargs="+", required=True, help="The names of different designs. e.g. TSS, pause_site")

    return parser.parse_args()


def main(args):
    """ 
    Filtering.
    Grouping.
    """
    ### Filter out proximal enhancers
    ### Filter out enhancers overlapping with other enhancers
    ### Filter out elements with GROcap-reads less than 5
    ref = pybedtools.BedTool(args.enh_file) ## info: chr, +60bp srt & end, TSS srt & end
    # tss = pybedtools.BedTool("/fs/cbsuhy01/storage/yz2676/ref/gencode/gencode.v37.annotation.1kb.TSS.sorted.bed.gz")
    junke_tss = pybedtools.BedTool("/fs/cbsuhy01/storage/jz855/Reference/hg38/proximal_v45/promoter_1kbp_protein_coding_TSS_centered_gencode_v45.bed")

    ### step 1: filter out proximal enhancers and save only those partial that have a corresponding full element
    full = pybedtools.BedTool(os.path.join(args.outdir, "data", "full", "srt_full_e.bed"))
    overlap = full.coverage(junke_tss, counts=True, f=0.9).to_dataframe(disable_auto_names=True, header=None)
    overlap["index"] = "full"+(overlap.index+1).astype(str)
    # print(overlap)
    overlap = overlap[overlap[16]==0].loc[:,[0,1,2,"index"]] ## indexes are the DESeq ids to save, info: chr, +60bp start and end
    overlap = pybedtools.BedTool.from_dataframe(overlap)
    
    ### step 2: filter out enhancers overlapping with other enhancers
    non_overlap = overlap.intersect(ref, C=True).to_dataframe(disable_auto_names=True, header=None)
    non_overlap = non_overlap[non_overlap[4] == 1].loc[:,[0,1,2,3]]
    # print(non_overlap)

    ### step 3: make sure full enhancers are actually enhancers -> browser shots manually check 
    bw1_file_path = "/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_pl.bw"
    bw2_file_path = "/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_mn.bw"
    

    ref = ref.to_dataframe()
    ### select those partial elements that have corresponding full elements
    for design in args.design:
        deseq_file = pd.read_csv(os.path.join(args.outdir, "DESeq", "DE_results_"+design+".txt"), sep="\t", index_col=0)
        ## transform all tss_p coordinates -> tss coordinates in ref -> +60 coordinates in ref -> merge those exist in overlap -> for all designs
        ## merge the output from last step with deseq-file -> add in design (tss_b), coordinates (for visualization)
        non_overlap.rename(columns={0:"chrom", 1:"start", 2:"end", 3:"full_index"},inplace="True")
        output = map_partial_to_full_coordinates(args.outdir, design, ref, non_overlap)
        # print(output)
        deseq_index = set(deseq_file.index.to_list())
        index_dict = list((set(output["partial_index"]) | set(output["full_index"])) & deseq_index)
        deseq_available = deseq_file.loc[index_dict].sort_index()

        deseq_merged_partial = deseq_available.merge(output, left_index=True, right_on="partial_index", how="left").dropna().reset_index(drop=True)
        deseq_merged_partial = deseq_merged_partial.loc[:,["log2FoldChange", "pvalue", "padj", "chrom", "start", "end", "thickStart", "thickEnd", "partial_index", "full_index"]]
        
        condition_mapping = {
            'TSS_b': (deseq_merged_partial['thickStart'] + 32, deseq_merged_partial['thickEnd'] - 32),
            'TSS_p': (deseq_merged_partial['start'], deseq_merged_partial['thickEnd'] - 32),
            'TSS_n': (deseq_merged_partial['thickStart'] + 32, deseq_merged_partial['end']),
            'pause_site_b': (deseq_merged_partial['thickStart'] - 15, deseq_merged_partial['thickEnd'] + 15),
            'pause_site_p': (deseq_merged_partial['start'], deseq_merged_partial['thickEnd'] + 15),
            'pause_site_n': (deseq_merged_partial['thickStart'] - 15, deseq_merged_partial['end']),
            'INR_b': (deseq_merged_partial['thickStart'] + 1, deseq_merged_partial['thickEnd'] - 1),
            'INR_p': (deseq_merged_partial['start'], deseq_merged_partial['thickEnd'] - 1),
            'INR_n': (deseq_merged_partial['thickStart'] + 1, deseq_merged_partial['end'])
        }

        for prefix, (start_transform, end_transform) in condition_mapping.items():
            mask = deseq_merged_partial['partial_index'].str.startswith(prefix)
            deseq_merged_partial.loc[mask, "designStart"] = start_transform[mask]
            deseq_merged_partial.loc[mask, "designEnd"] = end_transform[mask]

        deseq_merged_partial["orientation"] = deseq_merged_partial["partial_index"].str.split("_").apply(lambda x: x[-1][0])
        # print(deseq_merged_partial)
        ### check GROcap signal within the peak and define maximum TSS and minimum TSS
        for ori in ("p", "n"):
            one_side_deletion_df = deseq_merged_partial[deseq_merged_partial["orientation"]==ori].loc[:, ["chrom", "start", "end"]]
            forward = count_reads_from_bigwig(bw1_file_path, one_side_deletion_df)
            reverse = count_reads_from_bigwig(bw2_file_path, one_side_deletion_df)
            non_exch_partial = forward.merge(reverse, on=["chr", "start", "end"])
            if ori == "p":
                non_exch_partial["TSS_group"] = np.where(non_exch_partial["expression_x"]>non_exch_partial["expression_y"], "maxi", "mini")
            else:
                non_exch_partial["TSS_group"] = np.where(non_exch_partial["expression_x"]>non_exch_partial["expression_y"], "mini", "maxi")
                
            deseq_merged_partial.loc[deseq_merged_partial["orientation"] == ori, "orientation"] = non_exch_partial["TSS_group"].values

        deseq_merged_full = deseq_available.merge(output, left_index=True, right_on="full_index", how="left").dropna().drop_duplicates(subset=["log2FoldChange","full_index"]).reset_index(drop=True)
        deseq_merged_full["partial_index"] = "nan"
        deseq_merged_full = deseq_merged_full.loc[:,["log2FoldChange", "pvalue", "padj", "chrom", "start", "end", "thickStart", "thickEnd", "partial_index", "full_index"]]
        deseq_merged_full["designStart"] = deseq_merged_full["start"]
        deseq_merged_full["designEnd"] = deseq_merged_full["end"]
        deseq_merged_full["orientation"] = "nan"
        # print(deseq_merged_full)
        
        deseq_output = pd.concat([deseq_merged_full, deseq_merged_partial])
        # print(deseq_output)
        deseq_output.to_csv(args.outdir+"/DESeq/DE_results_"+design+"_annotated.txt", sep="\t", header=True, index=None)
    print("Finished filtering the pairwise elements.")
    pybedtools.cleanup(remove_all=True)


if __name__ == '__main__':
    args = parse_args()
    main(args)