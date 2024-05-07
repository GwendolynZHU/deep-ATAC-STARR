import pandas as pd
import numpy as np
import pybedtools
import os
from subprocess import check_call, call
from functools import reduce

def run_command(command, **args):
    print("Running command: " + command)
    return check_call(command, shell=True, **args)


def write_params(args, file):
    with open(file, 'w') as outfile:
        for arg in vars(args):
            outfile.write(arg + " " + str(getattr(args, arg)) + "\n")


def grace_period_bins(file, gSize, enh_typ):
    """
    Generate a sorted dataframe of binned fragments for respected bedfiles.
    Modified: Use 5bp grace period bins for all references. This means that the regions 
    with at most 5bp differences with the original regions are all counted as 
    possible bins.

    Parameter file: Dataframe, should be in the bed format of (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd) 
    or at least four required fields (chr, chromStart, chromEnd, name)
    Parameter gSize: Int, grace period length
    Parameter enh_typ: String, PINTS or EnhancerNet
    """
    assert enh_typ == "PINTS" or enh_typ == "EnhancerNet", "File input should be either PINTS-formatted or EnhancerNet-formatter"

    adjustments = [(i, j) for i in range(-gSize, gSize + 1) for j in range(-gSize, gSize + 1)]
    new_rows = []

    # Iterate over each row only once
    for _, row in file.iterrows():
        for i_sta, i_end in adjustments:
            new_row = {
                0: row[0], 
                1: row[1] + i_sta, 
                2: row[2] + i_end, 
                3: row[3], 
                4: row[4] if enh_typ == "PINTS" else np.nan, 
                5: row[5] if enh_typ == "PINTS" else np.nan, 
                6: row[6] if enh_typ == "PINTS" else np.nan, 
                7: row[7] if enh_typ == "PINTS" else np.nan
            }
            new_rows.append(new_row)
    
    # Create a DataFrame from all new rows
    new_df = pd.DataFrame(new_rows, columns=file.columns)

    print("Finished generating reference")
    return new_df.sort_values(by=[0, 1, 2])


def get_reference(path, type, idx):
    """ 
    Return the orientation-separated counts file in a dataframe for respective replicate.
    """
    if idx in [5,6,7]:
       df2 = pybedtools.BedTool(path+type+str(idx-4)+"/all/count.bed").to_dataframe(disable_auto_names=True, header=None)
    elif idx in [1,3]:
       df2 = pybedtools.BedTool(path+type+str(idx)+"/all/count.bed.gz").to_dataframe(disable_auto_names=True, header=None)
    else:
       df2 = pybedtools.BedTool(path+type+"1"+"/all/count.bed").to_dataframe(disable_auto_names=True, header=None)

    f = df2[df2[3] == "+"]
    r = df2[df2[3] == "-"]
    return f, r


### func for align reference with diff grace period bins
def align(ls):
    """
    Align the sorted STARR-seq segments with different grace period bins.
    Return a list that contains 100% sequence coverage file:
    (chr, start, end, counts)
    """
    input = pybedtools.BedTool.from_dataframe(ls[0]) # forward/reverse counts.bed.gz
    file = pybedtools.BedTool.from_dataframe(ls[1]) # binned reference
    # print(file)
    overlap = file.coverage(input, sorted=True, counts=True, f=1.0, r=True, n=48)
    overlap = overlap.to_dataframe(disable_auto_names=True, header=None)
    return overlap


def count_mapped_bins(UMI, counts, ori_ref, ref, data_type, idx, design, orientation, enh_typ, out_dir):
    # def count_mapped_bins(UMI, counts, ori_ref, ref): # for test
    """ 
    Need to go to the original orientation separated counts.bed.gz file 
    to sum up the overlapped bins read counts if necessary (RNA).
    Sum the reads from fragmented bins.

    Parameter UMI: type boolean; True for RNA, False for DNA
    Parameter counts: dataframe (aligned fragmented bins counts) from pybedtools.coverage
    Parameter ori_ref: dataframe (elements regions after design)
    Parameter ref: forward/reverse counts.bed.gz file
    Parameter data_type: string, DNA or RNA
    Parameter idx: string, 1-6 for DNA, 1-7 for RNA
    Parameter design: string, e.g. TSS_b
    Parameter orientation: string, f or r (forward or reverse)
    Parameter enh_typ: string, PINTS or EnhancerNet
    Parameter out_dir: string
    """
    assert enh_typ == "PINTS" or enh_typ == "EnhancerNet", "Please input PINTS or EnhancerNet for the input enhancer file"
    count_idx = 8 if enh_typ == "PINTS" else 4
    counts = counts.rename(columns={count_idx: "overlap"})
    ## subset only overlapped elements from pybedtools.coverage
    counts = counts[counts["overlap"]==1].reset_index(drop=True)
    ref = ref.rename(columns={4: "count"})

    if UMI: # RNA samples
        ## replace the reads since the dataset has already been deduplicated with UMI
        ## also merge the fragmented bins through element id
        # counts = counts.merge(ref, on=[0,1,2], how="inner")
        # print(counts)
        counts["count"] = ref.merge(counts, on=[0,1,2], how="inner").loc[:,["count"]]
        counts = counts.loc[:,[0,1,2,3,"count"]]
        counts = counts.groupby(3, as_index=True).agg({0: 'first', 1: 'first', 2: 'first', "count": 'sum'})


    else: # DNA samples
        ## merge the fragmented bins through element id
        ref["count"] = 1
        counts = counts.merge(ref, on=[0,1,2], how="inner").loc[:,[0,1,2,"3_x","count"]]
        counts = counts.groupby("3_x", as_index=True).agg({0: 'first', 1: 'first', 2: 'first', "count": 'sum'})

    ## replace the values of fragmented bins by original values retrieved from element id
    if not counts.empty:
        merged_df = counts.merge(ori_ref, left_index=True, right_on=3, how='left')
        
        output = merged_df.loc[:, ["0_y","1_y", "2_y", "count"]]
        # print(output)
    else:
        output = pd.DataFrame()
    
    # output.to_csv("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/partial/data/deep_ATAC_STARR/DNA/"+data_type+idx+"_"+orientation+"_"+design+".bed", sep='\t', index=False, header=False)
    output.to_csv(out_dir+"/"+data_type+str(idx)+"_"+orientation+"_"+design+".bed", sep='\t', index=False, header=False)


def append_file(out_dir, design, orie, dnas, rnas):
    """ 
    Load the files.
    """
    DNA_files = []
    RNA_files = []

    for idx in dnas:
        outdir = ("{0}/"+design+"/DNA/DNA{1}").format(out_dir, idx)
        DNA_files.append(outdir+"/DNA"+str(idx)+"_"+orie+"_"+design+".bed")
    for idx in rnas:
        outdir = ("{0}/"+design+"/RNA/RNA{1}").format(out_dir, idx)
        RNA_files.append(outdir+"/RNA"+str(idx)+"_"+orie+"_"+design+".bed")

    ls = DNA_files + RNA_files
    input_ls = [] 
    for index in range(len(ls)):
        file = pybedtools.BedTool(ls[index])
        if file == "": # deal with completely nan file
            print("File {} is empty".format(ls[index]))
            # need to append the coordinates and zeros as placeholders, or it cannot merge in the next step??
            df = pd.DataFrame({0:["chrN"],1:[1.0],2:[2.0],3:[0.0]})
            input_ls.append(df) 
        else:
            output_file = file.to_dataframe(disable_auto_names=True, header=None)
            input_ls.append(output_file)

    return input_ls


def safe_sort(input, output):
    """ 
    """
    call("sort -k1,1 -k2,2n -V " + input + " > " + output, shell=True)


### func for combining results from biological repeats + DNA/RNA
def combine(input_list, outdir, design, orie):
    """
    Combine the alignment files into a big matrix.
    """
    output_1 = input_list[0].merge(input_list[1], how="outer",left_on=[0,1,2],right_on=[0,1,2], suffixes=('_1', '_2'))
    output_2 = output_1.merge(input_list[2], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_3"})
    output_3 = output_2.merge(input_list[3], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_4"})
    output_4 = output_3.merge(input_list[4], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_5"})
    output_5 = output_4.merge(input_list[5], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_6"})
    output_6 = output_5.merge(input_list[6], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_7"})
    output_7 = output_6.merge(input_list[7], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_8"})
    output_8 = output_7.merge(input_list[8], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_9"})
    output_9 = output_8.merge(input_list[9], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_10"})
    output_10 = output_9.merge(input_list[10], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_11"})
    output_11 = output_10.merge(input_list[11], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_12"})
    output_12 = output_11.merge(input_list[12], how="outer",left_on=[0,1,2],right_on=[0,1,2]).rename(columns={3:"3_13"})
    output_12 = output_12.fillna(0)

    output_df = output_12[output_12[0] != "chrN"]
    # print(output_12)
    unsrt_path = outdir+"/"+design+"/"+design+"_"+orie+".bed"
    srt_path = outdir+"/"+design+"/srt_"+design+"_"+orie+".bed"
    output_df.to_csv(unsrt_path, sep='\t', index=False, header=False)
    safe_sort(unsrt_path, srt_path)
    cmds = ["rm " + unsrt_path]
    
    for cmd in cmds:
        os.system(cmd)


### func for subsetting both-orientations available elements
def select_pairs(outdir, design):
    """ 
    """
    forward = pybedtools.BedTool(outdir+"/"+design+"/srt_"+design+"_f.bed").to_dataframe(disable_auto_names=True, header=None)
    reverse = pybedtools.BedTool(outdir+"/"+design+"/srt_"+design+"_r.bed").to_dataframe(disable_auto_names=True, header=None)

    overlap = pd.merge(forward, reverse, how ='inner', on = [0, 1, 2])
    fd = overlap.iloc[:,0:16]
    rv = overlap.iloc[:,list(range(0, 3)) + list(range(16, 29))]
    # print(rv)
    # print("full: ", len(overlap))
    fd.to_csv(outdir+"/"+design+"/srt_"+design+"_f.bed", sep='\t', index=False, header=False)
    rv.to_csv(outdir+"/"+design+"/srt_"+design+"_r.bed", sep='\t', index=False, header=False)
    print("Pairwise elements saved.")


### func for merging forward/reverse elements
def merge_pairs(outdir, design):
    """ 
    """
    forward = pybedtools.BedTool(outdir+"/"+design+"/srt_"+design+"_f.bed").to_dataframe(disable_auto_names=True, header=None)
    reverse = pybedtools.BedTool(outdir+"/"+design+"/srt_"+design+"_r.bed").to_dataframe(disable_auto_names=True, header=None)

    ovl = pd.concat([forward,reverse],ignore_index=True).groupby([0,1,2]).sum().reset_index()
    print(len(ovl))
    print(ovl)
    
    ovl.to_csv(outdir+"/"+design+"/srt_"+design+".bed", sep='\t', index=False, header=False)
    print("Forward & reverse reads combined for {}".format(design))


### func for generating partial elements
def generate_partial_elements(input_file, outdir):
    """
    Generate bedfile for all designs.
    """
    full_enh = pybedtools.BedTool(input_file).to_dataframe()

    config = {
        ("pause_site", "b"): (-15,-15,"pause_site_b"),
        ("pause_site", "p"): (0,-15,"pause_site_p"),
        ("pause_site", "n"): (-15,0,"pause_site_n"),
        ("tss", "b"): (32,32,"TSS_b"),
        ("tss", "p"): (0,32,"TSS_p"),
        ("tss", "n"): (32,0,"TSS_n"),
    }

    for (condition, strand) in config.keys():
        start, end, name = config.get((condition, strand))

        partial = full_enh
        partial["start"] = full_enh["thickStart"] + start
        partial["end"] = full_enh["thickEnd"] - end
        partial = partial[partial["end"]-partial["start"]>=40]

        out_path = outdir + "/design_ref"
        os.makedirs(out_path, exist_ok=True)

        out_file = out_path + "/divergent_60bp_without_" + name + ".bed"
        partial.to_csv(out_file, sep="\t", header=False, index=False)
    
    print("Finished generating partial deletion references.")
