import os
import pandas as pd
import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

'''writes the info file &
a file with all abricate results from the different DBs
sorted by the contig they were found on & their starting position'''
def write_all_file(all_path, info_out_path, strain, dbs, assembly, logfile):
    first_file = True
    info_ls = []
    outdir = Path(all_path).parent

    ## runs abricate for each database in the DB list & saves all information 
    for db in dbs:
        abricate_out = "{}/{}_{}.tsv".format(outdir, strain, db)
        abricate(db, assembly, abricate_out, logfile)

        tsv_df = pd.read_csv(abricate_out, sep="\t")
        info_ls.append((db, len(tsv_df)))
        
        if not tsv_df.empty:
            if first_file:
                write_info_file(info_ls, info_out_path, assembly, first_file)

                all_df = tsv_df.copy()
                first_file = False
            else:
                all_df = pd.concat([all_df,tsv_df], ignore_index=True)

        ## after saving the information the abricate output file is deleted 
        remove_abricate(abricate_out)
    
    to_drop = ["COVERAGE_MAP","GAPS", "#FILE"]
    to_rename = {"SEQUENCE":"CONTIG"}
    all_df.drop(columns=to_drop, inplace=True)
    all_df.rename(columns=to_rename, inplace=True)
    all_df.sort_values(by=["CONTIG", "START"], inplace=True)

    ## write output files
    write_info_file(info_ls, info_out_path, assembly, first_file)
    all_df.to_csv(all_path, index=False)

'''runs abricate with the specified DB'''
def abricate(db, assembly, outpath, logfile):
    os.system(f"abricate --db {db} {assembly} > {outpath} 2>> {logfile}") #funktionierts mit 2 >?

'''writes a file with information
on the assembly file & DBs used for the abricate analysis
& the number of results found per DB'''
def write_info_file(info_ls, info_out_path, assembly, first):
    if first:
        info_file = open(info_out_path, "w")
        info_file.write("file: {}\n\n".format(assembly))
        info_file.close()
    else:
        info_df = pd.DataFrame.from_records(info_ls, columns=["database", "#results"])
        info_df.to_csv(info_out_path, mode="a", sep="\t", index=False)

'''deletes the abricate output file after the information is saved'''
def remove_abricate(file):
    os.system("rm {}".format(file))


config = snakemake.config

strain = snakemake.wildcards.strain[0]
assembly = snakemake.input[0] #"results/final_assemblies/115_copy/assembly.fasta" #str(snakemake.input)
all_path = snakemake.output.all[0] #"{}/{}_all.csv".format(outdir, strain)
info_out_path = snakemake.output.info[0] #"{}/{}_info.txt".format(outdir, strain)

#dbs = ["card", "ncbi", "resfinder", "argannot", "megares", "ecoli_vf", "plasmidfinder", "ecoh", "vfdb"]
dbs = config['abricate_dbs']

write_all_file(all_path, info_out_path, strain, dbs, assembly, snakemake.log[0])