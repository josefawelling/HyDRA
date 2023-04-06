import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

'''returns a list of all indices for ARGs in the provided region
which have a lower value for the provided filter (%identity or %coverage)
than the maximum for this reason'''
def get_non_max_index_list(same_region, filter_reason):
    max_val = same_region[filter_reason].max()
    non_max_ind_lst = same_region.loc[same_region[filter_reason]!= max_val].index.tolist()
    return(non_max_ind_lst)

'''removes the genes that are marked to be filtered from the same region df
and saves the genes filtered out with the reason for filtering'''
def filter_by(same_region, filter_reason, filtout_ind_lst, filtout_genes):
    filt_out = same_region.loc[filtout_ind_lst]
    filt_out["filtered_by"] = filter_reason
    filtout_genes.update(filt_out.to_dict('index'))
    return(same_region.drop(filtout_ind_lst))

config = snakemake.config
strain = snakemake.wildcards.strain[0]
all_path = snakemake.input[0]
args_path = snakemake.params.args
filtout_args_path = snakemake.output.filtout
dbs = config['abricate_dbs']

'''
strain = "115"
path = "renameForTest_results/analysis/abricate/{}".format(strain)
all_path = "{}/{}_all.csv".format(path, strain)
info_out_path = "{}/{}_info.txt".format(path, strain)
args_out = "{0}/{1}_ARGs.csv".format(path, strain)
filtout_args_out = "{0}/{1}_filtoutARGs.csv".format(path, strain)

dbs = ["card", "ncbi", "resfinder", "argannot", "megares", "ecoli_vf", "plasmidfinder", "ecoh", "vfdb"]
'''

all_df = pd.read_csv(all_path)
all_df.sort_values(by=["CONTIG", "START"], inplace=True)
all_args = {}
filtout_args = {}

for contig in all_df["CONTIG"].unique():
    
    df_per_contig = all_df.loc[all_df["CONTIG"]==contig]
    #contig_df = df_per_contig.drop(columns="CONTIG")
    contig_df = df_per_contig.sort_values(by=["START"])

    ## filter the ARGs until there is only one entry per position
    for start in contig_df["START"].unique():
        same_start = contig_df.loc[contig_df["START"]== start]

        ## test if there is a entry with this start in the sorted df anymore
        ## since entrys for already filtered regions are removed from the df
        if len(same_start) >= 1:
            is_filtered = "no"
            ## looking for other ARGs that start within the region of the first = same region genes
            end = same_start.iloc[0]["END"]
            inner = contig_df.loc[contig_df["START"] < end]
            same_region = pd.concat([same_start, inner]).drop_duplicates()
            contig_df.drop(index=same_region.index, inplace=True)

            ## if there are more than 1 ARGs in the same region they are filtered until only 1 is left
            ##  1. by max coverage
            if len(same_region) > 1:
                is_filtered = "yes"
                non_cov_max_ind_lst = get_non_max_index_list(same_region, '%COVERAGE')
                if len(non_cov_max_ind_lst) >= 1:
                    same_region = filter_by(same_region, '%COVERAGE', non_cov_max_ind_lst, filtout_args)

            ##  2. by max identity
            if len(same_region) > 1:
                non_ident_max_ind_lst = get_non_max_index_list(same_region, '%IDENTITY')
                if len(non_ident_max_ind_lst) >= 1:
                    same_region = filter_by(same_region, '%IDENTITY', non_ident_max_ind_lst, filtout_args)
            
            ##  3. by databases    
            if len(same_region) > 1:
                ind_lst = same_region.index.tolist()

                found_best = False
                for db in dbs:
                    if found_best:
                        break
                    for ind in ind_lst:
                        if same_region['DATABASE'][ind] == db:
                            found_best = True
                            ind_lst.remove(ind)
                            break  

                same_region = filter_by(same_region, 'DATABASE', ind_lst, filtout_args)
            
            ## save all found ARGs in one dictionary with the information if they were filtered
            same_region["filtered"] = is_filtered
            all_args.update(same_region.to_dict('index'))
    
args_df = pd.DataFrame.from_dict(all_args, orient='index')
args_df.to_csv(args_path, index=False)

filtout_args_df = pd.DataFrame.from_dict(filtout_args, orient='index')
filtout_args_df.to_csv(filtout_args_path, index=False)
    