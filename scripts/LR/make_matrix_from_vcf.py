import vcf
import pandas as pd
import glob
import re
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='create SNP matrix')

parser.add_argument('--input_directory', help='input directory file path',required=True,action='store')
parser.add_argument('--output_directory',help='output directory file path',required=True,action='store')
parser.add_argument('--prefix',help='prefix name',required=True,action='store')

args = parser.parse_args()

input_dir=args.input_directory+'/*.vcf'
barcodes_list=[re.split("/|\.",i)[-2] for i in glob.glob(pathname=input_dir)]

tmp_col_list=[]
for barcodes in barcodes_list:
    vcf_reader = vcf.Reader(open(args.input_directory+'/'+barcodes+'.vcf', 'r'))
    vcf_list=[{'CHROM':record.CHROM,'POS':record.POS,'REF':record.REF,'ALT':record.ALT,'samples':record.samples} for record in vcf_reader]
    if len(vcf_list) == 0:
        tmp_col=[{'variant_ID':"no_value",barcodes:0}]
    else:
        tmp_col=[{'variant_ID':i['CHROM']+":"+str(i['POS'])+":"+i['REF']+":"+str(i['ALT'][0]),barcodes:i['samples'][0]['DP4'][2]} for i in vcf_list]
    tmp_col_list.append(tmp_col)
    
tmp_df_dict_list=[]
for barcodes_dict in tmp_col_list:
    tmp_df=pd.DataFrame(np.array([i[list(barcodes_dict[0].keys())[1]] for i in barcodes_dict]),index=[i['variant_ID'] for i in barcodes_dict],columns=[list(barcodes_dict[0].keys())[1]]).T
    tmp_df_dict_list.append(tmp_df)
    
final_df=pd.concat(tmp_df_dict_list, axis=0, ignore_index=False)

final_df.fillna(0, inplace=True)

final_df=final_df.drop("no_value", axis='columns')

output_path_file=args.output_directory+"/matrix_SNP_"+args.prefix+".tsv"
final_df.to_csv(output_path_file,sep="\t",header=True,index=True)