import pandas as pd
import glob
import re
import argparse

parser = argparse.ArgumentParser(description='create SNP matrix')

parser.add_argument('--input_directory', help='input directory file path',required=True,action='store')
parser.add_argument('--output_directory',help='output directory file path',required=True,action='store')
parser.add_argument('--prefix',help='prefix name',required=True,action='store')

args = parser.parse_args()

input_dir=args.input_directory+'*.tsv'
print("input_dir :"+input_dir)
barcodes_list=[i.split("/")[-1].split(".tsv")[0] for  i in glob.glob(pathname=input_dir)]

tmp_col_list=[]
print("Get Variant Info\n")
for barcode in barcodes_list:
    with open(args.input_directory+barcode+'.tsv', 'r') as tmp_tsv:
        #get read lines and get number of line
        lines = tmp_tsv.readlines()
        line_count = len(lines)
        
        #if line count superior to zero iterate through lines
        if line_count > 0:
            for line in lines:
                #split line by tab separator
                line_split=line.split("\t")
                #start iterate to 5 and get each count by nucleotides
                for i in range(5,len(line_split)):
                    line_tmp=line_split[i].split(":")
                    tmp_col={"variant_ID":line_split[0]+":"+line_split[1]+":"+line_split[2]+":"+line_tmp[0],barcode:line_tmp[1]}
                    tmp_col_list.append(tmp_col)
        elif line_count == 0:
            tmp_col={"variant_ID":"no_value",barcode:0}
            tmp_col_list.append(tmp_col)
    tmp_tsv.close()

print("Creating matrix\n")
tmp_dict_list=[]
for barcode in barcodes_list:
    tmp_barcodes_dict={}
    for i in tmp_col_list:
        #get barcodes in tmp_col_list
        barcodes_tmp=list(i.keys())[1]
        #compare barcodes to barcodes_tmp if it is the same add value to tmp_barcodes_dict and append it
        if barcodes_tmp == barcode:
           tmp_barcodes_dict[i[list(i.keys())[0]]]=i[list(i.keys())[1]]
    tmp_dict_list.append(tmp_barcodes_dict)

final_df = pd.DataFrame(tmp_dict_list,index=barcodes_list)

#final_df.fillna(0, inplace=True)

final_df=final_df.drop("no_value", axis='columns')
print("Saving")
output_path_file=args.output_directory+"/matrix_SNP_"+args.prefix+".tsv"
final_df.to_csv(output_path_file, sep='\t')