import os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='split barcodes in 100 text file')

parser.add_argument("--barcodes_list", help="barcodes text file",required=True,action='store')
parser.add_argument("--output_dir", help="output dir path",required=True,action='store')
parser.add_argument("--sample_id", help="output dir path",required=True,action='store')
parser.add_argument("--batch_id", help="number between 0 to 99",required=True,action='store')

args = parser.parse_args()

os.makedirs(args.output_dir+"barcodes_batch",exist_ok=True)

list_of_barcodes=np.array_split(pd.read_csv(args.barcodes_list,sep="\t",header=None).iloc[:, 0].values,100)

file_name=args.output_dir+"barcodes_batch/"+args.sample_id+"_barcodes_list_"+str(args.batch_id)+".txt"
with open(file_name, 'w') as f:
    for line in list_of_barcodes[int(args.batch_id)]:
        f.write(f"{line}\n")
    f.close()

#with open(args.args.output_dir+"tmp/"+args.sample_id+"_barcodes_list_done.txt", 'w') as f:
#    pass
#    f.close()
 