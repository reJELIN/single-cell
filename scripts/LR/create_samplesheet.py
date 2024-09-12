import os
import pandas as pd
import argparse
#from pathlib import Path

parser = argparse.ArgumentParser(description='Create samplesheet for nextflow epi2me: wf-single-cell')

parser.add_argument("--sample_id", help="sample_id",required=True,action='store')
parser.add_argument("--output_path", help="output path",required=True,action='store')
parser.add_argument("--design_file_path", help="output path",required=True,action='store')

args = parser.parse_args()

design_file=pd.read_csv(args.design_file_path,sep=",")

os.makedirs(args.output_path+"/samplesheet",exist_ok=True)

output_file=args.output_path+"/samplesheet/"+args.sample_id+"_samplesheet.csv"

row = design_file.loc[design_file["sample_id"] == args.sample_id]

row[["sample_id", "kit_name","kit_version","expected_cells"]].to_csv(output_file,sep=",",index=False)