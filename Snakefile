from pathlib import Path
import snakemake.utils
import os
import datetime
import re
import glob
import zipfile
import json
import numpy
import pandas as pd
import re
import csv


__author__ = "Marine AGLAVE"

#using: snakemake --profile /mnt/beegfs/pipelines/single-cell/profiles/slurm -s /mnt/beegfs/pipelines/single-cell/Snakefile --configfile /mnt/beegfs/userdata/m_aglave/pipeline/test_new_data/Params.yaml

sys.stderr.write("\n############################################################# \n")
sys.stderr.write("\n\n\t Single-cell RNA-seq pipeline \n\n")
sys.stderr.write("\n############################################################# \n\n")

### parameters ###################################################################################################################################
sys.stderr.write("\n#################### Setting Parameters ####################\n\n")

#https://stackoverflow.com/questions/72268814/importing-python-function-from-outside-of-the-current-folder
#rootpath = os.path.join(os.getcwd(), '..')
#sys.path.append(rootpath)

sys.path.append('/mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/common')

from function_utils import check_delimiter

STEPS = config['Steps']
PIPELINE_FOLDER = workflow.snakefile
PIPELINE_FOLDER = PIPELINE_FOLDER.replace("/Snakefile", "")
GLOBAL_TMP = config['Tmp'] if 'Tmp' in config else "/tmp"
if os.path.normpath(GLOBAL_TMP) != "/tmp" :
    if os.path.exists(GLOBAL_TMP) :
        sys.stderr.write("Temporary directory is set to: " + GLOBAL_TMP + "\n")
    else :
        sys.stderr.write(GLOBAL_TMP + " doesn't exist! Temporary directory is set to /tmp \n")
        GLOBAL_TMP = "/tmp"
        
if config["Sequencing_type"] == "long-reads":
    DESIGN_FILE_GE = config['design.file']

if "Alignment_countTable_GE" in STEPS:
    ### Sample/Project
    if config["Sequencing_type"] != "short-reads": sys.exit("Error in sequencing_type: please use \'short-reads\'")
    if 'Alignment_countTable_GE' in config and 'sample.name.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No sample.name.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'input.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No input.dir.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'output.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No output.dir.ge in configfile (Alignment_countTable_GE)!")
    ALIGN_SAMPLE_NAME_GE_RAW = config['Alignment_countTable_GE']['sample.name.ge']
    ALIGN_INPUT_DIR_GE_RAW = os.path.normpath(config['Alignment_countTable_GE']['input.dir.ge'])
    ALIGN_OUTPUT_DIR_GE = os.path.normpath(config['Alignment_countTable_GE']['output.dir.ge'])
    ALIGN_INPUT_DIR_GE = os.path.normpath(GLOBAL_TMP + "/fastq/")
    ### Index
    KINDEX_GE = config['Alignment_countTable_GE']['kindex.ge'] if 'Alignment_countTable_GE' in config and 'kindex.ge' in config['Alignment_countTable_GE'] else sys.exit("Error: No kindex.ge in configfile (Alignment_countTable_GE)!")
    TR2GFILE_GE = config['Alignment_countTable_GE']['tr2g.file.ge'] if 'Alignment_countTable_GE' in config and 'tr2g.file.ge' in config['Alignment_countTable_GE'] else sys.exit("Error: No tr2g.file.ge in configfile (Alignment_countTable_GE)!")
    REF_TXT_GE = config['Alignment_countTable_GE']['reference.txt'] if 'reference.txt' in config['Alignment_countTable_GE'] else "<insert_you_reference_here>"
    ### File names
    ALIGN_SAMPLE_NAME_GE = []
    ALIGN_SYMLINK_FILES_GE = []
    ALIGN_SYMLINK_FILES_NAME_GE = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_GE_RAW),1):
        #check samples names and add "_GE" if needed
        ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i] + "_GE") if (ALIGN_SAMPLE_NAME_GE_RAW[i][len(ALIGN_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i])
        #ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "*_R2_*.f*q*"))
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_[1-4]_S*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_S[0-9]*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_[1-4]_S*_R2_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_S[0-9]*_R2_*.f*q*"))
        #files with path and extention
        ALIGN_SYMLINK_FILES_GE = ALIGN_SYMLINK_FILES_GE + [ os.path.normpath(ALIGN_INPUT_DIR_GE + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_GE_RAW[i], ALIGN_SAMPLE_NAME_GE[i])) for file in ORIG_FILES]
    #files without path and extention
    ALIGN_SYMLINK_FILES_NAME_GE = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_GE]

if "Alignment_countTable_LR_GE" in STEPS:
    if config["Sequencing_type"] != "long-reads": sys.exit("Error in sequencing_type: please use \'long-reads\'")
    if 'Alignment_countTable_LR_GE' in config and 'sample.name.ge' not in config['Alignment_countTable_LR_GE']: sys.exit("Error: No sample.name.ge in configfile (Alignment_countTable_LR_GE)!")
    if 'Alignment_countTable_LR_GE' in config and 'output.dir.ge' not in config['Alignment_countTable_LR_GE']: sys.exit("Error: No output.dir.ge in configfile (Alignment_countTable_LR_GE)!")
    if 'design.file' not in config: sys.exit("Error: No design.file in configfile !")
    if 'Alignment_countTable_LR_GE' in config and 'species' not in config['Alignment_countTable_LR_GE']: sys.exit("Error: No species in configfile (Alignment_countTable_LR_GE)!")
    
    ALIGN_OUTPUT_DIR_GE = os.path.normpath(config['Alignment_countTable_LR_GE']['output.dir.ge'])
    ALIGN_SAMPLE_NAME_GE = config['Alignment_countTable_LR_GE']['sample.name.ge']
    
    if check_delimiter(DESIGN_FILE_GE) != ',':
        raise SyntaxError("csv design file is not using commas delimiters")
    else:
        design_file=pd.read_csv(DESIGN_FILE_GE,sep=',')
        
    expected_columns=['sample_id', 'path_to_fastq', 'expected_cells', 'kit_name','kit_version']
    #check element in both lists : header of the design file and list of expected column
    header_check=all(x == y for x, y in 
                     zip(design_file.columns, expected_columns))

    if header_check is not True:
        raise SyntaxError("Your design file accept only this columns: sample_id, path_to_fastq, expected_cells, kit_name,kit_version")
        
    nrow=design_file.shape[0]
    if nrow < 0:
     raise ValueError("Your design file is empty")
    if nrow > len(config['Alignment_countTable_LR_GE']['sample.name.ge']):
        raise ValueError("Your design file has more sample than expected")
    
    design_file_dict=pd.read_csv(DESIGN_FILE_GE).to_dict('records')
    FASTQ_PATH_GE={i['sample_id']:i['path_to_fastq'] for i in design_file_dict if i['sample_id'] in ALIGN_SAMPLE_NAME_GE}
    
    

if "Alignment_countTable_ADT" in STEPS:
    ### Sample/Project
    if 'Alignment_countTable_ADT' in config and 'sample.name.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No sample.name.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'input.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No input.dir.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'output.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No output.dir.adt in configfile (Alignment_countTable_ADT)!")
    ALIGN_SAMPLE_NAME_ADT_RAW = config['Alignment_countTable_ADT']['sample.name.adt']
    ALIGN_INPUT_DIR_ADT_RAW = os.path.normpath(config['Alignment_countTable_ADT']['input.dir.adt'])
    ALIGN_OUTPUT_DIR_ADT = os.path.normpath(config['Alignment_countTable_ADT']['output.dir.adt'])
    ALIGN_INPUT_DIR_ADT = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/fastq/")
    ### Index
    KINDEX_ADT = config['Alignment_countTable_ADT']['kindex.adt'] if 'Alignment_countTable_ADT' in config and 'kindex.adt' in config['Alignment_countTable_ADT'] else sys.exit("Error: No kindex.adt in configfile (Alignment_countTable_ADT)!")
    TR2GFILE_ADT = config['Alignment_countTable_ADT']['tr2g.file.adt'] if 'Alignment_countTable_ADT' in config and 'tr2g.file.adt' in config['Alignment_countTable_ADT'] else sys.exit("Error: No tr2g.file.adt in configfile (Alignment_countTable_ADT)!")
    ### File names
    ALIGN_SAMPLE_NAME_ADT = []
    ALIGN_SYMLINK_FILES_ADT = []
    ALIGN_SYMLINK_FILES_NAME_ADT = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_ADT_RAW),1):
        #check samples names and add "_ADT" if needed
        ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i] + "_ADT") if (ALIGN_SAMPLE_NAME_ADT_RAW[i][len(ALIGN_SAMPLE_NAME_ADT_RAW[i])-4:] != "_ADT") else ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i])
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R2_*.f*q*"))
        #files with path and extention
        ALIGN_SYMLINK_FILES_ADT = ALIGN_SYMLINK_FILES_ADT + [ os.path.normpath(ALIGN_INPUT_DIR_ADT + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_ADT_RAW[i], ALIGN_SAMPLE_NAME_ADT[i])) for file in ORIG_FILES]
    #files without path and extention
    ALIGN_SYMLINK_FILES_NAME_ADT = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_ADT]

if "Alignment_annotations_TCR_BCR" in STEPS:
    ### Sample/Project
    if 'Alignment_annotations_TCR_BCR' in config and 'sample.name.tcr' not in config['Alignment_annotations_TCR_BCR'] and 'sample.name.bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No sample.name.tcr or sample.name.bcr in configfile (Alignment_annotations_TCR_BCR)!")
    if 'Alignment_annotations_TCR_BCR' in config and 'input.dir.tcr' not in config['Alignment_annotations_TCR_BCR'] and 'input.dir.bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No input.dir.tcr or input.dir.bcr in configfile (Alignment_annotations_TCR_BCR)!")
    if 'Alignment_annotations_TCR_BCR' in config and 'output.dir.tcr_bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No output.dir.tcr_bcr in configfile (Alignment_annotations_TCR_BCR)!")
    ALIGN_SAMPLE_NAME_TCR_RAW = config['Alignment_annotations_TCR_BCR']['sample.name.tcr'] if 'sample.name.tcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_SAMPLE_NAME_BCR_RAW = config['Alignment_annotations_TCR_BCR']['sample.name.bcr'] if 'sample.name.bcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_INPUT_DIR_TCR_RAW = os.path.normpath(config['Alignment_annotations_TCR_BCR']['input.dir.tcr'] + "/") if 'input.dir.tcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_INPUT_DIR_BCR_RAW = os.path.normpath(config['Alignment_annotations_TCR_BCR']['input.dir.bcr'] + "/") if 'input.dir.bcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_OUTPUT_DIR_TCR_BCR = os.path.normpath(config['Alignment_annotations_TCR_BCR']['output.dir.tcr_bcr'])
    ALIGN_INPUT_DIR_TCR_BCR = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/fastq/")
    ### Index
    CRINDEX_TCR_BCR=config['Alignment_annotations_TCR_BCR']['crindex.tcr_bcr'] if ('Alignment_annotations_TCR_BCR' in config and 'crindex.tcr_bcr' in config['Alignment_annotations_TCR_BCR']) else "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"
    ### File names
    #TCR
    ALIGN_SAMPLE_NAME_TCR = []
    ALIGN_ORIG_FILES_TCR = []
    ALIGN_SYMLINK_FILES_TCR = []
    ALIGN_SYMLINK_FILES_NAME_TCR = []
    if ALIGN_SAMPLE_NAME_TCR_RAW is not None:
        for i in range(0,len(ALIGN_SAMPLE_NAME_TCR_RAW),1):
            #check samples names and add "_TCR" if needed
            ALIGN_SAMPLE_NAME_TCR.append(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_TCR") if (ALIGN_SAMPLE_NAME_TCR_RAW[i][len(ALIGN_SAMPLE_NAME_TCR_RAW[i])-4:] != "_TCR") else ALIGN_SAMPLE_NAME_TCR.append(ALIGN_SAMPLE_NAME_TCR_RAW[i])
            ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_TCR_RAW, str(ALIGN_SAMPLE_NAME_TCR_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_TCR_RAW, str(ALIGN_SAMPLE_NAME_TCR_RAW[i]) + "*_R2_*.f*q*"))
             #files with path and extention
            ALIGN_ORIG_FILES_TCR = ALIGN_ORIG_FILES_TCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_TCR_RAW[i], ALIGN_SAMPLE_NAME_TCR[i])) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_TCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for TCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_TCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_TCR]
    else:
        ALIGN_SAMPLE_NAME_TCR_RAW = []
    #BCR
    ALIGN_SAMPLE_NAME_BCR = []
    ALIGN_ORIG_FILES_BCR = []
    ALIGN_SYMLINK_FILES_BCR = []
    ALIGN_SYMLINK_FILES_NAME_BCR = []
    if ALIGN_SAMPLE_NAME_BCR_RAW is not None:
        for i in range(0,len(ALIGN_SAMPLE_NAME_BCR_RAW),1):
            #check samples names and add "_BCR" if needed
            ALIGN_SAMPLE_NAME_BCR.append(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_BCR") if (ALIGN_SAMPLE_NAME_BCR_RAW[i][len(ALIGN_SAMPLE_NAME_BCR_RAW[i])-4:] != "_BCR") else ALIGN_SAMPLE_NAME_BCR.append(ALIGN_SAMPLE_NAME_BCR_RAW[i])
            ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_BCR_RAW, str(ALIGN_SAMPLE_NAME_BCR_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_BCR_RAW, str(ALIGN_SAMPLE_NAME_BCR_RAW[i]) + "*_R2_*.f*q*"))
            #files with path and extention
            ALIGN_ORIG_FILES_BCR = ALIGN_ORIG_FILES_BCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_BCR_RAW[i], ALIGN_SAMPLE_NAME_BCR[i])) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_BCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for BCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_BCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_BCR]
    else:
        ALIGN_SAMPLE_NAME_BCR_RAW = []
    #Fusion TCR/BCR
    ALIGN_SAMPLE_NAME_TCR_BCR_RAW = ALIGN_SAMPLE_NAME_TCR_RAW  + ALIGN_SAMPLE_NAME_BCR_RAW
    ALIGN_SAMPLE_NAME_TCR_BCR = ALIGN_SAMPLE_NAME_TCR + ALIGN_SAMPLE_NAME_BCR
    ALIGN_ORIG_FILES_TCR_BCR = ALIGN_ORIG_FILES_TCR + ALIGN_ORIG_FILES_BCR
    ALIGN_SYMLINK_FILES_TCR_BCR = ALIGN_SYMLINK_FILES_TCR + ALIGN_SYMLINK_FILES_BCR
    ALIGN_SYMLINK_FILES_NAME_TCR_BCR = ALIGN_SYMLINK_FILES_NAME_TCR + ALIGN_SYMLINK_FILES_NAME_BCR

if "Alignment_countTable_GE" in STEPS or "Alignment_countTable_ADT" in STEPS:
    # 10X Technology
    if 'Alignment_countTable_GE' in config and 'sctech' in config['Alignment_countTable_GE']:
        SCTECH = config['Alignment_countTable_GE']['sctech']
    elif 'Alignment_countTable_ADT' in config and 'sctech' in config['Alignment_countTable_ADT']:
        SCTECH = config['Alignment_countTable_ADT']['sctech']
    else:
        SCTECH = '10xv3' # '10xv2' '10xv3'
    if SCTECH == '10xv3' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/3M-february-2018.txt' # '737K-august-2016.txt' '3M-february-2018.txt'
    elif SCTECH == '10xv2' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/737K-august-2016.txt'
    else :
        sys.exit("Error: sctech doesn't exist! Only '10xv2' and '10xv3' are available.\n")

if "Alignment_countTable_GE" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS:
    # Fastq-screen Index
    if 'Alignment_countTable_GE' in config and 'fastqscreen_index' in config['Alignment_countTable_GE']:
        FASTQSCREEN_INDEX = config['Alignment_countTable_GE']['fastqscreen_index']
    elif 'Alignment_annotations_TCR_BCR' in config and 'fastqscreen_index' in config['Alignment_annotations_TCR_BCR']:
        FASTQSCREEN_INDEX = config['Alignment_annotations_TCR_BCR']['fastqscreen_index']
    else :
        FASTQSCREEN_INDEX = "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf"

if "Alignment_countTable_GE" in STEPS or "Alignment_countTable_ADT" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS:
    # Cutadapt parameters
    ADAPTERSEQ='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    MINBASEQ=28
    # Name of conda environment
    CONDA_ENV_QC_ALIGN_GE_ADT = PIPELINE_FOLDER + "/envs/conda/QC_Alignment.yml"

if "Droplets_QC_GE" in STEPS:
    ### Sample/Project
    if 'Sequencing_type' in config and config["Sequencing_type"] == "short-reads":
        if 'Droplets_QC_GE' in config and 'sample.name.ge' in config['Droplets_QC_GE'] and 'input.dir.ge' in config['Droplets_QC_GE']:
            QC_SAMPLE_NAME_GE_RAW = config['Droplets_QC_GE']['sample.name.ge']
            QC_INPUT_DIR_GE = config['Droplets_QC_GE']['input.dir.ge']
            #check samples names and add "_GE" if needed
            QC_SAMPLE_NAME_GE = []
            for i in range(0,len(QC_SAMPLE_NAME_GE_RAW),1):
                QC_SAMPLE_NAME_GE.append(QC_SAMPLE_NAME_GE_RAW[i] + "_GE") if (QC_SAMPLE_NAME_GE_RAW[i][len(QC_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else QC_SAMPLE_NAME_GE.append(QC_SAMPLE_NAME_GE_RAW[i])
        elif 'sample.name.ge' in config['Alignment_countTable_GE'] and 'input.dir.ge' in config['Alignment_countTable_GE']  and "Alignment_countTable_GE" in STEPS:
            sys.stderr.write("Note: No sample.name.ge or input.dir.ge find in Droplets_QC_GE section of configfile; sample.name.ge and input.dir.ge will be determine from Alignment_countTable_GE step for Droplets_QC_GE step!\n")
            QC_SAMPLE_NAME_GE = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
            QC_INPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x), "KALLISTOBUS") for x in ALIGN_SAMPLE_NAME_GE]
        else:
            sys.exit("Error: No sample.name.ge or/and input.dir.ge in configfile!\n")
        if 'Droplets_QC_GE' in config and 'output.dir.ge' in config['Droplets_QC_GE'] :
            QC_OUTPUT_DIR_GE = config['Droplets_QC_GE']['output.dir.ge']
        elif 'output.dir.ge' in config['Alignment_countTable_GE'] :
            QC_OUTPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x)) for x in ALIGN_SAMPLE_NAME_GE]
            sys.stderr.write("Note: No output.dir.ge find in Droplets_QC_GE section of configfile; output.dir.ge will be determine from Alignment_countTable_GE step for Droplets_QC_GE step!\n")
        else :
            sys.exit("Error: No output.dir.ge find in configfile!\n")
    if 'Sequencing_type' in config and config["Sequencing_type"] == "long-reads":
        if 'Droplets_QC_GE' in config and 'sample.name.ge' in config['Droplets_QC_GE'] and 'input.dir.ge' in config['Droplets_QC_GE']:
            QC_SAMPLE_NAME_GE_RAW = config['Droplets_QC_GE']['sample.name.ge']
            QC_INPUT_DIR_GE = config['Droplets_QC_GE']['input.dir.ge']
            QC_SAMPLE_NAME_GE= [i+"_GE" if i[-3::] != "_GE" else i for i in QC_SAMPLE_NAME_GE_RAW]
        if 'Droplets_QC_GE' in config and 'sample.name.ge' in config['Droplets_QC_GE'] and 'output.dir.ge' in config['Droplets_QC_GE']:
            QC_OUTPUT_DIR_GE = config['Droplets_QC_GE']['output.dir.ge']
        elif 'sample.name.ge' in config['Alignment_countTable_LR_GE'] and 'output.dir.ge' in config['Alignment_countTable_LR_GE']:
            sys.stderr.write("Note: No sample.name.ge or input.dir.ge find in Droplets_QC_GE section of configfile; sample.name.ge and input.dir.ge will be determine from Alignment_countTable_LR_GE step for Droplets_QC_GE step!\n")
            QC_SAMPLE_NAME_GE_RAW = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
            QC_OUTPUT_DIR_GE=QC_OUTPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x)) for x in ALIGN_SAMPLE_NAME_GE]
            QC_SAMPLE_NAME_GE= [i+"_GE" if i[-3::] != "_GE" else i for i in QC_SAMPLE_NAME_GE_RAW]
            QC_INPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x), "/",str(x),"/gene_raw_feature_bc_matrix/") for x in QC_SAMPLE_NAME_GE_RAW]
        print(QC_INPUT_DIR_GE)
    else :
        sys.exit("Error: No sample.name.ge or/and input.dir.ge in configfile!\n")
        
    QC_SPECIES = config['Droplets_QC_GE']['species'] if ('Droplets_QC_GE' in config and 'species' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['species'] != None) else "NULL"
    QC_AUTHOR_NAME = config['Droplets_QC_GE']['author.name'].replace(", ", ",").replace(" ", "_") if ('Droplets_QC_GE' in config and 'author.name' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['author.name'] != None) else "NULL"
    QC_AUTHOR_MAIL = config['Droplets_QC_GE']['author.mail'].replace(", ", ",") if ('Droplets_QC_GE' in config and 'author.mail' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['author.mail'] != None) else "NULL"
    ### Analysis Parameters
    # Emptydrops
    QC_EMPTYDROPS_FDR = config['Droplets_QC_GE']['emptydrops.fdr'] if ('Droplets_QC_GE' in config and 'emptydrops.fdr' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['emptydrops.fdr'] != None) else "NULL"
    QC_DROPLETS_LIMIT = config['Droplets_QC_GE']['droplets.limit'] if ('Droplets_QC_GE' in config and 'droplets.limit' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['droplets.limit'] != None) else "NULL"
    QC_EMPTYDROPS_RETAIN = config['Droplets_QC_GE']['emptydrops.retain'] if ('Droplets_QC_GE' in config and 'emptydrops.retain' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['emptydrops.retain'] != None) else "NULL"
    # Translate ENSG into Gene Symbol
    QC_TRANSLATION_BOOL = config['Droplets_QC_GE']['translation'] if ('Droplets_QC_GE' in config and 'translation' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['translation'] != None) else "NULL"
    # QC cell
    QC_PCMITO_MIN = config['Droplets_QC_GE']['pcmito.min'] if ('Droplets_QC_GE' in config and 'pcmito.min' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcmito.min'] != None) else "NULL"
    QC_PCMITO_MAX = config['Droplets_QC_GE']['pcmito.max'] if ('Droplets_QC_GE' in config and 'pcmito.max' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcmito.max'] != None) else "NULL"
    QC_PCRIBO_MIN = config['Droplets_QC_GE']['pcribo.min'] if ('Droplets_QC_GE' in config and 'pcribo.min' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcribo.min'] != None) else "NULL"
    QC_PC_RIBO_MAX = config['Droplets_QC_GE']['pcribo.max'] if ('Droplets_QC_GE' in config and 'pcribo.max' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcribo.max'] != None) else "NULL"
    QC_MIN_FEATURES = config['Droplets_QC_GE']['min.features'] if ('Droplets_QC_GE' in config and 'min.features' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.features'] != None) else "NULL"
    QC_MIN_COUNTS = config['Droplets_QC_GE']['min.counts'] if ('Droplets_QC_GE' in config and 'min.counts' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.counts'] != None) else "NULL"
    # QC gene
    QC_MIN_CELLS = config['Droplets_QC_GE']['min.cells'] if ('Droplets_QC_GE' in config and 'min.cells' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.cells'] != None) else "NULL"
    ### Databases
    # Metadata file
    QC_METADATA_FILE = config['Droplets_QC_GE']['metadata.file'].replace(", ", ",") if ('Droplets_QC_GE' in config and 'metadata.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['metadata.file'] != None) else "NULL"
    # QC
    QC_MT_FILE = config['Droplets_QC_GE']['mt.genes.file'] if ('Droplets_QC_GE' in config and 'mt.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['mt.genes.file'] != None) else "NULL"
    QC_RB_FILE = config['Droplets_QC_GE']['crb.genes.file'] if ('Droplets_QC_GE' in config and 'crb.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['crb.genes.file'] != None) else "NULL"
    QC_ST_FILE = config['Droplets_QC_GE']['str.genes.file'] if ('Droplets_QC_GE' in config and 'str.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['str.genes.file'] != None) else "NULL"
    # Translation into gene Symbols
    QC_TRANSLATION_FILE = config['Droplets_QC_GE']['translation.file'] if ('Droplets_QC_GE' in config and 'translation.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['translation.file'] != None) else "NULL"
    ### Snakefile parameters
    #check end paths (del "/" if necessary)
    for i in range(0,len(QC_INPUT_DIR_GE),1):
        QC_INPUT_DIR_GE[i] = os.path.normpath(QC_INPUT_DIR_GE[i])
        QC_OUTPUT_DIR_GE[i] = os.path.normpath(QC_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_SAMPLE_NAME_GE_INFO = {}
    for i in range(0,len(QC_SAMPLE_NAME_GE),1):
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]] = {}
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]]['QC_INPUT_DIR'] = QC_INPUT_DIR_GE[i]
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]]['QC_OUTPUT_DIR'] = QC_OUTPUT_DIR_GE[i]

#singularity
if "Droplets_QC_GE" in STEPS or "Filtering_GE" in STEPS or "Norm_DimRed_Eval_GE" in STEPS or "Clust_Markers_Annot_GE" in STEPS or "Adding_ADT" in STEPS or "Int_Clust_Markers_Annot_GE" in STEPS or "Grp_Norm_DimRed_Eval_GE" in STEPS or "Grp_Clust_Markers_Annot_GE" in STEPS:
    SINGULARITY_ENV = PIPELINE_FOLDER + "/envs/singularity/single_cell.simg"
if "Int_Norm_DimRed_Eval_GE" in STEPS :
    INT_SINGULARITY_ENV = PIPELINE_FOLDER + "/envs/singularity/single_cell_integration.simg"
if "Alignment_annotations_TCR_BCR" in STEPS or "Adding_TCR" in STEPS or "Adding_BCR" in STEPS or "Int_Adding_TCR" in STEPS or "Int_Adding_BCR" in STEPS or "Grp_Adding_TCR" in STEPS or "Grp_Adding_BCR" in STEPS:
    SINGULARITY_ENV_TCR_BCR = PIPELINE_FOLDER + "/envs/singularity/single_cell_TCR_BCR.simg"


### rule all ###################################################################################################################################
sys.stderr.write("\n########################### Run ############################\n\n")

include: "rules/Rule_all.smk"
rule all:
    input:
        **get_targets(STEPS)
    message:
        "Single-cell RNA-seq pipeline done!"


### real rules ###################################################################################################################################
if "Alignment_countTable_GE" in STEPS:
    include: "rules/Alignment_countTable_GE.smk"

if "Alignment_countTable_ADT" in STEPS:
    include: "rules/Alignment_countTable_ADT.smk"

if "Alignment_annotations_TCR_BCR" in STEPS:
    include: "rules/Alignment_annotations_TCR_BCR.smk"

if "Droplets_QC_GE" in STEPS:
    include: "rules/Droplets_QC_GE.smk"

if "Filtering_GE" in STEPS:
    include: "rules/Filtering_GE.smk"

if "Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Norm_DimRed_Eval_GE.smk"

if "Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Clust_Markers_Annot_GE.smk"

if "Adding_ADT" in STEPS:
    include: "rules/Adding_ADT.smk"

if "Adding_TCR" in STEPS:
    include: "rules/Adding_TCR.smk"

if "Adding_BCR" in STEPS:
    include: "rules/Adding_BCR.smk"

if "Cerebro" in STEPS:
    include: "rules/Cerebro.smk"

if "Int_Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Int_Norm_DimRed_Eval_GE.smk"

if "Int_Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Int_Clust_Markers_Annot_GE.smk"

if "Int_Adding_ADT" in STEPS:
    include: "rules/Int_Adding_ADT.smk"

if "Int_Adding_TCR" in STEPS:
    include: "rules/Int_Adding_TCR.smk"

if "Int_Adding_BCR" in STEPS:
    include: "rules/Int_Adding_BCR.smk"

if "Grp_Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Grp_Norm_DimRed_Eval_GE.smk"

if "Grp_Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Grp_Clust_Markers_Annot_GE.smk"

if "Grp_Adding_ADT" in STEPS:
    include: "rules/Grp_Adding_ADT.smk"

if "Grp_Adding_TCR" in STEPS:
    include: "rules/Grp_Adding_TCR.smk"

if "Grp_Adding_BCR" in STEPS:
    include: "rules/Grp_Adding_BCR.smk"
    
if "Alignment_countTable_LR_GE":
    include: "rules/Alignment_countTable_LR_GE.smk"