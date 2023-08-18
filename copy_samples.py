#!/bin/python3

# creates directory for each individual sample and makes links to fastq.gz files inside
## USAGE:
# python3 ./copy_samples.py <path/to/directory/with/fastq/files>

from pathlib import Path
import sys
import os

RIGHT_EXTENSION = ".fastq.gz"

## just for debug
# prelim_path = '/net/waxman-server/mnt/data/volume2/Waxman_Illumina_HiSeq_Raw_Data_02/G197/usftp21.novogene.com/01.RawData'
# prelim_path = '/net/waxman-server/mnt/data/volume2/Waxman_Illumina_HiSeq_Raw_Data_02/G196/usftp21.novogene.com/01.RawData'

## calculate hamming distance of paths, if the distance equal to 1 that means
## we found forward and reverse reads
def hamming_distance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

## calculate position of disimilarities
## I use this function to get the common part of the name and drop
## tail with extension and read information (need for lane dialog with user)
def hamming_indexes(str1,str2):
    ix = 0
    indexes = []
    for c1,c2 in zip(str1, str2):
        if c1 != c2:
            indexes.append(ix)
        ix=ix+1
    return(indexes)

## transforms initial dict of lists to dict of dicts
## each sample contains dict of lanes where each lane
## is ordered pair of reads
def transform_samples_dict(key):
    transformed_dict = {}
    for str1 in samples[key]:
        for str2 in samples[key]:
            str1 = str(str1)
            str2 = str(str2)
            
            if hamming_distance(str1, str2) == 1:
                hamming_ix = hamming_indexes(str1,str2)[0]

                ## chop _1.gz -> convert to path -> extract name -> use as a key
                hamming_prefix = Path(str1[:hamming_ix-1]).name
                hamming_prefix = str(hamming_prefix)
                
                if hamming_prefix not in transformed_dict:
                    transformed_dict[hamming_prefix] = sorted([str1,str2])
    
    return transformed_dict

def add_before_position(s, ix, val):
    return s[0:ix]+val+s[ix:] 

def create_symlink(fq_files, dst):

    files = fq_files

    
    ## set flag if it is require to change indexes _1,_2 -> R1,R2
    change_read_name = False
    f1,f2 = str(Path(files[0]).name),str(Path(files[1]).name)
    ix = hamming_indexes(f1,f2)[-1]
    
    if f1[ix-1] != "R":
        change_read_name = True
        
    
    for fastq_file in files:
        file_src = Path(fastq_file)
        fname = file_src.name

        ## change extension from fq.gz etc... to fastq.gz
        ## if it required
        if not str(file_src).endswith(RIGHT_EXTENSION):
            extension = "".join(file_src.suffixes)
            fname = fname.replace(extension, RIGHT_EXTENSION)

        ## change _1,_2 -> _R1,_R2
        ## if it is required
        if change_read_name:
            fname = add_before_position(fname,ix,"R")
           
        file_dst = Path(dst+'/'+fname)
        
        # if not file_dst.is_file():
        if not file_dst.is_symlink():
            os.symlink(fastq_file, file_dst)

def create_dir_and_make_links(samples, sample_id, lane_index):
    current_path = Path().absolute()
    # print(current_path)
    sample_path = current_path / Path(sample_id+'/fastq')
    # print(sample_path)

    ## create directory
    sample_path.mkdir(parents=True, exist_ok = True)

    ## create symlink
    lane_name = list(samples[sample_id].keys())[lane_index]
    fastq_files = samples[sample_id][lane_name]
    create_symlink(fastq_files, str(sample_path))


## are we going to copy the following samples dialog
def user_dialog_samples(sample_id):
    res = False
    while True:
        
        print(f"\nWould you like to copy/make link for {sample_id} (YyNn, default: N): ", end='')
        answer = input()
        if answer not in "Yynn":
            print("Wrong answer. Try it again")
            continue

        if answer in "Yy":
            res = True
            break
        break
    return res

        
## asks user which lane to use if script found multiple inside one
## sample directory (ex. L1, L4, etc)
## return index of selected lane
def user_dialog_lanes(d):

    selected_index = -1
    while True:
        print("-"*20)
        for index,lane in enumerate(d.keys()):
            print(f"{index}) {lane}")
        print(f"\nwhich lane would you like to copy {list(range(len(d)))}: ", end='')
        try:
            selected_index = int(input())
        except ValueError:
            print("Error: index should be integer\n")
            continue

        if selected_index not in range(len(d)):
            print(f"Error: index {selected_index} not in a range {list(range(len(d)))}\n")
        else:
            break
        
    return int(selected_index)

if __name__ == "__main__":

    if sys.argv[1] == "":
        print("Please provide proper path to directory with fastq files")
        exit(1)

    if not Path(sys.argv[1]).exists():
        print("Please provide proper path to directory with fastq files")
        exit(2)
        
    prelim_path = Path(sys.argv[1])
    print(prelim_path)
    samples = {}
    
    for path in Path(prelim_path).rglob('*.gz'):
        if path.parent.name not in samples:
            samples[path.parent.name]=[]
            samples[path.parent.name].append(path)
        else:
            samples[path.parent.name].append(path)

    for sample_id in samples.keys():
        samples[sample_id] = transform_samples_dict(sample_id)


    for ix,sample_id in enumerate(samples.keys()):

        ## if it not required to copy the following samples
        ## just jump to another one
        if not user_dialog_samples(sample_id):
            continue
        
        print(f"Processing sample #{ix+1} ({sample_id})")
        
        if len(samples[sample_id].keys()) == 1:
            print(sample_id)
            create_dir_and_make_links(samples, sample_id, 0)
        else:
            ## get index for lane first
            ix = user_dialog_lanes(samples[sample_id])
            create_dir_and_make_links(samples, sample_id, ix)

        print("#"*40)
        print(f"Directory for {sample_id} created successfully.")
        print("#"*40+"\n")
