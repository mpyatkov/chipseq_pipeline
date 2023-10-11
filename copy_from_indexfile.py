# copy samples which specified in index_file
# python3 ./copy_from_indexfile.py <index_file.csv>
import sys
from pathlib import Path
import os

RIGHT_EXTENSION = ".fastq.gz"

def read_index_file(fname):
    with open(fname, "r") as f:
        index_file = f.readlines()
        first_line = [x.strip() for x in index_file[0].split(",")]
        prj_name, prj_path = first_line[1], first_line[2]

        reads = []
        for line in index_file[1:]:
            tmp = [x.strip() for x in line.split(",")]
            name = tmp[0].strip()
            R1 = Path(prj_path) / Path(tmp[1].strip())
            R2 = Path(prj_path) / Path(tmp[2].strip())

            reads.append((name, R1, R2))

        return prj_name, prj_path, sorted(reads, key = lambda x: x[0])

def modify_extension(fname):
    ## change extension from fq.gz etc... to fastq.gz
    ## if it required
    outname = fname
    if not str(fname).endswith(RIGHT_EXTENSION):
        extension = "".join(fname.suffixes)
        outname = str(fname).replace(extension, RIGHT_EXTENSION)
    
    # change _1,_2 -> _R1,_R2 if it is required
    outname = outname.replace("_1", "_R1").replace("_2", "_R2")

    return Path(outname)

def create_dir_and_make_links(rd): #(name, R1, R2)
    current_path = Path().absolute()
    sample_path = current_path / Path(rd[0]+'/fastq')

    ## create sample dir
    sample_path.mkdir(parents=True, exist_ok = True)


    ## create symlinks
    R1_dist = sample_path / modify_extension(rd[1].name)
    if not R1_dist.is_symlink():
        os.symlink(rd[1], R1_dist)
        
    R2_dist = sample_path / modify_extension(rd[2].name)
    if not R2_dist.is_symlink():
        os.symlink(rd[2], R2_dist)

    
if __name__ == "__main__":
    fname_path = sys.argv[1]
    prj_name, prj_path, reads = read_index_file(fname_path)

    for rd in reads:
        print(f"Processing: {rd[0]}")
        create_dir_and_make_links(rd)

    
