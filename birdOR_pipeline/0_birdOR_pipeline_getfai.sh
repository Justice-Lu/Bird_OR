#!/bin/bash


: '
The purpose of this pipeline is to 
1) Get genome
2) Create blast database
3) Create fai reference file 

birdOR_pipeline_getfai.sh requires two arguments to run:
    - arg1={NCBI genome URL}
    - arg2=path/to/output/directory (Preferrably in the format of Species_db)

Example to call the script 
bash birdOR_pipeline_getfai.sh \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/639/285/GCA_017639285.1_CAU_Pekin_2.0/GCA_017639285.1_CAU_Pekin_2.0_genomic.fna.gz \
    /path/to/species/

bash birdOR_pipeline_getfai.sh \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/637/925/GCA_023637925.1_DD_ASM_B1a/GCA_023637925.1_DD_ASM_B1a_genomic.fna.gz \
    /data/jlu/Bird_OR/Misc/DdAsm
'

# To the run the script requires a argument for a NCBI genome URL 
url="$1"
genome=$(basename "$url")

# Check if there is second argument for the path for output 
if [ -z "$2" ]; then 
    echo "Error: 2nd Argument not provided. Please provide a output directory"
    exit 1
else 
    SPP_PATH="$2"
    spp=$(basename "$SPP_PATH")
    spp_db="$SPP_PATH"/"$spp"_db


    # Check if the directory doesn't exist
    if [ ! -d "$spp_db" ]; then
        mkdir -p "$spp_db"
    fi
fi 


# Change directory 
cd $SPP_PATH
# source /data/jlu/anaconda3/etc/profile.d/conda.sh
source activate bird

# Download genome url 
echo "=== Downloading from $genome ==="
wget $url


# Unzip file
echo "=== gunzip $genome ==="
gunzip "$genome"
genome="${genome%.gz}"

# Build blast database
echo "=== Build blast database $genome ==="
samtools faidx $genome
# Create fai reference file 
echo "=== Create fai reference $genome ==="
makeblastdb -in $genome -out "$spp" -dbtype 'nucl' -hash_index


# Organize output files 
mv "${SPP_PATH}"/"$spp".* $spp_db
mv "$genome" "$spp"_genome.fna
mv "$genome".fai "$spp"_genome.fna.fai

echo -e "=== Completed ===\nOutput:\n - $SPP_PATH"


# rm $genome
