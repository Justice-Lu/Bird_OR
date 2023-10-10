#!/bin/bash

: '
Example script 
bash birdOR_pipeline_TEST.sh -spp turtle \
    -o /data/jlu/Bird_OR/Misc \
    -url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/637/925/GCA_023637925.1_DD_ASM_B1a/GCA_023637925.1_DD_ASM_B1a_genomic.fna.gz
'


# TODO argv for directory, genome, spp, run_clustalw, 
: '
Instantiating flags and arguements first 
'

run_iqtree=false
run_get_genome=false
run_command_only=false

# Loop through the command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -spp|--species_name)
            # Check if an argument is provided after -a
            if [ -n "$2" ]; then
                spp="$2"
                shift 2  # Consume both the flag and its value
            else
                echo "Error: Argument missing for -a"
                exit 1
            fi
            ;;
        -o|--output_path)
            # Check if an argument is provided after -b
            if [ -n "$2" ]; then
                output_path="$2"
                shift 2  # Consume both the flag and its value
            else
                echo "Error: Argument missing for -b"
                exit 1
            fi
            ;;

        -url|--genome_url)
            # Check if an argument is provided after -b
            if [ -n "$2" ]; then
                genome_url="$2"
                run_get_genome=true
                shift 2  # Consume both the flag and its value
            else
                echo "Error: Argument missing for -b"
                exit 1
            fi
            ;;

        -iqtree)
            run_iqtree=true
            ;;
        
        -run_command_only)
            run_command_only=true
            ;;

        -h|--help)
            # Print -help statements. 
            : '
            TODO 
            Finish specifying the print statements required for the scripts to run. 
            '
            echo "Error: Argument missing for -b"
            exit 1
            ;;
        -*)
            # Handle unknown flags
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            # Report an error for non-flag arguments
            echo "Error: Unexpected argument '$1'"
            exit 1
            ;;
    esac
done



: '
>>> Pipeline starts here
'

# Argument processing 
SPP_PATH="$output_path/$spp/"
# Extract the base name of url as the genome_filename
# genome_filename=$(basename "$genome_url")  # Extract the filename
# genome_filename="${genome_filename%.gz}"
genome="$spp"_genome.fna


echo -e "======================== birdOR pipeline for $SPP ========================\n"

# Check if the directory exists
if [ ! -d "$SPP_PATH" ]; then
    echo -e "Creating $SPP_PATH"
    mkdir -p $SPP_PATH  # Create the directory, including parent directories if needed
else
    echo -e "$SPP_PATH Directory already exists."
fi

# DONE
if [ "$run_get_genome" = true ]; then
    echo -e "======================== Download Genome ========================\n\n"
    bash 0_birdOR_pipeline_getfai.sh \
        $genome_url \
        "$SPP_PATH"
    echo -e "======================== Complete Genome ========================\n\n"
else 
    echo -e "======================== SKIP Genome ========================\nNo -url provided. Skip downloading genome\n\n"
fi

# DONE
echo -e "======================== BEGIN Blast ========================\n"
perl 1_ab_tblastN_ORs.pl \
    "$SPP_PATH"/"$spp"_db/"$spp" \
    OR_query_mini_ORN_gharial.fasta \
    $SPP_PATH
echo -e "======================== END Blast ========================\n\n"



# DONE
echo -e "======================== BEGIN Rscript Filter_ORs.R ========================\n\n"
Rscript 2_Filter_ORs.R \
        $SPP_PATH \
        "$genome".fai 

# rm "$SPP_PATH""$genome".fai 
echo -e "======================== END Rscript Filter_ORs.R ========================\n\n"


# DONE
bedtools getfasta -s -name -bed \
    "$SPP_PATH""$spp"_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed \
    -fi "$SPP_PATH""$genome" \
    -fo "$SPP_PATH"step_2_results.olfacUniqueBlastHits.stranded.fasta


perl 3_Find_ORF-ABmod_GetORFCoords.pl \
    "$SPP_PATH"step_2_results.olfacUniqueBlastHits.stranded.fasta \
    810


echo -e "======================== BEGIN einsi ========================\n\n"
cat ./Human_OR2J3.fasta \
    "$SPP_PATH"ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta > \
    "$SPP_PATH"ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta

einsi --preservecase --thread 22 --inputorder \
    "$SPP_PATH"ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta > \
    "$SPP_PATH""$spp"_step_4_result_mafftAlignment.fasta


echo -e "============ filter pseudogene ============\n\n"
python 4_filter_pseudogene.py \
    --save_fasta True \
    --query_file_path "$SPP_PATH""$spp"_step_4_result_mafftAlignment.fasta \
    --domain_file_path ./HumanORTMD.txt \
    --output_path "$SPP_PATH"
echo -e "============ einsi Re-align  ============\n\n"
einsi --preservecase --thread 22 --inputorder \
    "$SPP_PATH""$spp"_step_4_result_mafftAlignment_filtered.fasta > \
    "$SPP_PATH""$spp"_step_4_result_mafftAlignment_removed.fasta
echo -e "\n======================== END einsi ========================\n\n"



# The following chunk reads into re-aligned fasta file and find the TM1's first amino acid location and store in bash 
output=$(python 4_filter_pseudogene.py -get_tm1_pos \
    --query_file_path "$SPP_PATH""$spp"_step_4_result_mafftAlignment_removed.fasta \
    --domain_file_path ./HumanORTMD.txt)
# Extract the last 5 characters
tm1_pos=$(echo "$output" | rev | cut -c 1-5 | rev)
# Filter only the numbers
tm1_pos=$(echo "$tm1_pos" | grep -o '[0-9]*')


echo -e "\n======================== BEGIN 5 M pick up ========================\n\n"
perl 5_Start_M_pick_up.120.pl \
    "$SPP_PATH" \
    $tm1_pos

sed -i'' 's/>//g' "$SPP_PATH"step_5_result.pickedM.CoordinatesOfBestStart.txt
sed -i'' 's/Human_OR2J3 	//g' "$SPP_PATH"step_5_result.pickedM.CoordinatesOfBestStart.txt
echo -e "\n======================== END 5 M pick up ========================\n\n"


echo -e "======================== BEGIN einsi ========================\n\n"
cat ./niimura.genbank.outgroupSeqs.fasta \
    ./zebra_query.fasta \
    "$SPP_PATH"step_5_result.pickedM.fasta > \
    "$SPP_PATH"step_5_result.pickedM.wOutgroups.wRepresentatives.fasta 
einsi --preservecase --thread 22 --inputorder \
    "$SPP_PATH"step_5_result.pickedM.wOutgroups.wRepresentatives.fasta > \
    "$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.fasta 
echo -e "======================== END einsi ========================\n\n"



echo -e "======================== BEGIN clustalw ========================\n\n"

clustalw -CLUSTERING=NJ -bootstrap=1000 -KIMURA -TOSSGAPS -BOOTLABELS=node -quiet \
    -infile="$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.fasta \
    -outfile="$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.phb 

echo -e "======================== END clustalw ========================\n\n"

echo -e "======================== BEGIN R ========================\n\n"
Rscript 6_Create_final_bedfile.R \
    "$SPP_PATH"
echo -e "======================== END R ========================\n\n"

bedtools getfasta -s -name -bed \
    "$SPP_PATH""$spp"_finalORCoordinates.all.0based.bed \
    -fi "$SPP_PATH""$genome" \
    -fo "$SPP_PATH""$spp"_step_10_getFinalGenomeCoordsAndSeqs_PerRuf.step_10_results.finalORntSeqs.all.fna

if [ "$run_iqtree" = true ]; then
    echo -e "======================== BEGIN iqtree ========================\n\n"
    # Filters out >Human sequence from the step_6_result for iqtree
    awk '/^>Human/ {p = 0; next} /^>/ {p = 1} p' \
        "$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.fasta > \
        "$SPP_PATH"2step_6_result_mafftAlignment.wOutgroups.fasta
    mv "$SPP_PATH"2step_6_result_mafftAlignment.wOutgroups.fasta \
        "$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.fasta

    iqtree -s "$SPP_PATH"step_6_result_mafftAlignment.wOutgroups.fasta \
        -st AA -pre "$SPP_PATH""$spp"_einsi.iqtree \
        -m TEST -bb 1000 -alrt 1000 -nt 22 
    echo -e "======================== END iqtree ========================\n\n"
else 
    echo -e "======================== SKIP iqtree ========================\nTo run iqtree use -iqtree flag. \n"
fi


echo -e "\n\n======================== birdOR pipeline complete SUCCESSFULLY ========================\n
             ====================================== "$spp" ======================================\n"
