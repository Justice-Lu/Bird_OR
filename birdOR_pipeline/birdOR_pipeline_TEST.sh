#!/bin/bash




perl Step_1_ab_tblastN_ORs.pl


# head -{first half of file} {.sum} > HEAD.sum
# tail -{second half of file} {.sum} > TAIL.sum

# R

# bedtools intersect -v -a {HEAD.bed} -b {TAIL.bed} 

# bedtools getfasta -s -name -bed {species}_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed -fi {NCBI genome} -fo step_2_results.olfacUniqueBlastHits.stranded.fasta

# perl step_3_a_Find_ORF-ABmod.pl
# perl step_3_b_Find_ORF-ABmod_withPrintStatements.pl
# perl step_3_d_GetORFCoords.pl

# cat ./Human_OR2J3.fasta ./ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta > ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta
# einsi --preservecase --thread 22 --inputorder ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta > {species}_step_4_result_mafftAlignment.fasta

# python filter_pseudogene.py \
# --save_fasta True \
# --query_file_path ./pathtofile/MyQuery.fasta \
# --domain_file_path ./pathtofile/DomainReference.txt \

# perl step_5_a_Start_M_pick_up.120.pl
# perl step_5_d_WriteOutCoords.120.pl
# sed -i'' 's/>//g' step_5_result.pickedM.CoordinatesOfBestStart.txt
# sed -i'' 's/Human_OR2J3 	//g' step_5_result.pickedM.CoordinatesOfBestStart.txt

# cat ./niimura.genbank.outgroupSeqs.fasta ./zebra_query.fasta ./step_5_result.pickedM.fasta > step_5_result.pickedM.wOutgroups.wRepresentatives.fasta 
# einsi --preservecase --thread 22 --inputorder step_5_result.pickedM.wOutgroups.wRepresentatives.fasta > step_6_result_mafftAlignment.wOutgroups.fasta &

# clustalw -CLUSTERING=NJ -bootstrap=1000 -KIMURA -TOSSGAPS -BOOTLABELS=node -infile=step_6_result_mafftAlignment.wOutgroups.fasta -outfile=step_6_result_mafftAlignment.wOutgroups.phb &

# R
# bedtools getfasta -s -name -bed ./{species}_2finalORCoordinates.all.0based.bed -fi /{NCBI genome} -fo ./{species}_step_10_getFinalGenomeCoordsAndSeqs_PerRuf.step_10_results.finalORntSeqs.all.fna

# iqtree -s step_6_result_mafftAlignment.wOutgroups.fasta -st AA -pre {species}_einsi.iqtree -m TEST -bb 1000 -alrt 1000 -nt 22 &