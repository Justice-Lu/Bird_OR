# Bird_OR


## Tools 

### Filter Pseudogene
FilterPseudoGene class is a tool to iterate through a `.fasta` file and determine if the OR sequences have mismatches against Human_ref sequence based on pre-defined criteria: 
- Within a span of an transmembrane domain `<= 5` mismatches
- Query sequence should not have prematurely ended C-terminus

Example to call the script
```bash
python filter_pseudogene.py \
--save_fasta True \
--query_file_path ./pathtofile/MyQuery.fasta \
--domain_file_path ./pathtofile/DomainReference.txt \
```

## Steps

Step 1: Set up. Get genome. Create blast database. Create fai reference file.
```bash
wget {NCBI genome URL}
gunzip {NCBI genome}
samtools faidx {NCBI genome}
makeblastdb -in {NCBI genome} -out outputDB -dbtype 'nucl' -hash_index
```

Step 2: Blast. Input: blast query, genome. Output: blast results (.sum file)
```bash
perl Step_1_ab_tblastN_ORs.pl
```

Step 2b: If .sum file is too large (output from step 2), divide into 2 separate sum files. Input: blast results (.sum file). Output: two .sum files. Remember, paste header of .sum file to TAIL file.
```bash
head -{first half of file} {.sum} > HEAD.sum
tail -{second half of file} {.sum} > TAIL.sum
```

Step 3: Filter ORs using R. Input: blast results (.sum file). Output: bed file
```bash
R
```

Step 4a: If two bed files, determine non-overlapping values to add to one single bed file. Input: two bed files. Output: one bed file. 
```bash
bedtools intersect -v -a {HEAD.bed} -b {TAIL.bed} 
```

Step 4b: Convert bed file to fasta. Input: bed file. Output: fasta file
```bash
bedtools getfasta -s -name -bed {species}_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed -fi {NCBI genome} -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
```

Step 5: determine ORFs. 
```bash
perl step_3_a_Find_ORF-ABmod.pl
perl step_3_b_Find_ORF-ABmod_withPrintStatements.pl
perl step_3_d_GetORFCoords.pl
```

Step 6: Link to human OR, align. Input: "ORF" fasta. Output: alignment.
```bash
cat ./Human_OR2J3.fasta ./ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta > ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta
einsi --preservecase --thread 22 --inputorder ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta > {species}_step_4_result_mafftAlignment.fasta
```

Step 7: Run Justice script eliminating "pseudogenes"
```bash
python filter_pseudogene.py \
--save_fasta True \
--query_file_path ./pathtofile/MyQuery.fasta \
--domain_file_path ./pathtofile/DomainReference.txt \
```
Step 8: Edit step 5 perl scripts to include first transmembrane domain number. (Start of TM1 valine.) Then select appropriate start codon (perl script). From output: manually remove human OR sequence.
```bash
perl step_5_a_Start_M_pick_up.120.pl
perl step_5_d_WriteOutCoords.120.pl
sed -i'' 's/>//g' step_5_result.pickedM.CoordinatesOfBestStart.txt
sed -i'' 's/Human_OR2J3 	//g' step_5_result.pickedM.CoordinatesOfBestStart.txt
```

Step 9: Link to non-OR outgroups, zebra finch. Align again.
```bash
cat ./niimura.genbank.outgroupSeqs.fasta ./zebra_query.fasta ./step_5_result.pickedM.fasta > step_5_result.pickedM.wOutgroups.wRepresentatives.fasta 
einsi --preservecase --thread 22 --inputorder step_5_result.pickedM.wOutgroups.wRepresentatives.fasta > step_6_result_mafftAlignment.wOutgroups.fasta &
```
Step 10: Build NJ tree
```bash
clustalw -CLUSTERING=NJ -bootstrap=1000 -KIMURA -TOSSGAPS -BOOTLABELS=node -infile=step_6_result_mafftAlignment.wOutgroups.fasta -outfile=step_6_result_mafftAlignment.wOutgroups.phb &
```

Step 11: Get filtered bed file and nucleotide sequences.
```bash
R
bedtools getfasta -s -name -bed ./{species}_2finalORCoordinates.all.0based.bed -fi /{NCBI genome} -fo ./{species}_step_10_getFinalGenomeCoordsAndSeqs_PerRuf.step_10_results.finalORntSeqs.all.fna
```

Step 12: Build ML tree
```bash
iqtree -s step_6_result_mafftAlignment.wOutgroups.fasta -st AA -pre {species}_einsi.iqtree -m TEST -bb 1000 -alrt 1000 -nt 22 &
```




# TODO 
- CORP score for pseudogene (https://genome.weizmann.ac.il/horde/CORP/)
- Check if the common AA is presered 
- Domain_reference - get domain from alphafold 
- MSA alignment ? 
