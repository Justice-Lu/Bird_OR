# Bird_OR


## Tools 

### Filter Pseudogene
FilterPseudoGene class is a tool to iterate through a `.fasta` file and determine if the OR sequences have mismatches against Human_ref sequence based on pre-defined criteria: 
- Within a span of an transmembrane domain `<= 5` mismatches
- Query sequence should not have prematurely ended C-terminus

Exmaple to call the script
```bash
python filter_pseudogene.py \
--save_fasta True \
--query_file_path ./pathtofile/MyQuery.fasta \
--domain_file_path ./pathtofile/DomainReference.txt \
```
