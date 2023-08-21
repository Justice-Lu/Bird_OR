import pandas as pd 
import numpy as np 
import re 


def read_fasta(filename):
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence_name = None
        sequence = ''
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    data.append((sequence_name, sequence))
                    sequence = ''
                sequence_name = line[1:]
            else:
                sequence += line
        if sequence_name is not None:
            data.append((sequence_name, sequence))
    data = pd.DataFrame(data, columns = ['name', 'sequence']) 
    return data

def get_domain_info(domain_reference: pd.DataFrame,
                    human_seq:str 
                    ):

    for q_motif in domain_reference.motif:
        
        pattern = '-*'+('-*'.join(list(q_motif)))
        matches = re.finditer(pattern, human_seq)
        
        for match in matches:
            start_pos = match.start()
            end_pos = match.end() - 1  # Adjust for inclusive end position    
        
            
        domain_reference.loc[domain_reference['motif'] == q_motif,['start']] = start_pos
        domain_reference.loc[domain_reference['motif'] == q_motif,['end']] = end_pos
    
    
    domain_reference[['start', 'end']] = domain_reference[['start', 'end']].astype(int)
    
    return domain_reference

# def calculate_indel(query, 
#                     domain_reference, 
#                     domain = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']):

#     # assert ['start', 'end'] in domain_reference.columns
    
#     for index, row in query.iterrows():
#         for domain in domains:
#             start = int(domain_reference[domain_reference['domain'] == domain].start)
#             end = int(domain_reference[domain_reference['domain'] == domain].end)
#             indel_count = int(row.sequence[start:end].count('-'))
#             query.loc[query['name'] == row['name'], domain] = indel_count
    
#     return query
        
    
def seq_mismatch(ref_seq: str, 
                 q_seq:str, 
                 consecutive_mismatch_threshold = 5
                ):
    
    consecutive_mismatch_count = 0
    mismatch = False
    for ref_char, q_char in zip(ref_seq, q_seq):
        if (ref_char == '-') and (q_char == '-'):
            continue

        if (ref_char != '-' and q_char == '-') or (ref_char == '-' and q_char != '-'):
            consecutive_mismatch_count += 1
            if (consecutive_mismatch_count >= consecutive_mismatch_threshold) | (q_char == '*'):
                mismatch = True    
                break
        else:
            consecutive_mismatch_count = 0

    return mismatch
        