#PARAMETERS

#Sequences of all indices, constant parts (CPs) and controls.

BC_length = 18
ROI_length = 8

index_dict = {
        'M1': 'AGCGAGCT', 'M2': 'CTGCACGT', 'N1': 'CCTATGGT', 'N2': 'AACGTCGT',
        'E1': 'ACAATTCG', 'E2': 'TACTTGTC'
        }
        
mapping_indices = [
        'M1', 'M2'
        ]
        
norm_or_expr_indices = [
        'N1', 'N2', 'E1', 'E2' 
        ]
        
library_list = [
        '29-36','37-44'
        ]
        
replicate_type_list = ['normalization', 'expression', 'mapping']

neCP1 = 'CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTC'
neCP2 = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGC'
neCP2control = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCG'

mCP1 = 'GACACTCGAGGATCGAG'
mCP2 = {
        '29-36':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATT',
        '37-44':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCAT',
            }
mCP3 = {
        '29-36': 'GATTATCTTTAACGTACGTCACAAT', '37-44': 'TTAACGTACGTCACAAT'
        }

wt_ROI = {
        '29-36': 'TTACGCAT', '37-44': 'GATTATCT'
        }

control_BCs = {
        'forward': {
                'bc_1': 'TTCCAAGTGCAGGTTAGGCG', 'bc_2': 'TGTGTACGGCTTGCTCTCAA',
                'bc_3': 'GAGCCCGGATCCACTCCAAG', 'bc_4': 'TGTCACGTCAGCTAACCCAC'
                },
        'reversed':{
                'bc_1': 'CGCCTAACCTGCACTTGGAA', 'bc_2': 'TTGAGAGCAAGCCGTACACA',
                'bc_3': 'CTTGGAGTGGATCCGGGCTC', 'bc_4': 'GTGGGTTAGCTGACGTGACA'
                }
             }

control_ROIs = {'wt': 'TTACGCAT', 'del_c': 'TTAGCATG'}

#Default mismatch levels, lengths of last CPs, minimum read counts.

index_mismatch_8bp = 1.1
BC_mismatch = 2.2
ROI_mismatch = 1.1

neCP1_mismatch = 4.4
neCP2_mismatch = 2.2
mCP1_mismatch = 2.2

neCP2_length = 20

mCP2_mismatch = {
        '29-36': 8.8, '37-44': 8.8
        }
mCP3_length = {
        '29-36': 10, '37-44': 10
        }
mCP3_mismatch = {
        '29-36': 1.1, '37-44': 1.1
        }


minimum_read_count = {'normalization': 3, 'expression': 0}
major_ROI_portion = 0.9
