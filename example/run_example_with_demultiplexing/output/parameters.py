#PARAMETERS
#Sequences of all indexes, constant parts (CPs), 1-0s and PHRED quality.
BC_length = 18
min_BC_length = 16
max_BC_length = 20

ROI_length = 8
min_ROI_length = 7
max_ROI_length = 9

phred_quality = 10

index_dict = {'M1': 'AGCGAGCT', 'M2': 'CTGCACGT', 'N1': 'CCTATGGT', 'N2': 'AACGTCGC',
        'E1': 'ACAATTCG', 'E2': 'TACTTGTC'}

mapping_indexes = ['M1', 'M2']

norm_or_expr_indexes = ['N1', 'N2', 'E1', 'E2']

library_list = ['ref-libr','example-libr']

replicate_type_list = ['normalization', 'expression', 'mapping']

neCP1 = 'CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTC'
neCP2 = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGC'
neCP2reference = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCG'

mCP1 = 'GACACTCGAGGATCGAG'
mCP2 = {'ref-libr':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATT',
        'example-libr':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCAT'}

mCP3 = {'ref-libr': 'GATTATCTTTAACGTACGTCACAAT', 'example-libr': 'TTAACGTACGTCACAAT'}

wt_ROI = {'ref-libr': 'TTACGCAT', 'example-libr': 'GATTATCT'}

reference_BCs = {'forward': {
                'control1': 'TTCCAAGTGCAGGTTAGGCG', 'control2': 'TGTGTACGGCTTGCTCTCAA',
                'experiment1': 'GAGCCCGGATCCACTCCAAG', 'experiment2': 'TGTCACGTCAGCTAACCCAC'},
        'reverse':{
                'control1': 'CGCCTAACCTGCACTTGGAA', 'control2': 'TTGAGAGCAAGCCGTACACA',
                'experiment1': 'CTTGGAGTGGATCCGGGCTC', 'experiment2': 'GTGGGTTAGCTGACGTGACA'}}
#If you don't have reference samples, please remove names of BCs as shown below:
#reference_BCs = {'forward': {}, 'reverse':{}}

reference_ROIs = {'control': 'TTACGCAT', 'experiment': 'TTAGCATG'}

#Default mismatches levels, lengths of last CPs, minimum read counts.

index_mismatch_8bp = 1
BC_mismatch = 2
ROI_mismatch = 0

neCP1_mismatch = 4
neCP2_mismatch = 2
mCP1_mismatch = 2

neCP2_length = 20

mCP2_mismatch = {'ref-libr': 7, 'example-libr': 8}

mCP3_length = {'ref-libr': 10, 'example-libr': 10}

mCP3_mismatch = {'ref-libr': 1, 'example-libr': 1}

minimum_read_count = {'normalization': 3, 'expression': 0, 'mapping': 1}

genuine_ROI_portion = 0.9

#ne_mode = 'reverse' if BC sequences in normalization and expression samples are reversed, otherwise ne_mode = 'forward'
ne_mode = 'reverse' 
#ne_mode = 'forward'

#map_mode = 'forward' if BC and ROI sequences in mapping samples are forward, otherwise map_mode = 'reverse'
map_mode = 'forward'
#map_mode = 'reverse' 

hash_length = 6