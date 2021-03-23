File parameters.py
Please specify in this file next configurations before starting the analysis:
BC_length, min_BC_length and max_BC_length - length of barcodes (BCs) in your MPRA analysis (including minimum and maximum allowable BC length)
ROI_length, min_ROI_length and max_ROI_length - length of regions of interest (ROIs) in your MPRA analysis (including minimum and maximum allowable ROI length)
phred_quality - minimum allowable level of PHRED quality in your reads
index_dict - name (e.g. M1, E1, N1 etc) and sequences of your indexes (e.g. AGCGAGCT)
mapping_indexes - please, specify here indexes corresponding to your mapping samples (e.g. M1, M2)
norm_or_expr_indexes - please, specify here indexes corresponding to your normalization or expression samples (e.g. N1, E1)
library_list - names of libraries you want to analyze 
replicate_type_list - types of your MPRA samples, e.g. normalization, expression, mapping
neCP1, neCP2 and neCP2reference - please, specify here sequences of constant parts 1/2 of your normalization and expression samples, including reference sample 
mCP1, mCP2, mCP3 - please, specify here sequences of constant parts of your mapping samples
wt_ROI - sequence of ROIs in original (wild-type) sample
reference_BCs - sequences of BCs in reference samples. NOTE: if you don't have reference samples, please remove names of BCs as shown here: reference_BCs = {'forward': {}, 'reverse':{}}
reference_ROIs - sequences of ROIs in reference samples
index_mismatch_8bp, BC_mismatch and ROI_mismatch - allowable mismatch level in indexes, BCs and ROIs
neCP1_mismatch, neCP2_mismatch, mCP1_mismatch, mCP2_mismatch, mCP3_mismatch - allowable mismatch level in constant parts 1/2/3 of normalization, expression and mapping samples
neCP2_length and mCP3_length - length of constant parts 2/3 of normalization and expression samples
minimum_read_count - please, specify here minimum allowable read count for your normalization, expression and mapping samples
genuine_ROI_portion - minimum portion of ROI to be a genuine ROI associated with a genuine BC
ne_mode, map_mode - mode is reverse if BC sequences in samples are reversed, otherwise mode is forward
hash_length - the length of k-mer is calculated as: BC length // (BC_mismatch + 1)

To start analysis you need to put files parameters.py and MPRAdecoder.py in one folder and run MPRAdecoder.py.
Upon the end of the run please put files parameters.py and MPRAdecoder_stat.py in output directory and run MPRAdecoder_stat.py.
