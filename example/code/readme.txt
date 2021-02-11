File parameters.py

BC_length - length of barcodes (BCs) in your MPRA analysis
ROI_length - length of regions of interest (ROIs) in your MPRA analysis
index_dict - name (e.g. M1, E1, N1 etc) and sequences of your indices (e.g. AGCGAGCT)
mapping_indices - specify here indices corresponding to your mapping samples (e.g. M1, M2)
norm_or_expr_indices - specify here indices corresponding to your normalization or expression samples (e.g. N1, E1)
library_list - names of libraries you want to analyze
replicate_type_list - types of your MPRA samples
neCP1, neCP2 and neCP2control - specify here sequences of constant parts 1/2 of your normalization and expression samples, including control sample
mCP1, mCP2, mCP3 - specify here sequences of constant parts of your mapping samples
wt_ROI - sequence of ROI in original (wild-type) sample
control_BCs - sequences of BCs in control samples
index_mismatch_8bp, BC_mismatch and ROI_mismatch - allowable mismatch level in indices, BCs and ROIs
neCP1_mismatch, neCP2_mismatch, mCP1_mismatch, mCP2_mismatch, mCP3_mismatch - allowable mismatch level in constant parts 1/2/3 of your normalization, expression and mapping samples
neCP2_length and mCP3_length - length of constant parts 2/3 of your normalization and expression samples
minimum_read_count - please specify minimum read count for your normalization and expression samples
major_ROI_portion - min portion of ROI to be a major ROI associated with a certain BC
