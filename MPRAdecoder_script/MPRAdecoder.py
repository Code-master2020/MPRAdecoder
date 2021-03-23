#!/usr/bin/env python
# -*- coding: utf-8 -*-


import gzip
import os
from Bio import SeqIO, pairwise2, motifs
from Bio.Seq import Seq
from datetime import datetime
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import json
from random import randint
import math
import seaborn as sns
import parameters as param



# Some functions to collect experiment information.
#Ask user how many files or libraries will be analyzed.
def input_count(name, additional):
	input_number = input(
			'How many %s%s would you like to analyze? Type number: '
			%(name, additional)
			)
	try:
		response = int(input_number)
		print('%s %s will be analyzed.' %(input_number, name))
		return response
	except ValueError:
		print('You entered an incorrect value: %s. Please try again.'
				% input_number
				)
		return input_count(name, additional)

# Check if the path to file or library name exist.
def check_name(n, name_2):
	if n == 0:
		next_name = input('Enter first %s: ' % name_2)
	elif n == 1:
		next_name = input('Enter second %s: ' % name_2)
	elif n == 2:
		next_name = input('Enter third %s: ' % name_2)
	else:
		next_name = input('Enter next %s: ' % name_2)

	if name_2 == 'absolute path to file':
		if next_name.endswith('.fastq.gz'):
			try:
				with gzip.open(next_name, 'rt'):
					pass
			except FileNotFoundError:
				print('No such file or directory: %s. Please, try again.'
						% next_name
						)
				check_name(n, name_2)
			else:
				return next_name
		else:
			print('Wrong file format. File name should end with .fastq.gz.\nPlease, try again.')
			return check_name(n, name_2)


	elif name_2 == 'library name':
		if next_name in param.library_list:
			return next_name
		else:
			print('No such library: %s. \nPlease, try again.'
					% next_name
					)
			return check_name(n, name_2)

# Validate index.
def check_index(
		rep_type, index_name_list, file_dictionary, name_file, libr_dictionary,
		libr_name
		):
	valid_indexes = []
	for index_name in index_name_list:
		if index_name in param.index_dict:
			if (index_name in param.mapping_indexes and rep_type == 'mapping') or (index_name in param.norm_or_expr_indexes and rep_type == 'expression') or (index_name in param.norm_or_expr_indexes and rep_type == 'normalization'):
				if index_name not in file_dictionary[name_file]:
					file_dictionary[name_file][index_name] = {'read count': 0}
					valid_indexes.append(index_name)
					print('Index %s was successfully added.' % index_name)
				else:
					print('Index %s was already used in file %s.' %(index_name, name_file))
			else:
				print('Index %s is not appropriate for %s.' % (index_name, rep_type))
		else:
			print('There is no such index: %s.' % index_name)
	if valid_indexes:
		libr_dictionary[libr_name][rep_type][name_file].append(valid_indexes)
	return file_dictionary, libr_dictionary

# Continue or restart step.
def yes_or_no():
	answer = input(
			'Is everything correct? Type Y to continue or N to restart previous step. '
			).upper()
	if answer == 'Y' or answer == 'N':
		return answer
	else:
		print('You entered an incorrect value: %s. Try again.' % answer)
		return yes_or_no()

# Form dictionary of input files.
def add_file(name_plural, addit, name_single):
	file_dict = {}
	file_cout = input_count(name_plural, addit)

	for i in range(file_cout):
		new_file_name = check_name(i, name_single)
		if new_file_name not in file_dict:
			file_dict[new_file_name] = {}
		else:
			print('File %s has already been added.' % new_file_name)

	print('You added %s files:' % str(len(file_dict)))
	for added_file in file_dict:
		print(added_file)

	ans = yes_or_no()
	if ans == 'N':
		add_file(name_plural, addit, name_single)
	else:
		return file_dict
	
# Add path to output file, make a text file.
def add_out_path():
	out_path = input(
			'Enter absolute path to your output folder (e.g. /mnt/d/User/output/): '
			)
	new_results = out_path
	run_info = new_results + 'run_info.txt'
	try:
		os.makedirs(new_results, exist_ok=True)
		with open(run_info, 'w') as file:
			file.write('Hello, user!\n')
		return new_results, run_info
	except OSError:
		print('Creation of the directory %s failed' % out_path)
		return add_out_path()

# Check presence of normalization, expression or mapping replicates. 
def have_any(rep_name, lib_name, file_name):
	answ = input(
			'Do you have any %s replicates for library %s in file %s? Type Y or N. '
			%(rep_name, lib_name, file_name)
			)
	if answ.upper() == 'Y' or answ.upper() == 'N':
		return answ.upper()
	else:
		print('You entered an incorrect value: %s. Try again.' % answ)
		return have_any(rep_name, lib_name, file_name)		

# Form dictionary of libraries in format: {
#	'29-36': {'normalization': {1.fastq.gz: [[A11, A13], [A15]]},
#			  'expression': {1.fastq.gz: [[A12], [A14]]},
#			  'mapping': {1.fastq.gz: [[A1]]},
#			 }
#	 'H' : {1.fastq.gz: [A2, A3, A4]}
#	}
def add_lib(name_plur, file_diction, name_sing):
	library_dict = {}

	for file in file_diction:
		lib_count = input_count(name_plur, ' in ' + file)

		for _ in range(lib_count):
			library_name = check_name(_, name_sing)

			if library_name != 'H':
				if library_name not in library_dict:
					library_dict[library_name] = {
							'normalization': {}, 'expression': {},
							'mapping': {}
							}

				for replicate_type in param.replicate_type_list:
					response = have_any(replicate_type, library_name, file)
					if response == 'Y':
						library_dict[library_name][replicate_type][file] = []
						print('Enter indexes for %s replicates of library %s in file %s.\nIf there are more than one index per replicate, print one index per line. Finally type N.'
								%(replicate_type, library_name, file
								  ))
						ind_list = input().upper().strip().split()
						while ind_list != ['N']:
							file_diction, library_dict = check_index(
									replicate_type, ind_list, file_diction,
									file, library_dict, library_name
									)
							ind_list = input().upper().strip().split()
			else:
				print('Enter indexes for chimera test in file %s.\nIf there are more than one index per replicate, print one index per line. Finally type N.'
						% file
						)
				if library_name not in library_dict:
					library_dict['H'] = {file: []}
				else:
					library_dict['H'][file] = []

				ind = input().upper().strip()
				while ind != 'N':
					if ind in param.mapping_indexes:
						if ind not in file_diction[file]:
							file_diction[file][ind] = {'read count': 0}
							library_dict['H'][file].append(ind)
							print('Index %s was successfully added.' % ind)
						else:
							print('Index %s was already used in file %s.'
								  %(ind, file)
								  )
					else:
						print('Index %s is not appropriate for chimera test.'
							  % ind
							  )
					ind = input().upper().strip()
	
	for lib_name in library_dict:
		if library_dict[library_name] == {
				'normalization': {}, 'expression': {}, 'mapping': {}
				}:
			raise Exception('Error! Indices for library %s was not added!'
							% lib_name)

	print('You added %s libraries:' % str(len(library_dict)))
	for library in library_dict:
		print('Library %s.' % library)
		if library != 'H':
			for rep_type in library_dict[library]:
				for file in library_dict[library][rep_type]:
					for replicate in library_dict[library][rep_type][file]:
						print('{0} replicate with indexes {1} in file {2}.'.format(
								rep_type.capitalize(), replicate, file
								))
		else:
			for file in library_dict[library]:
				print('Chimera test indexes {0} in file {1}.'.format(
						library_dict[library][file], file
						))
					
	resp = yes_or_no()
	if resp == 'Y':
		return file_diction, library_dict
	else:
		for file in file_diction:
			file_diction[file] = {}
		return add_lib(name_plur, file_diction, name_sing)

# Asks the user if he wants to change the level of mismatches.
def change_or_not():
	ans = input('Would you like to change something? Type Y or N: ').upper()
	if ans == 'Y' or ans == 'N':
		return ans
	else:
		print('You entered an incorrect value: %s. Try again.' % ans)
		return change_or_not()

def mismatch_max(input_message, max_value, libr_name=''):
	try:
		current_value = int(
				input('Enter ' + input_message + libr_name + ': ')
				)

		if current_value > max_value:
			print('Max ' + input_message + libr_name + ' is ' + str(max_value))
			current_value = max_value
		return current_value
	except ValueError:
		print('Invalid input for %s. Only integers are acceptable! Try again.'
			  % input_message)
		return mismatch_max(input_message, max_value, libr_name='')


#Chenges values in dictionary minimum_read_count
def change_read_count(rep):
	try:
		min_read_count = int(input('Enter minimum read count for ' + rep + ': '))
		if rep == 'normalization' and min_read_count <= 0:
			print('Minimum read count for normalization can not be equal 0. Please, try again.')
			change_read_count(rep)
		return min_read_count
	except ValueError:
		print('You entered incorrect value. Try again.')
		change_read_count(rep)


#Change the level of mismatches.
def change_mismatches(libr_dict):
	param.index_mismatch_8bp, param.BC_mismatch, param.ROI_mismatch, param.neCP1_mismatch, param.neCP2_mismatch, param.neCP2_length, param.mCP1_mismatch, param.mCP2_mismatch, param.mCP3_length, param.mCP3_mismatch, param.minimum_read_count
	
	param.index_mismatch_8bp = mismatch_max(
			'mismatch value for indexes', 4
			)
	param.BC_mismatch = mismatch_max(
			'mismatch value for BCs', param.BC_length // 2 - 1
			)
	param.ROI_mismatch = mismatch_max(
			'mismatch value for ROIs', param.ROI_length // 2 - 1
			)
	param.neCP1_mismatch = mismatch_max(
			'mismatch value for neCP1',
			len(param.neCP1) // 2
			)
	param.neCP2_length = mismatch_max(
			'length of neCP2',
			len(param.neCP2reference)
			)
	param.neCP2_mismatch = mismatch_max(
			'mismatch value for neCP2',
			param.neCP2_length // 2
			)
	param.mCP1_mismatch = mismatch_max(
			'mismatch value for mCP1',
			len(param.mCP1) // 2
			)
		
	for libr in libr_dict:
		if libr == 'H':
			param.mCP2_mismatch['25-32'] = mismatch_max(
					'mismatch value for chimera test CP2',
					len(param.mCP2['25-32']) // 2
					)
			param.mCP3_length['25-32'] = mismatch_max(
					'length of CP3 for chimera test',
					len(param.mCP3['25-32'])
					)
			param.mCP3_mismatch['25-32'] = mismatch_max(
					'mismatch level for chimera test CP3',
					param.mCP3_length['25-32'] // 2
					)
		else:
			if  libr_dict[libr]['mapping']:
				param.mCP2_mismatch[libr] = mismatch_max(
						'mismatch value for mCP2 in library ',
						len(param.mCP2[libr]) // 2, libr
						)
				
				param.mCP3_length[libr] = mismatch_max(
						'length of CP3 in library ',
						len(param.mCP3[libr]), libr
						)
				param.mCP3_mismatch[libr] = mismatch_max(
						'mismatch value for mCP3 in library ',
						param.mCP3_length[libr] // 2, libr
						)				   

	param.mCP2_mismatch['ref-libr'] = mismatch_max(
			'mismatch level for reference mCP2',
			len(param.mCP2['ref-libr']) // 2
			)
	param.mCP3_length['ref-libr'] = mismatch_max(
			'length of reference mCP3',
			len(param.mCP3['ref-libr'])
			)
	param.mCP3_mismatch['ref-libr'] = mismatch_max(
			'mismatch level for reference mCP3',
			param.mCP3_length['ref-libr'] // 2
			)
	
	for rep_type in param.minimum_read_count:
		param.minimum_read_count[rep_type] = change_read_count(rep_type)
	
	report = yes_or_no()
	if report == 'N':
		change_mismatches(libr_dict)		

# Display and change mismatches levels
def mismatches_count(libr_dict, run_inf_text):
	print('Please, check acceptable mismatch levels.\n{0} mismatches in indexes.'.format(
			str(param.index_mismatch_8bp)
			))
	print('{0} mismatches in BCs, {1} in ROIs.'.format(
			str(param.BC_mismatch), str(param.ROI_mismatch)
			))
	print('{0} mismatches in neCP1;\n{1} mismatches in neCP2.'.format(
			str(param.neCP1_mismatch), str(param.neCP2_mismatch)
			))
	print('{0} mismatches in mCP1'.format(
			str(param.mCP1_mismatch)
			))

	for lib in libr_dict:
		if lib == 'H':
			print('For chimera formation test {0} mismatches for mCP2;\n{1} mismatches for chimera test CP3 {2} bp long.'.format(
					str(param.mCP2_mismatch['25-32']),
					str(param.mCP3_mismatch['25-32']),
					str(param.mCP3_length['25-32'])
					))
		elif libr_dict[lib]['mapping']:
			print('For library {0} {1} mismatches in mCP2;\n{2} mismatches for mCP3.'.format(
					lib, str(param.mCP2_mismatch[lib]),
					str(param.mCP3_mismatch[lib])
					))

	if len(param.reference_BCs['forward']) > 0:
		print('For mapping reference {0} mismatches for mCP2;\n{1} mismatches for mCP3.'.format(
				str(param.mCP2_mismatch['ref-libr']),
				str(param.mCP3_mismatch['ref-libr'])
				))

	for rep_type in param.minimum_read_count:
		print('Minimum read count for ' + rep_type + ' sample is ' +
			  str(param.minimum_read_count[rep_type]) + '.')

	respons = change_or_not()

	if respons == 'Y':
		change_mismatches(libr_dict)
		
	lines_of_text = [
			'{0} mismatch in indexes.\n'.format(str(param.index_mismatch_8bp)),
			'{0} mismatches in BCs, {1} in ROIs.\n'.format(str(param.BC_mismatch), str(param.ROI_mismatch)),
			'{0} mismatches in neCP1;\n{1} mismatches in neCP2 {2} bp long.\n'.format(str(param.neCP1_mismatch), str(param.neCP2_mismatch), str(param.neCP2_length)),
			'{0} mismatches in mCP1\n'.format(str(param.mCP1_mismatch)),
			]
	for lib in libr_dict:
		if lib == 'H':
			lines_of_text.append(
					'For chimera test {0} mismatches for mCP2;\n{1} mismatches for chimera test CP3 {2} bp long.\n'.format(str(param.mCP2_mismatch['25-32']), str(param.mCP3_mismatch['25-32']), str(param.mCP3_length['25-32']))
					)
		elif libr_dict[lib]['mapping']:
			lines_of_text.append(
					'For library {0} {1} mismatches in mCP2;\n{2} mismatches for mCP3 {3} bp long.\n'.format(lib, str(param.mCP2_mismatch[lib]), str(param.mCP3_mismatch[lib]), str(param.mCP3_length[lib]))
					)
	with open(run_inf_text, 'a') as info_file:
		info_file.writelines(lines_of_text)

# 16 Form dictionary of barcodes and ROIs in format:
#{'1.fastq.gz':
#			  {A11: {'read count': 1}, {'eff read count': 1}, {'bc1': ['seq1']}},
#			  {A13: {'read count': 1}, {'eff read count': 1}, {'bc2': ['seq2']}},
#			  {A15: {'read count': 0}} {'eff read count': 0},,
#			  {A12: {'read count': 1}, {'eff read count': 1}, {'bc3': ['seq3']}},
#			  {A14: {'read count': 1}, {'eff read count': 1}, {'bc4': ['seq4']}},
#			  {A1: {'read count': 1}, {'eff read count': 1}, {'bc4': {'mut4': ['seq8']}}},
#			  {A2: {'read count': 1}, {'eff read count': 1}, {'bc1': {'mut1': ['seq5']}}},
#			  {A3: {'read count': 1}, {'eff read count': 1}, {'bc2': {'mut2': ['seq6']}}},
#			  {A4: {'read count': 1}, {'eff read count': 1}, {'bc3': {'mut3': ['seq7']}}}
#}

# Greeting user and collect experiment information.
def collect_info():
	print('Hello, user!')
	
	fastq_file_dict = add_file('FASTQ files', '', 'absolute path to file')
	
	output_path, output_file = add_out_path()
	
	fastq_file_dict, lib_dict = add_lib(
			'libraries', fastq_file_dict, 'library name'
			)

	mismatches_count(lib_dict, output_file)
	
	with open(output_path + '/' + 'fastq_file_dictionary.json', 'w+') as f:
		json.dump(fastq_file_dict, f)
		
	with open(output_path + '/' + 'lib_diction.json', 'w+') as file:
		json.dump(lib_dict, file)

	return fastq_file_dict, lib_dict, output_path, output_file


fastq_file_dictionary, lib_diction, output_folder, run_info_text = collect_info()

# Align two sequences
def count_score(query, ref):
	return pairwise2.align.globalxs(query, ref, -1, -1, score_only=True)

# Find the best alignment
def aligner(read, reference, mismatch):
	insertion = 0
	deletion = 0
	max_score = count_score(read[:len(reference)], reference)
	best_align_length = len(reference)
	if mismatch > 0:
		while insertion < mismatch:
			insertion += 1
			insertion_score = count_score(
					read[:len(reference) + insertion], reference
					)
			if insertion_score > max_score:
				max_score = insertion_score
				best_align_length = len(reference) + insertion
			else:
				insertion -= 1
				break
		if best_align_length == len(reference):
			while deletion < mismatch:
				deletion += 1
				deletion_score = count_score(
						read[:len(reference) - deletion], reference
						)
				if deletion_score >= max_score:
					max_score = deletion_score
					best_align_length = len(reference) - deletion
				else:
					break
	minimum_score = len(reference) - mismatch + insertion
	return max_score, best_align_length, minimum_score
	if mismatch == 0:
		minimum_score = len(reference)
	return max_score, best_align_length, minimum_score
	


# Check if read quality upper than 10.
def check_quality(seq, start, end):
	quality_list = seq.letter_annotations['phred_quality']
	for quality in range(start, end):
		if quality_list[quality] < param.phred_quality:
			return False
	return True

# Find barcodes
def find_barcode(
		rec, start, const1, const1_mismatch, const2, const2_mismatch, mode
		):
	score1, const1_length, min_score1 = aligner(
			rec.seq[start:], const1, const1_mismatch
			)
	if type(score1) == float and score1 >= min_score1:
		for control_barcode in param.reference_BCs[mode]:
			if len(param.reference_BCs[mode]) > 0:
				score_conbc, conbc_length, min_score_conbc = aligner(
						rec.seq[(start + const1_length):],
						param.reference_BCs[mode][control_barcode], param.BC_mismatch
						)
				if score_conbc >= min_score_conbc:
					break

		else:
			best_score2, const2_length, best_min_score2 = aligner(
				rec.seq[(start + const1_length + param.BC_length):],
				const2, const2_mismatch
				)
			best_BC_length = param.BC_length
			
			insertion = 0
			deletion = 0
			while True:
				insertion += 1
				ins_score2, const2_length, min_score2 = aligner(
						rec.seq[(start + const1_length + param.BC_length + insertion):],
						const2, const2_mismatch
						)
				if ins_score2 > best_score2:
					best_score2 = ins_score2
					best_min_score2 = min_score2
					best_BC_length = param.BC_length + insertion
				else:
					break
					
			if best_BC_length == param.BC_length:
				while True:
					deletion += 1
					del_score2, const2_length, min_score2 = aligner(
							rec.seq[(start + const1_length + param.BC_length - deletion):],
							const2, const2_mismatch
							)
					if del_score2 > best_score2:
						best_score2 = del_score2
						best_min_score2 = min_score2
						best_BC_length = param.BC_length - deletion
					else:
						break
			if best_BC_length <= param.BC_length + param.BC_mismatch and best_BC_length >= param.BC_length - param.BC_mismatch:
				if best_score2 >= best_min_score2:
					if check_quality(rec, (start + const1_length), (start + const1_length + best_BC_length)):
						return rec.seq[start + const1_length:start + const1_length + best_BC_length]
					return 'low quality'
		if len(param.reference_BCs['forward']) > 0 and check_quality(rec, (start + const1_length), (start + const1_length + conbc_length)):
			score2, const2_length, min_score2 = aligner(
					rec.seq[(start + const1_length + conbc_length):],
					const2, const2_mismatch
					)
			if score2 >= min_score2:
				return control_barcode
		return 'low quality'
	return 'undetermined'

#Find barcodes and ROIs for chimera test	
def hybrid_test_analysis(rec, start):
	seq = rec.seq
	score1, const1_length, min_score1 = aligner(
			seq[start:], param.mCP1, param.mCP1_mismatch
			)
		
	if type(score1) == float and score1 >= min_score1:
		for hybrid_test_barcode in param.hybrid_test_barcodes:
			score_hbc, hbc_length, min_score_hbc = aligner(
					seq[start + const1_length:],
					param.hybrid_test_barcodes[hybrid_test_barcode],
					param.BC_mismatch
					)
			if score_hbc >= min_score_hbc:
				break
		else:
			return 'undetermined', 'undetermined'
		if check_quality(
				rec,
				(start + const1_length),
				(start + const1_length + hbc_length)
				):
			if hybrid_test_barcode == 'bc15':
				score2, const2_length, min_score2 = aligner(
						seq[start + const1_length + hbc_length:],
						param.mCP2['hybrid_test'],
						param.mCP2_mismatch['25-32']
						)
			else:
				score2, const2_length, min_score2 = aligner(
						seq[start + const1_length + hbc_length:],
						param.mCP2['25-32'],
						param.mCP2_mismatch['25-32']
						)
			if score2 >= min_score2:
				for hybrid_test_ROI in param.hybrid_test_ROIs:
					score_hmut, hmut_length, min_score_hmut = aligner(
							seq[start + const1_length + hbc_length + const2_length:],
							param.hybrid_test_ROIs[hybrid_test_ROI],
							param.ROI_mismatch
							)
					if score_hmut >= min_score_hmut:
						break
				else:
					return hybrid_test_barcode, 'undetermined'
				if check_quality(
						rec,
						(start + const1_length + hbc_length + const2_length),
						(start + const1_length + hbc_length + const2_length + hmut_length)
						):
					score3, const3_length, min_score3 = aligner(
							seq[start + const1_length + hbc_length + const2_length + hmut_length:],
							param.mCP3['25-32'][:param.mCP3_length['25-32']],
							param.mCP3_mismatch['25-32']
							)
					if score3 >= min_score3:
						return hybrid_test_barcode, hybrid_test_ROI
				return hybrid_test_barcode, 'low quality'
			return hybrid_test_barcode, 'undetermined'
		return 'low quality', 'undetermined'
	return 'undetermined', 'undetermined' 

# Find BCs and ROIs for mapping
def find_barcode_ROI(rec, ind, start, lib_dictionar, file, mode):
	seq = rec.seq
	for lib in lib_dictionar:
		if lib == 'H':
			if file in lib_dictionar['H']:
				if ind in lib_dictionar['H'][file]:
					barcode, ROI = hybrid_test_analysis(rec, start)
					return barcode, ROI
		else:
			if file in lib_dictionar[lib]['mapping']:
				for index_list in lib_dictionar[lib]['mapping'][file]:
					if ind in index_list:
						librar = lib
						break
				else:
					continue
				break
	score1, const1_length, min_score1 = aligner(
			seq[start:], param.mCP1, param.mCP1_mismatch
			)
	if type(score1) == float and score1 >= min_score1:
		for control_barcode in param.reference_BCs['forward']:
			if len(param.reference_BCs['forward']) > 0:
				score_conbc, conbc_length, min_score_conbc = aligner(
						seq[(start + const1_length):],
						param.reference_BCs['forward'][control_barcode],
						param.BC_mismatch
						)
				if score_conbc >= min_score_conbc:
					cont_barcode = control_barcode
					break
		else:
			best_score2, best_const2_length, best_min_score2 = aligner(
					seq[start + const1_length + param.BC_length:],
					param.mCP2[librar], param.mCP2_mismatch[librar]
					)
			best_BC_length = param.BC_length
			

			insertion = 0
			deletion = 0
			while True:
				insertion += 1
				ins_score2, const2_length, min_score2 = aligner(
						seq[start + const1_length + param.BC_length + insertion:],
						param.mCP2[librar],
						param.mCP2_mismatch[librar]
						)
				if ins_score2 > best_score2:
					best_score2 = ins_score2
					best_const2_length = const2_length
					best_min_score2 = min_score2
					best_BC_length = param.BC_length + insertion
				else:
					break
			if best_BC_length == param.BC_length:
				while True:
					deletion += 1
					del_score2, const2_length, min_score2 = aligner(
							seq[(start + const1_length + param.BC_length - deletion):],
							param.mCP2[librar],
							param.mCP2_mismatch[librar]
							)
					if del_score2 > best_score2:
						best_score2 = del_score2
						best_const2_length = const2_length
						best_min_score2 = min_score2
						best_BC_length = param.BC_length - deletion
					else:
						break
			if best_BC_length <= param.BC_length + param.BC_mismatch and best_BC_length >= param.BC_length - param.BC_mismatch:
				if best_score2 >= best_min_score2:
					if check_quality(
							rec, (start + const1_length),
							(start + const1_length + best_BC_length)
							):
						best_score3, const3_length, best_min_score3 = aligner(
								seq[start + const1_length + best_BC_length + best_const2_length + param.ROI_length:],
								param.mCP3[librar][:param.mCP3_length[librar]],
								param.mCP3_mismatch[librar]
								)
						best_ROI_length = param.ROI_length
						insertion = 0
						deletion = 0
						while True:
							insertion += 1
							ins_score3, const3_length, min_score3 = aligner(
									seq[start + const1_length + best_BC_length + best_const2_length + param.ROI_length + insertion:],
									param.mCP3[librar][:param.mCP3_length[librar]],
									param.mCP3_mismatch[librar]
									)
							if ins_score3 > best_score3:
								best_score3 = ins_score3
								best_min_score3 = min_score3
								best_ROI_length = param.ROI_length + insertion
							else:
								break
					
						if best_ROI_length == param.ROI_length:
							while True:
								deletion += 1
								del_score3, const3_length, min_score3 = aligner(
										seq[start + const1_length + best_BC_length + best_const2_length + param.ROI_length - deletion:],
										param.mCP3[librar][:param.mCP3_length[librar]],
										param.mCP3_mismatch[librar]
										)
								if del_score3 > best_score3:
									best_score3 = del_score3
									best_min_score3 = min_score3
									best_ROI_length = param.ROI_length - deletion
								else:
									break
						if best_ROI_length >= param.ROI_length - param.ROI_mismatch and best_ROI_length <= param.ROI_length + param.ROI_mismatch:
							if best_score3 >= best_min_score3:
								if check_quality(
										rec,
										(start + const1_length + best_BC_length + best_const2_length),
										(start + const1_length + best_BC_length + best_const2_length + best_ROI_length)
										):
									return seq[
											start + const1_length:start + const1_length + best_BC_length
											], seq[start + const1_length + best_BC_length + best_const2_length:start + const1_length + best_BC_length + best_const2_length + best_ROI_length]
								return 'low quality ROI', 'low quality ROI'
						return 'undetermined', 'undetermined'
					return 'low quality barcode', 'low quality barcode'
			return 'undetermined', 'undetermined'
		if len(param.reference_BCs['forward']) > 0 and check_quality(
				rec, (start + const1_length),
				(start + const1_length + conbc_length)
				):
			score2, const2_length, min_score2 = aligner(
					seq[(start + const1_length + conbc_length):],
					param.mCP2['ref-libr'],
					param.mCP2_mismatch['ref-libr']
					)
			if score2 >= min_score2:
				if cont_barcode == 'control1' or cont_barcode == 'control2':
					score_cmut, cmut_length, min_score_cmut = aligner(
							seq[start + const1_length + conbc_length + const2_length:],
							param.reference_ROIs['control'],
							param.ROI_mismatch
							)
					cont_ROI = 'control'
				else:
					score_cmut, cmut_length, min_score_cmut = aligner(
							seq[start + const1_length + conbc_length + const2_length:],
							param.reference_ROIs['experiment'],
							param.ROI_mismatch
							)
					cont_ROI = 'experiment'
				if score_cmut >= min_score_cmut:
					if check_quality(
							rec,
							(start + const1_length + conbc_length + const2_length),
							(start + const1_length + conbc_length + const2_length + cmut_length)):
						score3, const3_length, min_score3 = aligner(
								seq[start + const1_length + conbc_length + const2_length + cmut_length:],
								param.mCP3['ref-libr'][:param.mCP3_length['ref-libr']],
								param.mCP3_mismatch['ref-libr']
								)
						if score3 >= min_score3:
							return cont_barcode, cont_ROI
					return 'low quality ROI', 'low quality ROI'
			return 'undetermined', 'undetermined'
		return 'low quality barcode', 'low quality barcode'					
	return 'undetermined', 'undetermined'

   
#Form dictionary of barcodes and ROIs in format:
#{'1.fastq.gz':
#			  {A11: {'read count': 1}, {'eff read count': 1}, {'bc1': ['seq1']}},
#			  {A13: {'read count': 1}, {'eff read count': 1}, {'bc2': ['seq2']}},
#			  {A15: {'read count': 0}} {'eff read count': 0},,
#			  {A12: {'read count': 1}, {'eff read count': 1}, {'bc3': ['seq3']}},
#			  {A14: {'read count': 1}, {'eff read count': 1}, {'bc4': ['seq4']}},
#			  {A1: {'read count': 1}, {'eff read count': 1}, {'bc4': {'mut4': ['seq8']}}},
#			  {A2: {'read count': 1}, {'eff read count': 1}, {'bc1': {'mut1': ['seq5']}}},
#			  {A3: {'read count': 1}, {'eff read count': 1}, {'bc2': {'mut2': ['seq6']}}},
#			  {A4: {'read count': 1}, {'eff read count': 1}, {'bc3': {'mut3': ['seq7']}}}
#}

def fastq_parsing(ind_dictionary, lib_dictionary):
	start_time = datetime.now()
	print('Analysis was started ' + str(start_time))
	for inp_file in ind_dictionary:
		undetermined_indexes = 0
		low_quality_indexes = 0
		total_read_count = 0
		with gzip.open(inp_file, 'rt') as input_file:
			for record in SeqIO.parse(input_file, 'fastq'):
				record_sequence = str(record.seq)
				total_read_count += 1
				if total_read_count % 1000000 == 0:
					print('%s reads were analyzed.' % str(total_read_count))
					print(datetime.now())
				for ind in ind_dictionary[inp_file]:
					if ind == 'R701':
						score, index_length, min_score = aligner(
								record.seq, param.index_dict['R701'],
								param.index_mismatch_R701
								)
					else:
						score, index_length, min_score = aligner(
								record.seq, param.index_dict[ind],
								param.index_mismatch_8bp
								)
					if score >= min_score:
						index = ind
						break
				else:
					undetermined_indexes += 1
					continue
				if check_quality(record, 0, index_length):
					ind_dictionary[inp_file][index]['read count'] += 1
					if index in param.norm_or_expr_indexes:
						barcode = find_barcode(
								record, index_length,
								param.neCP1,
								param.neCP1_mismatch,
								param.neCP2[:param.neCP2_length],
								param.neCP2_mismatch, param.ne_mode
								)
						if type(barcode) is not str and param.ne_mode == 'reverse':
							barcode = str(barcode.reverse_complement())
						if type(barcode) is not str and param.ne_mode == 'forward':
							barcode = str(barcode)	
						if barcode not in ind_dictionary[inp_file][index]:
							ind_dictionary[inp_file][index][barcode] = [record_sequence]
						else:
							ind_dictionary[inp_file][index][barcode].append(record_sequence)
					else:
						barcode, ROI = find_barcode_ROI(
								record, index, index_length,
								lib_dictionary, inp_file, param.map_mode
								)
						if type(barcode) is not str and param.map_mode == 'forward':
							barcode = str(barcode)
							ROI = str(ROI)
						if type(barcode) is not str and param.map_mode == 'reverse':
							barcode = str(barcode.reverse_complement())
							ROI = str(ROI.reverse_complement())
						if barcode not in ind_dictionary[inp_file][index]:
							ind_dictionary[inp_file][index][barcode] = {ROI: [record_sequence]}
						else:
							if ROI not in ind_dictionary[inp_file][index][barcode]:
								ind_dictionary[inp_file][index][barcode][ROI] = [record_sequence]
							else:
								ind_dictionary[inp_file][index][barcode][ROI].append(record_sequence)
				else:
					low_quality_indexes += 1
		
		lines_of_text = [
				'Total read count in file %s is %s.\n' %(inp_file, str(total_read_count)),
				'In file ' + inp_file + ' was found:\n'
				]
		for index in ind_dictionary[inp_file]:
			if index in param.mapping_indexes:
				try:
					ind_dictionary[inp_file][index]['low quality count'] = len(
							ind_dictionary[inp_file][index]['low quality barcode']['low quality barcode']
							) + len(ind_dictionary[inp_file][index]['low quality ROI']['low quality ROI'])
					delete = 6
				except KeyError:
					ind_dictionary[inp_file][index]['low quality count'] = len(
							ind_dictionary[inp_file][index]['low quality']['undetermined']
							) + len(
									ind_dictionary[inp_file][index]['bc15']['low quality']
									) + len(ind_dictionary[inp_file][index]['bc21']['low quality'])
					delete = 5
				ind_dictionary[inp_file][index]['eff read count'] = ind_dictionary[inp_file][index]['read count'] - len(ind_dictionary[inp_file][index]['undetermined']['undetermined']) - ind_dictionary[inp_file][index]['low quality count']
			else:
				ind_dictionary[inp_file][index]['low quality count'] = len(
						ind_dictionary[inp_file][index]['low quality']
						)
				ind_dictionary[inp_file][index]['eff read count'] = ind_dictionary[inp_file][index]['read count'] - len(ind_dictionary[inp_file][index]['undetermined']) - ind_dictionary[inp_file][index]['low quality count']
				delete = 5
			lines_of_text.append('%s total reads, %s effective reads, %s low quality reads and %s unique barcodes for index %s.\n'
				  %(str(ind_dictionary[inp_file][index]['read count']),
					str(ind_dictionary[inp_file][index]['eff read count']),
					str(ind_dictionary[inp_file][index]['low quality count']),
					str(len(ind_dictionary[inp_file][index]) - delete), index))
		lines_of_text += [
				str(undetermined_indexes) + ' reads without any index.\n',
				str(low_quality_indexes) + ' reads with low quality index.\n'
				]
		
		print(*lines_of_text)
		
		with open(run_info_text, 'a') as info_file:
			info_file.writelines(lines_of_text)
	
	print('Analysis was finished ' + str(datetime.now()))
	with open(output_folder + '/' + 'fastq_file_dictionary.json', 'w') as f:
		json.dump(ind_dictionary, f)
	return ind_dictionary

#no demultiplexing
def fastq_parsing1(ind_dictionary, lib_dictionary):
	start_time = datetime.now()
	print('Analysis was started ' + str(start_time))
	for inp_file in ind_dictionary:
		total_read_count = 0
		index_length = 0
		undetermined_indexes = 0
		low_quality_index = 0
		with gzip.open(inp_file, 'rt') as input_file:
			for record in SeqIO.parse(input_file, 'fastq'):
				record_sequence = str(record.seq)
				total_read_count += 1
				if total_read_count % 1000000 == 0:
					print('%s reads were analyzed.' % str(total_read_count))
					print(datetime.now())
				for ind in ind_dictionary[inp_file]:
					index = ind
				if check_quality(record, 0, 0):
					ind_dictionary[inp_file][index]['read count'] += 1
					if index in param.norm_or_expr_indexes:
						barcode = find_barcode(
								record, index_length,
								param.neCP1,
								param.neCP1_mismatch,
								param.neCP2[:param.neCP2_length],
								param.neCP2_mismatch, param.ne_mode
								)
						if type(barcode) is not str and param.ne_mode == 'reverse':
							barcode = str(barcode.reverse_complement())
						if type(barcode) is not str and param.ne_mode == 'forward':
							barcode = str(barcode)
						if barcode not in ind_dictionary[inp_file][index]:
							ind_dictionary[inp_file][index][barcode] = [record_sequence]
						else:
							ind_dictionary[inp_file][index][barcode].append(record_sequence)
					else:
						barcode, ROI = find_barcode_ROI(
								record, index, index_length,
								lib_dictionary, inp_file, param.map_mode
								)
						if type(barcode) is not str and param.map_mode == 'forward':
							barcode = str(barcode)
							ROI = str(ROI)
						if type(barcode) is not str and param.map_mode == 'reverse':
							barcode = str(barcode.reverse_complement())
							ROI = str(ROI.reverse_complement())
						if barcode not in ind_dictionary[inp_file][index]:
							ind_dictionary[inp_file][index][barcode] = {ROI: [record_sequence]}
						else:
							if ROI not in ind_dictionary[inp_file][index][barcode]:
								ind_dictionary[inp_file][index][barcode][ROI] = [record_sequence]
							else:
								ind_dictionary[inp_file][index][barcode][ROI].append(record_sequence)
				else:
					low_quality_index += 1
		
		lines_of_text = [
				'Total read count in file %s is %s.\n' %(inp_file, str(total_read_count)),
				'In file ' + inp_file + ' was found:\n'
				]
		for index in ind_dictionary[inp_file]:
			if index in param.mapping_indexes:
				try:
					ind_dictionary[inp_file][index]['low quality count'] = len(
							ind_dictionary[inp_file][index]['low quality barcode']['low quality barcode']
							) + len(ind_dictionary[inp_file][index]['low quality ROI']['low quality ROI'])
					delete = 6
				except KeyError:
					ind_dictionary[inp_file][index]['low quality count'] = len(
							ind_dictionary[inp_file][index]['low quality']['undetermined']
							) + len(
									ind_dictionary[inp_file][index]['bc15']['low quality']
									) + len(ind_dictionary[inp_file][index]['bc21']['low quality'])
					delete = 5
				ind_dictionary[inp_file][index]['eff read count'] = ind_dictionary[inp_file][index]['read count'] - len(ind_dictionary[inp_file][index]['undetermined']['undetermined']) - ind_dictionary[inp_file][index]['low quality count']
			else:
				ind_dictionary[inp_file][index]['low quality count'] = len(
						ind_dictionary[inp_file][index]['low quality']
						)
				ind_dictionary[inp_file][index]['eff read count'] = ind_dictionary[inp_file][index]['read count'] - len(ind_dictionary[inp_file][index]['undetermined']) - ind_dictionary[inp_file][index]['low quality count']
				delete = 5
			lines_of_text.append('%s total reads, %s effective reads, %s low quality reads and %s unique barcodes for index %s.\n'
				  %(str(ind_dictionary[inp_file][index]['read count']),
					str(ind_dictionary[inp_file][index]['eff read count']),
					str(ind_dictionary[inp_file][index]['low quality count']),
					str(len(ind_dictionary[inp_file][index]) - delete), index))
				
		print(*lines_of_text)
		
		with open(run_info_text, 'a') as info_file:
			info_file.writelines(lines_of_text)
	
	print('Analysis was finished ' + str(datetime.now()))
	with open(output_folder + '/' + 'fastq_file_dictionary.json', 'w') as f:
		json.dump(ind_dictionary, f)
	return ind_dictionary			
			
def demultiplexing_or_not():
	ans = input('Do you need to perform demultiplexing of your fastq file(s)? Type Y or N: ').upper()
	if ans == 'Y' or ans == 'N':
		return ans
	else:
		print('You entered an incorrect value: %s. Try again.' % ans)
		return demultiplexing_or_not()

reply = demultiplexing_or_not()

if reply == 'Y':
	fastq_parsing(fastq_file_dictionary, lib_diction)
if reply == 'N':
	fastq_parsing1(fastq_file_dictionary, lib_diction)
		

#print(lib_diction)
#n = 0
#for read in fastq_file_dictionary['C:/Users/user/1_S1_L001_R1_001.fastq.gz']['R701']['bc15']['undetermined']:
#	print(read)
#	n += 1
#	if n > 50:
#		break
#for barcode in fastq_file_dictionary['C:\\Users\\user\\1_S1_L001_R1_001_2018_10_05.fastq.gz']['A1']:
#	print(barcode)
#print(fastq_file_dictionary['C:\\Users\\user\\1_S1_L001_R1_001_2018_11_20.fastq.gz']['A1']['CCCGGAGAGGGGGGGAAG']['CATCTTAT'])

#print(len(fastq_file_dictionary['C:/Users/user/1_S1_L001_R1_001.fastq.gz']['A4']['bc21']['undetermined']))


#Form dictionary of barcodes-ROIs, found in at least two replicates.
#gold_mapping = {
#				'barcode1': {'mut1': mut1_count, 'mut2': mut2_count},
#				'barcode2': {'mut3': mut3_count, 'mut4': mut4_count}
#			   }
#mapping_count = {'file name': {'A1': {
#									  'barcode1': {
#												   'mut1': mut1_count1,
#												   'mut2': mut2_count1
#												  },
#									  'barcode2': {
#												   'mut3': mut3_count1,
#												   'mut4': mut4_count1
#												   }
#									 },
#							   'A2': {
#									  'barcode1': {
#												   'mut1': mut1_count2,
#												   'mut2': mut2_count2
#												   },
#									  'barcode2': {
#												   'mut3': mut3_count2,
#												   'mut4': mut4_count2
#												  }
#									 }}}
def form_gold_mapping(ind_dictionary, lib_dictionary, lib):
	gold_mapping = {}
	#This includes barcodes that are not found in other replicates.
	single_mapping = {}
	#This is temporal dictionary. Includes reads, found only in this replicate.
	rep_mapping = {}
	#This dictionary includes all barcodes, ROIs and ROI counts for each replicate.
	mapping_count = {}
	reference_BCs_dict = {}
	if len(lib_dictionary[lib]['mapping']) == 1:
		for file in lib_dictionary[lib]['mapping']:
			mapping_count[file] = {}
			reference_BCs_dict[file] = {}
			if len(lib_dictionary[lib]['mapping'][file]) == 1:
				for index_list in lib_dictionary[lib]['mapping'][file]:
					reference_BCs_dict[file][', '.join(index_list)] = {}
					for index in index_list:
						for barcode in ind_dictionary[file][index]:
							#Check that it is not 'read count' or control barcode.
							if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode not in [
									'low quality count',
									'low quality barcode',
									'low quality ROI'
									]:
								for ROI in ind_dictionary[file][index][barcode]:
									if len(ROI) >= param.min_ROI_length and len(ROI) <= param.max_ROI_length:
										ROI_count = len(
												ind_dictionary[file][index][barcode][ROI]
												) * 1000000 / ind_dictionary[file][index]['eff read count']
										if barcode in gold_mapping:
											if ROI in gold_mapping[barcode]:
												gold_mapping[barcode][ROI] += ROI_count
											else:
												gold_mapping[barcode][ROI] = ROI_count
										else:
											gold_mapping[barcode] = {
													ROI: ROI_count
													}
							elif barcode in param.reference_BCs['forward']:
								if barcode == 'control1' or barcode == 'control2':
									ROI = 'control'
								else:
									ROI = 'experiment'
								if barcode in reference_BCs_dict[file][', '.join(index_list)]:
									reference_BCs_dict[file][', '.join(index_list)][barcode] += len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
								else:
									reference_BCs_dict[file][', '.join(index_list)][barcode] = len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
					mapping_count[file][', '.join(index_list)] = gold_mapping
			else:	
				for i in range(len(lib_dictionary[lib]['mapping'][file])):
					rep_name = ', '.join(
							lib_dictionary[lib]['mapping'][file][i]
							)
					mapping_count[file][rep_name] = {}
					reference_BCs_dict[file][rep_name] = {}
					for index in lib_dictionary[lib]['mapping'][file][i]:
						for barcode in ind_dictionary[file][index]:
							#Check that it is not 'read count' or control barcode.
							if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode not in [
									'low quality count',
									'low quality barcode',
									'low quality ROI'
									]:
								for ROI in ind_dictionary[file][index][barcode]:
									if len(ROI) >= param.min_ROI_length and len(ROI) <= param.max_ROI_length:
										ROI_count = len(
												ind_dictionary[file][index][barcode][ROI]
												) * 1000000 / ind_dictionary[file][index]['eff read count']
										if barcode in mapping_count[file][rep_name]:
											if ROI in mapping_count[file][rep_name][barcode]:
												mapping_count[file][rep_name][barcode][ROI] += ROI_count
											else:
												mapping_count[file][rep_name][barcode][ROI] = ROI_count
										else:
											mapping_count[file][rep_name][barcode] = {
													ROI: ROI_count
													}
										if barcode in gold_mapping:
											if ROI in gold_mapping[barcode]:
												gold_mapping[barcode][ROI] += ROI_count
											else:
												for n in range(i + 1, len(
														lib_dictionary[lib]['mapping'][file]
														)):
													for ind in lib_dictionary[lib]['mapping'][file][n]:
														if barcode in ind_dictionary[file][ind]:
															if ROI in ind_dictionary[file][ind][barcode]:
																gold_mapping[barcode][ROI] = ROI_count
																break
													else:
														continue
													break
										else:
											for n in range(i + 1, len(
													lib_dictionary[lib]['mapping'][file]
													)):
												for ind in lib_dictionary[lib]['mapping'][file][n]:
													if barcode in ind_dictionary[file][ind]:
														if ROI in ind_dictionary[file][ind][barcode]:
															gold_mapping[barcode] = {
																	ROI: ROI_count
																	}
															break
												else:
													continue
												break
							elif barcode in param.reference_BCs['forward']:
								if barcode == 'control1' or barcode == 'control2':
									ROI = 'control'
								else:
									ROI = 'experiment'
								if barcode in reference_BCs_dict[file][rep_name]:
									reference_BCs_dict[file][rep_name][barcode] += len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
								else:
									reference_BCs_dict[file][rep_name][barcode] = len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
	else:
		for file in lib_dictionary[lib]['mapping']:
			mapping_count[file] = {}
			reference_BCs_dict[file] = {}
			if len(lib_dictionary[lib]['mapping'][file]) == 1:
				for index_list in lib_dictionary[lib]['mapping'][file]:
#This is temporal dictionary. Includes reads, found only in this replicate.
					rep_mapping = {}
					rep_name = ', '.join(index_list)
					mapping_count[file][rep_name] = {}
					reference_BCs_dict[file][rep_name] = {}
					for index in index_list:
						for barcode in ind_dictionary[file][index]:
							#Check that it is not 'read count' or control barcode.
							if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode not in [
									'low quality count',
									'low quality barcode',
									'low quality ROI'
									]:
								for ROI in ind_dictionary[file][index][barcode]:
									if len(ROI) >= param.min_ROI_length and len(ROI) <= param.max_ROI_length:
										ROI_count = len(
													ind_dictionary[file][index][barcode][ROI]
													) * 1000000 / ind_dictionary[file][index]['eff read count']
										if barcode in mapping_count[file][rep_name]:
											if ROI in mapping_count[file][rep_name][barcode]:
												mapping_count[file][rep_name][barcode][ROI] += ROI_count
											else:
												mapping_count[file][rep_name][barcode][ROI] = ROI_count
										else:
											mapping_count[file][rep_name][barcode] = {
													ROI: ROI_count
													}
										if barcode in gold_mapping:
											if ROI in gold_mapping[barcode]:
												gold_mapping[barcode][ROI] += ROI_count
											else:
												if barcode in single_mapping:
													if ROI in single_mapping[barcode]:
														gold_mapping[barcode][ROI] = ROI_count
														single_mapping[barcode].pop(
																ROI
																)
													else:
														if barcode in rep_mapping:
															if ROI in rep_mapping:
																rep_mapping[barcode][ROI] += ROI_count
															else:
																rep_mapping[barcode][ROI] = ROI_count
														else:
															rep_mapping[barcode] = {
																	ROI: ROI_count
																	}
												else:
													if barcode in rep_mapping:
														if ROI in rep_mapping:
															rep_mapping[barcode][ROI] += ROI_count
														else:
															rep_mapping[barcode][ROI] = ROI_count
													else:
														rep_mapping[barcode] = {
																ROI: ROI_count
																}
										else:
											if barcode in single_mapping:
												if ROI in single_mapping[barcode]:
													gold_mapping[barcode] = {
															ROI: ROI_count
															}
													single_mapping[barcode].pop(
															ROI
															)
												else:
													if barcode in rep_mapping:
														if ROI in rep_mapping:
															rep_mapping[barcode][ROI] += ROI_count
														else:
															rep_mapping[barcode][ROI] = ROI_count
													else:
														rep_mapping[barcode] = {
																ROI: ROI_count
																}
											else:
												if barcode in rep_mapping:
													if ROI in rep_mapping:
														rep_mapping[barcode][ROI] += ROI_count
													else:
														rep_mapping[barcode][ROI] = ROI_count
												else:
													rep_mapping[barcode] = {
															ROI: ROI_count
															}
							elif barcode in param.reference_BCs['forward']:
								if barcode == 'control1' or barcode == 'control2':
									ROI = 'control'
								else:
									ROI = 'experiment'
								if barcode in reference_BCs_dict[file][rep_name]:
									reference_BCs_dict[file][rep_name][barcode] += len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
								else:
									reference_BCs_dict[file][rep_name][barcode] = len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
					for barcode in rep_mapping:
						for ROI in rep_mapping[barcode]:
							ROI_count = rep_mapping[barcode][ROI]
							if barcode in single_mapping:
								single_mapping[barcode][ROI] = ROI_count
							else:
								single_mapping[barcode] = {
										ROI: ROI_count
										}
			else:
				for i in range(len(lib_dictionary[lib]['mapping'][file])):
#This is temporal dictionary. Includes reads, found only in this replicate.
					rep_mapping = {}
					rep_name = ', '.join(
							lib_dictionary[lib]['mapping'][file][i]
							)
					mapping_count[file][rep_name] = {}
					reference_BCs_dict[file][rep_name] = {}
					for index in lib_dictionary[lib]['mapping'][file][i]:
						for barcode in ind_dictionary[file][index]:
							#Check that it is not 'read count' or control barcode.
							if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode not in [
									'low quality count',
									'low quality barcode',
									'low quality ROI'
									]:
								for ROI in ind_dictionary[file][index][barcode]:
									if len(ROI) >= param.min_ROI_length and len(ROI) <= param.max_ROI_length:
										ROI_count = len(
												ind_dictionary[file][index][barcode][ROI]
												) * 1000000 / ind_dictionary[file][index]['eff read count']
										if barcode in mapping_count[file][rep_name]:
											if ROI in mapping_count[file][rep_name][barcode]:
												mapping_count[file][rep_name][barcode][ROI] += ROI_count
											else:
												mapping_count[file][rep_name][barcode][ROI] = ROI_count
										else:
											mapping_count[file][rep_name][barcode] = {
													ROI: ROI_count
													}
										if barcode in gold_mapping:
											if ROI in gold_mapping[barcode]:
												gold_mapping[barcode][ROI] += ROI_count
											else:
												for n in range(i + 1, len(
														lib_dictionary[lib]['mapping'][file]
														)):
													for ind in lib_dictionary[lib]['mapping'][file][n]:
														if barcode in ind_dictionary[file][ind]:
															if ROI in ind_dictionary[file][ind][barcode]:
																gold_mapping[barcode][ROI] = ROI_count
																break
													else:
														continue
													break
												else:
													if barcode in rep_mapping:
														if ROI in rep_mapping:
															rep_mapping[barcode][ROI] += ROI_count
														else:
															rep_mapping[barcode][ROI] = ROI_count
													else:
														rep_mapping[barcode] = {
																ROI: ROI_count
																}
										else:
											for n in range(i + 1, len(
													lib_dictionary[lib]['mapping'][file]
													)):
												for ind in lib_dictionary[lib]['mapping'][file][n]:
													if barcode in ind_dictionary[file][ind]:
														if ROI in ind_dictionary[file][ind][barcode]:
															gold_mapping[barcode] = {
																	ROI: ROI_count
																	}
															break
												else:
													continue
												break
											else:
												if barcode in rep_mapping:
													if ROI in rep_mapping:
														rep_mapping[barcode][ROI] += ROI_count
													else:
														rep_mapping[barcode][ROI] = ROI_count
												else:
													rep_mapping[barcode] = {
															ROI: ROI_count
															}
							elif barcode in param.reference_BCs['forward']:
								if barcode == 'control1' or barcode == 'control2':
									ROI = 'control'
								else:
									ROI = 'experiment'
								if barcode in reference_BCs_dict[file][rep_name]:
									reference_BCs_dict[file][rep_name][barcode] += len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
								else:
									reference_BCs_dict[file][rep_name][barcode] = len(
											ind_dictionary[file][index][barcode][ROI]
											) * 1000000 / ind_dictionary[file][index]['eff read count']
					for barcode in rep_mapping:
						for ROI in rep_mapping[barcode]:
							ROI_count = rep_mapping[barcode][ROI]
							if barcode in single_mapping:
								if ROI in single_mapping:
									if barcode in gold_mapping:
										gold_mapping[barcode][ROI] = ROI_count
									else:
										gold_mapping[barcode] = {
												ROI: ROI_count
												}
								else:
									single_mapping[barcode][ROI] = ROI_count
							else:
								single_mapping[barcode] = {
										ROI: ROI_count
										}
	if not gold_mapping:
		raise Exception('Error! There was not found any valid BC-ROI pair for mapping of library %s'
						% lib)

	return gold_mapping, mapping_count, reference_BCs_dict


#Form dictionary of barcodes, found in all replicates and
#dictionary of read count per million effective reads for each replicate.
#gold_barcode_list = {
#					 'barcode1': (barcode1_count1 + barcode1_count2),
#					 'barcode2': (barcode2_count1 + barcode2_count2)
#					}
#barcode_count_dict = {'file name': {'A10': {
#											'barcode1': barcode1_count1,
#											'barcode2': barcode2_count1},
#									'A11': {
#											'barcode1': barcode1_count2,
#											'barcode2': barcode2_count2
#										   }}}
def form_gold_norm_or_expr(ind_dictionary, lib_dictionary, lib, rep_type):
	gold_barcode_dict = {}
	#Includes all barcodes of first replicate.
	first_barcode_dict = {}
	first_file_dict = {}
	#This dictionary includes all BCs and their counts for each replicate.
	barcode_count_dict = {}
	reference_BCs_dict = {}
	if len(lib_dictionary[lib][rep_type]) == 1:
		for file in lib_dictionary[lib][rep_type]:
			barcode_count_dict[file] = {}
			reference_BCs_dict[file] = {}
			if len(lib_dictionary[lib][rep_type][file]) == 1:
				for index_list in lib_dictionary[lib][rep_type][file]:
					reference_BCs_dict[file][', '.join(index_list)] = {}
					rep_eff_read_count = 0
					for index in index_list:
						rep_eff_read_count += ind_dictionary[file][index]['eff read count']
						for barcode in ind_dictionary[file][index]:
							#Check that it is not 'read count' or reference BC.
							if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode != 'low quality count':
								barcode_count = len(
										ind_dictionary[file][index][barcode]
										)
								if barcode in gold_barcode_dict:
									first_barcode_dict[barcode] += barcode_count
								else:
									first_barcode_dict[barcode] = barcode_count
							elif barcode in param.reference_BCs['forward']:
								if barcode in reference_BCs_dict[file][', '.join(index_list)]:
									reference_BCs_dict[file][', '.join(index_list)][barcode] += len(
											ind_dictionary[file][index][barcode]
											)
								else:
									reference_BCs_dict[file][', '.join(index_list)][barcode] = len(
											ind_dictionary[file][index][barcode]
											)
					for control_barcode in reference_BCs_dict[file][', '.join(index_list)]:
						reference_BCs_dict[file][', '.join(index_list)][control_barcode] = reference_BCs_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
				for barcode in first_barcode_dict:
					#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
					if first_barcode_dict[barcode] >= param.minimum_read_count[rep_type]:
						gold_barcode_dict[barcode] = first_barcode_dict[barcode] / rep_eff_read_count * 1000000
				barcode_count_dict[file][', '.join(index_list)] = gold_barcode_dict
			else:
				first_rep_eff_read_count = 0
				reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])] = {}
				for index in lib_dictionary[lib][rep_type][file][0]:
					first_rep_eff_read_count += ind_dictionary[file][index]['eff read count']
					for barcode in ind_dictionary[file][index]:
						#Check that it is not 'read count' or reference BC.
						if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode != 'low quality count':
							barcode_count = len(
									ind_dictionary[file][index][barcode]
									)
							if barcode in first_barcode_dict:
								first_barcode_dict[barcode] += barcode_count
							else:
								first_barcode_dict[barcode] = barcode_count
						elif barcode in param.reference_BCs['forward']:
							if barcode in reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])]:
								reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][barcode] += len(
										ind_dictionary[file][index][barcode]
										)
							else:
								reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][barcode] = len(
										ind_dictionary[file][index][barcode]
										)
				for control_barcode in reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])]:
					reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][control_barcode] = reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][control_barcode] / first_rep_eff_read_count * 1000000
				for barcode in first_barcode_dict:
					first_barcode_dict[barcode] = first_barcode_dict[barcode] / first_rep_eff_read_count * 1000000
				barcode_count_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])] = first_barcode_dict.copy()
				
				for n in range(1, len(lib_dictionary[lib][rep_type][file])):
					barcode_count_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])] = {}
					reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])] = {}
					rep_eff_read_count = 0
					for ind in lib_dictionary[lib][rep_type][file][n]:
						rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
						for control_barcode in param.reference_BCs['forward']:
							if len(param.reference_BCs['forward']) > 0:
								if control_barcode in ind_dictionary[file][ind]:
									if control_barcode in reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])]:
										reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] += len(
												ind_dictionary[file][ind][control_barcode]
												)
									else:
										reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] = len(
												ind_dictionary[file][ind][control_barcode]
												)
					for control_barcode in reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])]:
						reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] = reference_BCs_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] / rep_eff_read_count * 1000000
				
				for barcode in first_barcode_dict:
					#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
					if round(first_barcode_dict[barcode] / 1000000 * first_rep_eff_read_count) >= param.minimum_read_count[rep_type]:
						for n in range(1, len(lib_dictionary[lib][rep_type][file])):
							rep_eff_read_count = 0
							rep_barcode_count = 0
							for ind in lib_dictionary[lib][rep_type][file][n]:
								rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
								if barcode in ind_dictionary[file][ind]:
									rep_barcode_count += len(
											ind_dictionary[file][ind][barcode]
											)
							#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
							if rep_barcode_count < param.minimum_read_count[rep_type]:
								break
							else:
								rep_barcode_count = rep_barcode_count / rep_eff_read_count * 1000000
								first_barcode_dict[barcode] += rep_barcode_count
								barcode_count_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][barcode] = rep_barcode_count
						else:
							gold_barcode_dict[barcode] = first_barcode_dict[barcode]
	else:
		for file in lib_dictionary[lib][rep_type]:
			barcode_count_dict[file] = {}
			reference_BCs_dict[file] = {}
		first_rep_eff_read_count = 0
		for frst_file in lib_dictionary[lib][rep_type]:
			first_file = frst_file
			break
		
		reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])] = {}
		for index in lib_dictionary[lib][rep_type][first_file][0]:
			first_rep_eff_read_count += ind_dictionary[first_file][index]['eff read count']
			for barcode in ind_dictionary[first_file][index]:
			#Check that it is not 'read count' or reference barcode.
				if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode != 'low quality count':
					barcode_count = len(
							ind_dictionary[first_file][index][barcode]
							)
					if barcode in first_barcode_dict:
						first_barcode_dict[barcode] += barcode_count
					else:
						first_barcode_dict[barcode] = barcode_count
				elif barcode in param.reference_BCs['forward']:
					if barcode in reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])]:
						reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][barcode] += len(
								ind_dictionary[first_file][index][barcode]
								)
					else:
						reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][barcode] = len(
								ind_dictionary[first_file][index][barcode]
								)
		for control_barcode in reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])]:
			reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][control_barcode] = reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][control_barcode] / first_rep_eff_read_count * 1000000
		
		for barcode in first_barcode_dict:
			first_barcode_dict[barcode] = first_barcode_dict[barcode] / first_rep_eff_read_count * 1000000
		barcode_count_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])] = first_barcode_dict.copy()
			
		for n in range(1, len(lib_dictionary[lib][rep_type][first_file])):
			barcode_count_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = {}
			reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = {}
			rep_eff_read_count = 0
			for ind in lib_dictionary[lib][rep_type][first_file][n]:
				rep_eff_read_count += ind_dictionary[first_file][ind]['eff read count']
				for control_barcode in param.reference_BCs['forward']:
					if len(param.reference_BCs['forward']) > 0:	
						if control_barcode in ind_dictionary[first_file][ind]:
							if control_barcode in reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])]:
								reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] += len(
										ind_dictionary[first_file][ind][control_barcode]
										)
			for control_barcode in reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])]:
				reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = reference_BCs_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] / rep_eff_read_count * 1000000
		
		for barcode in first_barcode_dict:
			#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
			if round(first_barcode_dict[barcode] * first_rep_eff_read_count / 1000000) >= param.minimum_read_count[rep_type]:
				for n in range(1, len(lib_dictionary[lib][rep_type][first_file])):
					rep_eff_read_count = 0
					rep_barcode_count = 0
					for ind in lib_dictionary[lib][rep_type][first_file][n]:
						rep_eff_read_count += ind_dictionary[first_file][ind]['eff read count']
						if barcode in ind_dictionary[first_file][ind]:
							rep_barcode_count += len(
									ind_dictionary[first_file][ind][barcode]
									)
					#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
					if rep_barcode_count < param.minimum_read_count[rep_type]:
						break
					else:
						rep_barcode_count = rep_barcode_count / rep_eff_read_count * 1000000
						first_barcode_dict[barcode] += rep_barcode_count
						barcode_count_dict[lib_dictionary[lib][rep_type][0]][', '.join(lib_dictionary[lib][rep_type][0][n])][barcode] = rep_barcode_count
				else:
					first_file_dict[barcode] = first_barcode_dict[barcode]
			
		for file in lib_dictionary[lib][rep_type]:
			if file != first_file:
				for index_list in lib_dictionary[lib][rep_type][file]:
					barcode_count_dict[file][', '.join(index_list)] = {}
					reference_BCs_dict[file][', '.join(index_list)] = {}
					rep_eff_read_count = 0
					for ind in index_list:
						rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
						for control_barcode in param.reference_BCs['forward']:
							if len(param.reference_BCs['forward']) > 0:	
								if control_barcode in ind_dictionary[file][ind]:
									if control_barcode in reference_BCs_dict[file][', '.join(index_list)]:
										reference_BCs_dict[file][', '.join(index_list)][control_barcode] += len(
												ind_dictionary[file][ind][control_barcode]
												)
									else:
										reference_BCs_dict[file][', '.join(index_list)][control_barcode] = len(
												ind_dictionary[file][ind][control_barcode]
												)
					for control_barcode in reference_BCs_dict[file][', '.join(index_list)]:
						reference_BCs_dict[file][', '.join(index_list)][control_barcode] = reference_BCs_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
						
		for barcode in first_file_dict:
			for file in lib_dictionary[lib][rep_type]:
				if file != first_file:
					for index_list in lib_dictionary[lib][rep_type][file]:
						rep_eff_read_count = 0
						rep_barcode_count = 0
						for ind in index_list:
							rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
							if barcode in ind_dictionary[file][ind]:
								rep_barcode_count += len(
										ind_dictionary[file][ind][barcode]
										)
						#Check that number of reads and number of reads per million effective reads are both greater than minimum_read_count.
						if rep_barcode_count < param.minimum_read_count[rep_type]:
							break
						else:
							rep_barcode_count = rep_barcode_count / rep_eff_read_count * 1000000
							first_file_dict[barcode] += rep_barcode_count
							barcode_count_dict[file][', '.join(index_list)][barcode] = rep_barcode_count
					else:
						continue
					break
			else:
				gold_barcode_dict[barcode] = first_file_dict[barcode]

	if not gold_barcode_dict:
		raise Exception('Error! BCs were not found for normalization of library %s'
						% lib)
	
	return gold_barcode_dict, barcode_count_dict, reference_BCs_dict


#This function works only if minimum read count for expression = 0.
def form_rep_expr(ind_dictionary, lib_dictionary, lib):
	barcode_count_dict = {}
	reference_BCs_dict = {}
	total_expression = {}
	total_eff_read_count = 0
	for file in lib_dictionary[lib]['expression']:
		barcode_count_dict[file] = {}
		reference_BCs_dict[file] = {}
		for index_list in lib_dictionary[lib]['expression'][file]:
			barcode_count_dict[file][', '.join(index_list)] = {}
			reference_BCs_dict[file][', '.join(index_list)] = {}
			rep_eff_read_count = 0
			for index in index_list:
				rep_eff_read_count += ind_dictionary[file][index]['eff read count']
				total_eff_read_count += ind_dictionary[file][index]['eff read count']
				for barcode in ind_dictionary[file][index]:
					#Check that it is not 'read count' or reference BC.
					if len(barcode) >= param.min_BC_length and len(barcode) <= param.max_BC_length and barcode != 'low quality count':
						if barcode in barcode_count_dict[file][', '.join(index_list)]:
							barcode_count_dict[file][', '.join(index_list)][barcode] += len(
									ind_dictionary[file][index][barcode]
									)
						else:
							barcode_count_dict[file][', '.join(index_list)][barcode] = len(
									ind_dictionary[file][index][barcode]
									)
						if barcode in total_expression:
							total_expression[barcode] += len(
									ind_dictionary[file][index][barcode]
									)
						else:
							total_expression[barcode] = len(
									ind_dictionary[file][index][barcode]
									)
					elif barcode in param.reference_BCs['forward']:
						if barcode in reference_BCs_dict[file][', '.join(index_list)]:
							reference_BCs_dict[file][', '.join(index_list)][barcode] += len(
									ind_dictionary[file][index][barcode]
									)
						else:
							reference_BCs_dict[file][', '.join(index_list)][barcode] = len(
									ind_dictionary[file][index][barcode]
									)
				for barcode in barcode_count_dict[file][', '.join(index_list)]:
					barcode_count_dict[file][', '.join(index_list)][barcode] = barcode_count_dict[file][', '.join(index_list)][barcode] / rep_eff_read_count * 1000000
				for control_barcode in reference_BCs_dict[file][', '.join(index_list)]:
					reference_BCs_dict[file][', '.join(index_list)][control_barcode] = reference_BCs_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
	for barcode in total_expression:
		total_expression[barcode] = total_expression[barcode] / total_eff_read_count * 1000000
	return barcode_count_dict, reference_BCs_dict, total_expression


#Return True if query is a mutant variant of ref.
def align_mut_variants(ref, query, mode):
	if mode == 'barcode':
		mismatch = param.BC_mismatch
		length = param.BC_length
	else:
		mismatch = param.ROI_mismatch
		length = param.ROI_length
	align_score = pairwise2.align.globalms(ref, query, 1, -1, -1, -1, score_only=True)
	if align_score >= length - mismatch:
		return True
	return False


def not_empty(barc_dict, rep_types, library):
	if not barc_dict:
		raise Exception('Error! BCs were not found in %s of library %s.'
						%(rep_types, library))


#Form gold dictionary of unique barcodes
#(only works if minimum_read_count > 0).
#gold_dictionary = {
#				 '29-36': {
#						   'barcodes': {'barcode1': barcode1_count1, 'barcode2': barcode1_count2},
#						   'barcodes-ROIs': {
#												   'barcode1': {'mut1': mut1_count, 'mut2': mut2_count},
#												   'barcode2': {'mut3': mut3_count, 'mut4': mut4_count}
#												   }
#						  }
#				 '25-32': {
#						   'barcodes': {'barcode3': barcode1_count3, 'barcode4': barcode1_count4},
#						   'barcodes-ROIs': {}
#						   }
#				   }
#rep_dictionary = {'25-32':{
#				  'normalization': {'file name': {
#									'A10': {
#											'barcode1': barcode1_count1,
#											'barcode2': barcode2_count1},
#									'A11': {
#											'barcode1': barcode1_count2,
#											'barcode2': barcode2_count2
#										   }
#												  }},
#				   'expression': {},
#				   'mapping': {'file name': {
#								 'A1': {
#									  'barcode1': {
#												   'mut1': mut1_count1,
#												   'mut2': mut2_count1
#												  },
#									  'barcode2': {
#												   'mut3': mut3_count1,
#												   'mut4': mut4_count1
#												   }
#									 },
#							   'A2': {
#									  'barcode1': {
#												   'mut1': mut1_count2,
#												   'mut2': mut2_count2
#												   },
#									  'barcode2': {
#												   'mut3': mut3_count2,
#												   'mut4': mut4_count2
#												  }
#									 }
#										  }
#							  }}}
def form_gold_dictionary(ind_dictionary, lib_dictionary):
	gold_dictionary = {}
	rep_dictionary = {}
	cont_dictionary = {}
	total_rep_dictionary = {}
	for lib in lib_dictionary:
		if lib != 'H':
			gold_dictionary[lib] = {'barcodes': {}, 'barcodes-ROIs': {}}
			gold_mapping = {}
			gold_normalization = {}
			gold_expression = {}
			rep_dictionary[lib] = {
					'normalization': {}, 'expression': {}, 'mapping': {}
					}
			cont_dictionary[lib] = {
					'normalization': {}, 'expression': {}, 'mapping': {}
					}
			lines_of_text = []
			total_expr = {}
			total_rep_dictionary[lib] = {}
			
			if lib_dictionary[lib]['mapping']:
				gold_mapping, rep_dictionary[lib]['mapping'], cont_dictionary[lib]['mapping'] = form_gold_mapping(
						ind_dictionary, lib_dictionary, lib
						)
				lines_of_text.append(
						str(len(gold_mapping)) + ' BCs for mapping of library ' + lib + '.\n'
						)
				print(lines_of_text[len(lines_of_text) - 1])
					
			if lib_dictionary[lib]['normalization']:
				gold_normalization, rep_dictionary[lib]['normalization'], cont_dictionary[lib]['normalization'] = form_gold_norm_or_expr(
						ind_dictionary, lib_dictionary, lib,
						'normalization'
						)
				lines_of_text.append(
						str(len(gold_normalization)) + ' BCs for normalization of library ' + lib + '.\n'
						)
				print(lines_of_text[len(lines_of_text) - 1])
					
			if lib_dictionary[lib]['expression']:
				if param.minimum_read_count['expression'] > 0:
					gold_expression, rep_dictionary[lib]['expression'], cont_dictionary[lib]['expression'] = form_gold_norm_or_expr(
							ind_dictionary, lib_dictionary, lib, 'expression'
							)
					lines_of_text.append(
							str(len(gold_expression)) + ' BCs for expression of library ' + lib + '.\n'
							)
					print(lines_of_text[len(lines_of_text) - 1])
				else:
					rep_dictionary[lib]['expression'], cont_dictionary[lib]['expression'], total_expr = form_rep_expr(
							ind_dictionary, lib_dictionary, lib
							)

			if gold_mapping and gold_normalization and gold_expression:
				for barcode in gold_mapping:
					if barcode in gold_normalization and barcode in gold_expression:
						barcode_count = gold_normalization[barcode] + gold_expression[barcode]
						for ROI in gold_mapping[barcode]:
							barcode_count += gold_mapping[barcode][ROI]
						gold_dictionary[lib]['barcodes'][barcode] = barcode_count
						gold_dictionary[lib]['barcodes-ROIs'][barcode] = gold_mapping[barcode]
				not_empty(
						gold_dictionary[lib]['barcodes'],
						'mapping, normalization and expression',
						lib
						)
			elif gold_mapping and gold_normalization:
				for barcode in gold_mapping:
					if barcode in gold_normalization:
						barcode_count = gold_normalization[barcode]
						for ROI in gold_mapping[barcode]:
							barcode_count += gold_mapping[barcode][ROI]
						if barcode in total_expr:
							barcode_count += total_expr[barcode]
						gold_dictionary[lib]['barcodes'][barcode] = barcode_count
						gold_dictionary[lib]['barcodes-ROIs'][barcode] = gold_mapping[barcode]
				not_empty(
						gold_dictionary[lib]['barcodes'],
						'mapping and normalization',
						lib
						)
			elif gold_mapping and gold_expression:
				for barcode in gold_mapping:
					if barcode in gold_expression:
						barcode_count = gold_expression[barcode]
						for ROI in gold_mapping[barcode]:
							barcode_count += gold_mapping[barcode][ROI]
						gold_dictionary[lib]['barcodes'][barcode] = barcode_count
						gold_dictionary[lib]['barcodes-ROIs'][barcode] = gold_mapping[barcode]
				not_empty(
						gold_dictionary[lib]['barcodes'],
						'mapping and expression',
						lib
						)
			elif gold_normalization and gold_expression:
				for barcode in gold_normalization:
					if barcode in gold_expression:
					  barcode_count = gold_normalization[barcode] + gold_expression[barcode]
					  gold_dictionary[lib]['barcodes'][barcode] = barcode_count
				not_empty(
						gold_dictionary[lib]['barcodes'],
						'normalization and expression',
						lib
						)
			elif gold_mapping:
				for barcode in gold_mapping:
					barcode_count = 0
					for ROI in gold_mapping[barcode]:
						barcode_count += gold_mapping[barcode][ROI]
					if barcode in total_expr:
						barcode_count += total_expr[barcode]
					gold_dictionary[lib]['barcodes'][barcode] = barcode_count
					gold_dictionary[lib]['barcodes-ROIs'][barcode] = gold_mapping[barcode]
			else:
				for barcode in gold_normalization:
					gold_dictionary[lib]['barcodes'][barcode] = gold_normalization[barcode]
				for barcode in gold_expression:
					gold_dictionary[lib]['barcodes'][barcode] = gold_expression[barcode]
				for barcode in gold_dictionary[lib]['barcodes']:
					if barcode in total_expr:
						gold_dictionary[lib]['barcodes'][barcode] += total_expr[barcode]
							
			if gold_dictionary[lib]['barcodes']:
				lines_of_text.append(
						'%s unique BCs for library %s.\n' %(str(len(gold_dictionary[lib]['barcodes'])), lib)
						)
				print(lines_of_text[len(lines_of_text) - 1])
				with open(run_info_text, 'a') as info_file:
					info_file.writelines(lines_of_text)
					
			total_rep_dictionary[lib] = {
					'normalization': gold_normalization,
					'expression': gold_expression,
					'mapping': gold_mapping
					}
		else:
			file_count = 0
			for file in lib_dictionary['H']:
				file_count += 1
				for index in lib_dictionary['H'][file]:
					df_dict = {'Barcode': [], 'ROI': [], 'Count': []}
					for barcode in ind_dictionary[file][index]:
						if barcode not in ['read count', 'low quality count', 'eff read count']:
							for ROI in ind_dictionary[file][index][barcode]:
								df_dict['Barcode'].append(barcode)
								df_dict['ROI'].append(ROI)
								df_dict['Count'].append(len(ind_dictionary[file][index][barcode][ROI]))
					df = pd.DataFrame(data=df_dict)
					table_path = output_folder + '/' + 'Chimera_formation_test' + index + str(file_count) + '.xlsx'
					table = df.to_excel(table_path, index = None)
			
	with open(output_folder + '/gold_dictionary.json', 'w') as f:
		json.dump(gold_dictionary, f)
		
	with open(output_folder + '/rep_dictionary.json', 'w') as fi:
		json.dump(rep_dictionary, fi)
		
	with open(output_folder + '/cont_dictionary.json', 'w') as fil:
		json.dump(cont_dictionary, fil)
		
	with open(output_folder + '/total_rep_dictionary.json', 'w') as file:
		json.dump(total_rep_dictionary, file)
			
	return gold_dictionary, rep_dictionary, cont_dictionary, total_rep_dictionary


gold_dict, rep_dict, cont_dict, total_rep_dict = form_gold_dictionary(fastq_file_dictionary, lib_diction)

#Form gold dictionary of genuine barcodes
#(only works if minimum_read_count > 0).
#gold_dictionary = {
#				 '29-36': {
#						   'barcodes': {'barcode1': barcode1_count1, 'barcode2': barcode1_count2},
#						   'barcodes-ROIs': {
#												   'barcode1': {'mut1': mut1_count, 'mut2': mut2_count},
#												   'barcode2': {'mut3': mut3_count, 'mut4': mut4_count}
#												   }
#						  }
#				 '25-32': {
#						   'barcodes': {'barcode3': barcode1_count3, 'barcode4': barcode1_count4},
#						   'barcodes-ROIs': {}
#						   }
#				   }



#Form gold dictionary of unique barcodes
#uniq_gold_dictionary = {
#				 '29-36': {
#						   'barcodes': {'barcode1': ['barcode2']},
#						   'barcodes-ROIs': {
#											'barcode1':{
#												   'barcode1': ['mut1', 'mut2'],
#												   'barcode2': ['mut3', 'mut4'],
#												   'main ROI' : {
#															  'mut1': ['mut3']
#																	  }
#												   'total reads': n,
#												   'chimeric reads': k
#												   }
#												 }
#						   }
#				 '25-32': {
#						   'barcodes': {'barcode3': [], 'barcode4': []},
#						   'barcodes-ROIs': {}
#						   }
#						 }
def form_uniq_gold_dictionary(gold_dictionary):
	print(datetime.now())
	uniq_gold_dictionary = {}
	for lib in gold_dictionary:
		if gold_dictionary[lib]['barcodes']:
			uniq_gold_dictionary[lib] = {
					'barcodes': {}, 'barcodes-ROIs': {}
					}
			hash_dict = {}
			uniq_barcodes_list = []
			mut_variant_dict = {}
			
#Form hash dictionary:
#hash_dict = {
#			 'ATGCAG': [
#						('TGACTGGATCGAATGCAG', 5),
#						('ATGCAGAAGGCTCTTGAT', 3),
#						('ATGCAGAAGGCTCTTGAC', 1),
#					   ],
#			  'TGCAGA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'GCAGAA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'CAGAAG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'AGAAGG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'GAAGGC': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'AAGGCT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'AGGCTC': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'GGCTCT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'GCTCTT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'CTCTTG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'TCTTGA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#			  'CTTGAT': [('ATGCAGAAGGCTCTTGAT', 3)],
#			  'CTTGAC': [('ATGCAGAAGGCTCTTGAC', 1)],
#			  'TGACTG': [('TGACTGGATCGAATGCAG', 5)],
#			  'GACTGG': [('TGACTGGATCGAATGCAG', 5)],
#			  'ACTGGA': [('TGACTGGATCGAATGCAG', 5)],
#			  'CTGGAT': [('TGACTGGATCGAATGCAG', 5)],
#			  'TGGATC': [('TGACTGGATCGAATGCAG', 5)],
#			  'GGATCG': [('TGACTGGATCGAATGCAG', 5)],
#			  'GATCGA': [('TGACTGGATCGAATGCAG', 5)],
#			  'ATCGAA': [('TGACTGGATCGAATGCAG', 5)],
#			  'TCGAAT': [('TGACTGGATCGAATGCAG', 5)],
#			  'CGAATG': [('TGACTGGATCGAATGCAG', 5)],
#			  'GAATGC': [('TGACTGGATCGAATGCAG', 5)],
#			  'AATGCA': [('TGACTGGATCGAATGCAG', 5)],
#			  'ATGCAG': [('TGACTGGATCGAATGCAG', 5)]
#			 }
			for barcode in gold_dictionary[lib]['barcodes']:
				for i in range(len(barcode) - param.hash_length):
					new_hash = barcode[i:i + param.hash_length]
					if new_hash not in hash_dict:
						hash_dict[new_hash] = [(
								barcode,
								gold_dictionary[lib]['barcodes'][barcode]
								)]
					else:
						hash_dict[new_hash].append((
								barcode,
								gold_dictionary[lib]['barcodes'][barcode]
								))
							
#Form list of barcodes and their mutant variants:
#uniq_barcodes_list = [
#					  ['TGACTGGATCGAATGCAG', 5],
#					  ['ATGCAGAAGGCTCTTGAT', 3, 'ATGCAGAAGGCTCTTGAC'],
#					  ['ATGCAGAAGGCTCTTGAC', 1]
#					 ]
			for new_hash in hash_dict:
				hash_dict[new_hash].sort(
						key=lambda x: (x[1], x[0]), reverse=True
						)
				for i in range(len(hash_dict[new_hash]) - 1):
					high_count_barcode = hash_dict[new_hash][i][0]
					
					if not uniq_barcodes_list:
						high_count_barcode_list = [
								high_count_barcode, hash_dict[new_hash][i][1]
								]
						uniq_barcodes_list.append(high_count_barcode_list)
					else:
						for barcode_list in uniq_barcodes_list:
							if high_count_barcode == barcode_list[0]:
								high_count_barcode_list = barcode_list
								break
						else:
							high_count_barcode_list = [
									high_count_barcode,
									hash_dict[new_hash][i][1]
									]
							uniq_barcodes_list.append(high_count_barcode_list)
							
					for n in range(i + 1, len(hash_dict[new_hash])):
						low_count_barcode = hash_dict[new_hash][n][0]
						if low_count_barcode not in high_count_barcode_list:
							if align_mut_variants(
									high_count_barcode,
									low_count_barcode,
									'barcode'
									):
								high_count_barcode_index = uniq_barcodes_list.index(
										high_count_barcode_list
										)
								uniq_barcodes_list[high_count_barcode_index].append(
										low_count_barcode
										)
			
			uniq_barcodes_list.sort(key=lambda x: (x[1], x[0]), reverse=True)
			
#Form dictionary of unique barcodes:
#uniq_gold_dictionary = {'29-36': {
#								  'barcodes': {
#											 'TGACTGGATCGAATGCAG': [],
#											 'ATGCAGAAGGCTCTTGAT':
#														['ATGCAGAAGGCTCTTGAC']
#											  },
#								  'barcodes-ROIs': {}
#								  }
#					   }
			for barcode_list in uniq_barcodes_list:
				if barcode_list[0] not in mut_variant_dict:
					uniq_gold_dictionary[lib]['barcodes'][barcode_list[0]] = []
					if len(barcode_list) > 2:
						for i in range(2, len(barcode_list)):
							uniq_gold_dictionary[lib]['barcodes'][barcode_list[0]].append(
									barcode_list[i]
									)
							if barcode_list[i] in mut_variant_dict:
								mut_variant_dict[barcode_list[i]].append(
										barcode_list[0]
										)
							else:
								mut_variant_dict[barcode_list[i]] = [barcode_list[0]]

#Check if mutant barcodes are the same for different unique barcodes.
#If so, check ROIs.
#If such unique barcodes have similar ROIs, they will be combined.
			for mut_variant in mut_variant_dict:
				if len(mut_variant_dict[mut_variant]) > 1:
					if gold_dictionary[lib]['barcodes-ROIs']:
						same_barcodes = []
						for barcode in mut_variant_dict[mut_variant]:
							for ROI in gold_dictionary[lib]['barcodes-ROIs'][mut_variant]:
								for mut in gold_dictionary[lib]['barcodes-ROIs'][barcode]:
									if align_mut_variants(
											ROI, mut, 'ROI'
											):
									   same_barcodes.append([
											   barcode,
											   gold_dictionary[lib]['barcodes'][barcode]
											   ])
									   break
								else:
									continue
								break
							else:
								if mut_variant in uniq_gold_dictionary[lib]['barcodes'][barcode]:
									uniq_gold_dictionary[lib]['barcodes'][barcode].remove(
											mut_variant
											)
						if len(same_barcodes) > 1:
							same_barcodes.sort(
									key=lambda x: (x[1], x[0]), reverse=True
									)
							uniq_barcode = same_barcodes[0][0]
							for i in range(1, len(same_barcodes)):
								if same_barcodes[i] not in uniq_gold_dictionary[lib]['barcodes'][uniq_barcode]:
									uniq_gold_dictionary[lib]['barcodes'][uniq_barcode].append(
											same_barcodes[i]
											)
								if same_barcodes[i] in uniq_gold_dictionary[lib]['barcodes']:
									for mut_same_barcode in uniq_gold_dictionary[lib]['barcodes'][same_barcodes[i]]:
										if mut_same_barcode not in uniq_gold_dictionary[lib]['barcodes'][uniq_barcode]:
											uniq_gold_dictionary[lib]['barcodes'][uniq_barcode].append(
													mut_same_barcode
													)
									uniq_gold_dictionary[lib]['barcodes'].pop(
											same_barcodes[i]
											)
					else:
						for barcode in mut_variant_dict[mut_variant]:
							if mut_variant in uniq_gold_dictionary[lib]['barcodes'][barcode]:
								uniq_gold_dictionary[lib]['barcodes'][barcode].remove(
										mut_variant
										)

#Form dictionary 'barcodes-ROIs'.
#Find main ROI for each unique barcode.
			barcodes_to_remove = []
			if gold_dictionary[lib]['barcodes-ROIs']:
				for uniq_barcode in uniq_gold_dictionary[lib]['barcodes']:
					uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode] = {
							uniq_barcode: []
							}
					ROI_lists = []
					for ROI in gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]:
						ROI_lists.append([
								ROI,
								gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode][ROI]
								])
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode][uniq_barcode].append(
								ROI
								)
					for barcode_mut_var in uniq_gold_dictionary[lib]['barcodes'][uniq_barcode]:
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode][barcode_mut_var] = []
						for mut in gold_dictionary[lib]['barcodes-ROIs'][barcode_mut_var]:
							for ROI_list in ROI_lists:
								if mut == ROI_list[0]:
									ROI_list[1] += gold_dictionary[lib]['barcodes-ROIs'][barcode_mut_var][mut]
									break
							else:
								ROI_lists.append([
										mut,
										gold_dictionary[lib]['barcodes-ROIs'][barcode_mut_var][mut]
										])
					ROI_lists.sort(key=lambda x: (x[1], x[0]), reverse=True)
					total_ROI_count = ROI_lists[len(ROI_lists) - 1][1]
					for i in range(len(ROI_lists) - 1):
						total_ROI_count += ROI_lists[i][1]
						if ROI_lists[i][0] != 'mut var':
							for n in range(i + 1, len(ROI_lists)):
								if align_mut_variants(
										ROI_lists[i][0],
										ROI_lists[n][0],
										'ROI'):
									ROI_lists[i][1] += ROI_lists[n][1]
									ROI_lists[i].append(ROI_lists[n][0])
									ROI_lists[n][0] = 'mut var'
					ROI_lists.sort(key=lambda x: (x[1], x[0]), reverse=True)
					main_ROI = ROI_lists[0][0]
					if main_ROI in gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode] and ROI_lists[0][1] >= total_ROI_count * param.genuine_ROI_portion:
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['main ROI'] = {
								main_ROI: []
								}
						if len(ROI_lists[0]) > 2:
							for k in range(2, len(ROI_lists[0])):
								uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['main ROI'][main_ROI].append(
										ROI_lists[0][k]
										)
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['total reads'] = total_ROI_count
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['chimeric reads'] = round(total_ROI_count - ROI_lists[0][1], 6)
					else:
						uniq_gold_dictionary[lib]['barcodes-ROIs'].pop(
								uniq_barcode
								)
						barcodes_to_remove.append(uniq_barcode)
				for barcode_to_remove in barcodes_to_remove:
					uniq_gold_dictionary[lib]['barcodes'].pop(
							barcode_to_remove
							)
			lines_of_text = [
					str(len(uniq_gold_dictionary[lib]['barcodes'])) + ' genuine BCs for library ' + lib + '.\n',
					str(len(barcodes_to_remove)) + ' BCs were removed from analysis due to high percent of chimeric molecules.\n'
					]
			print(lines_of_text)
			with open(run_info_text, 'a') as info_file:
				info_file.writelines(lines_of_text)
			print(datetime.now())
			
	if not os.path.isfile(output_folder + '/' + 'uniq_gold_dictionary.json'):
		with open(output_folder + '/uniq_gold_dictionary.json', 'w') as f:
			json.dump(uniq_gold_dictionary, f)
		
	return uniq_gold_dictionary
							

gold_dict = form_uniq_gold_dictionary(gold_dict)