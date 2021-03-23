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

#Reading json from a File 
with open('uniq_gold_dictionary.json') as json_file:
	uniq_gold_dictionary = json.load(json_file)
with open('rep_dictionary.json') as json_file:
	rep_dictionary = json.load(json_file)
with open('cont_dictionary.json') as json_file:
	cont_dictionary = json.load(json_file)
with open('total_rep_dictionary.json') as json_file:
	total_rep_dictionary = json.load(json_file)
with open('gold_dictionary.json') as json_file:
	gold_dictionary = json.load(json_file)


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

#Form gold dictionary of unique barcodes
#uniq_gold_dictionary = {
#				 '29-36': {
#						   'barcodes': {'barcode1': ['barcode2']},
#						   'barcodes-ROIs': {
#											'barcode1':{
#												   'barcode1': ['mut1', 'mut2'],
#												   'barcode2': ['mut3', 'mut4'],
#												   'major ROI' : {
#															  'mut1': ['mut3']
#																	  }
#												   'total reads': n,
#												   'hybrid reads': k
#												   }
#												 }
#						   }
#				 '25-32': {
#						   'barcodes': {'barcode3': [], 'barcode4': []},
#						   'barcodes-ROIs': {}
#						   }
#						 }
def form_uniq_gold_dictionary(gold_dictionary):
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

#Check if mutant barcode are the same for different unique barcodes.
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
#Find major ROI for each unique barcode.
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
					major_ROI = ROI_lists[0][0]
					if major_ROI in gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode] and ROI_lists[0][1] >= total_ROI_count * param.genuine_ROI_portion:
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['major ROI'] = {
								major_ROI: []
								}
						if len(ROI_lists[0]) > 2:
							for k in range(2, len(ROI_lists[0])):
								uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['major ROI'][major_ROI].append(
										ROI_lists[0][k]
										)
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['total reads'] = total_ROI_count
						uniq_gold_dictionary[lib]['barcodes-ROIs'][uniq_barcode]['hybrid reads'] = round(total_ROI_count - ROI_lists[0][1], 6)
					else:
						uniq_gold_dictionary[lib]['barcodes-ROIs'].pop(
								uniq_barcode
								)
						barcodes_to_remove.append(uniq_barcode)
				for barcode_to_remove in barcodes_to_remove:
					uniq_gold_dictionary[lib]['barcodes'].pop(
							barcode_to_remove
							)
	if not os.path.isfile('uniq_gold_dictionary.json'):
		with open('uniq_gold_dictionary.json', 'w') as f:
			json.dump(uniq_gold_dictionary, f)
		
	return uniq_gold_dictionary
							

gold_dictionary = form_uniq_gold_dictionary(gold_dictionary)


# Function to calculate correlation coefficients between two arrays
def corr(x, y, **kwargs):
	# Calculate the values
	coef, p_val = stats.pearsonr(x, y)
	sp, sp_p_val = stats.spearmanr(x, y)
	# Make the label
	label = r'$\rho$ = %s, pearson p-val = %s, sp = %s, spearman p-val = %s.' %(str(round(coef, 2)), str(round(p_val, 2)), str(round(sp, 2)), str(round(sp_p_val, 2)))
	
	# Add the label to the plot
	ax = plt.gca()
	ax.annotate(label, xy = (0.05, 0.95), size = 20, xycoords = ax.transAxes)


def make_tables(
		uniq_gold_dictionary,
		rep_dictionary,
		cont_dictionary,
		total_rep_dictionary
		):
	mut_dict = {}
	mut_log2_dict = {}
	wt_mut_log2 = {}
	log_2_wt_gaus = {}
	lines_of_text = []
	for lib in rep_dictionary:
		df_cont_dict = {
				'Barcode': ['control1', 'control2', 'experiment1', 'experiment2'],
				'ROI': ['control1', 'control2', 'experiment1', 'experiment2']
				}
		
		if uniq_gold_dictionary[lib]['barcodes']:
			if uniq_gold_dictionary[lib]['barcodes-ROIs']:
				df_lib_dict = {'Barcode': [], 'ROI': []}
				for barcode in uniq_gold_dictionary[lib]['barcodes-ROIs']:
					df_lib_dict['Barcode'].append(barcode)
					for ROI in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['major ROI']:
						df_lib_dict['ROI'].append(ROI)
			else:
				df_lib_dict = {'Barcode': []}
				for barcode in uniq_gold_dictionary[lib]['barcodes']:
					df_lib_dict['Barcode'].append(barcode)
			unique_barcode_count = len(uniq_gold_dictionary[lib]['barcodes'])
			
			for rep_type in rep_dictionary[lib]:
				file_count = 0
				if rep_type != 'mapping' and rep_dictionary[lib][rep_type]:
					df_cont_dict[rep_type + ' mean'] = [0, 0, 0, 0]
				
				if rep_type != 'mapping' and rep_dictionary[lib][rep_type]:
					df_lib_dict[rep_type + ' mean'] = [
						0 for x in range(unique_barcode_count)
						]
					replicate_count = 0
					replicate_list = []
					
					for file in rep_dictionary[lib][rep_type]:
						file_count += 1
						for replicate in rep_dictionary[lib][rep_type][file]:
							if rep_type != 'mapping':
								replicate_list.append(
									rep_type + ' ' + replicate
									)
							df_dict = {
									'Genuine Barcode': [],
									'Variants of Barcode': [],
									'Count': [],
									'Total Count': []
									}
							df_lib_dict[rep_type + ' ' + replicate] = [
									'' for a in range(unique_barcode_count)
									]
							
							if rep_type != 'mapping':
								df_cont_dict[rep_type + ' ' + replicate] = [
										'', '', '', ''
										]
							barcode_number = 0
							for cont_barcode in df_cont_dict['Barcode']:
								if rep_type != 'mapping' and cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
									df_cont_dict[rep_type + ' ' + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
									df_cont_dict[rep_type + ' mean'][barcode_number] += cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
								barcode_number += 1
								
							for barcode in uniq_gold_dictionary[lib]['barcodes']:
								barcode_index = df_lib_dict['Barcode'].index(
										barcode
										)
								variants_of_barcode = barcode + '\n'
								if barcode in rep_dictionary[lib][rep_type][file][replicate]:
									barcode_count = str(
											rep_dictionary[lib][rep_type][file][replicate][barcode]
											) + '\n'
									total_count = rep_dictionary[lib][rep_type][file][replicate][barcode]
								else:
									barcode_count = '0\n'
									total_count = 0
								for bar_variant in uniq_gold_dictionary[lib]['barcodes'][barcode]:
									variants_of_barcode += bar_variant + '\n'
									if bar_variant in rep_dictionary[lib][rep_type][file][replicate]:
										barcode_count += str(rep_dictionary[lib][rep_type][file][replicate][bar_variant]) + '\n'
										total_count += rep_dictionary[lib][rep_type][file][replicate][bar_variant]
									else:
										barcode_count += '0\n'
										
								df_dict['Genuine Barcode'].append(barcode)
								df_dict['Variants of Barcode'].append(
										variants_of_barcode
										)
								df_dict['Count'].append(barcode_count)
								df_dict['Total Count'].append(total_count)
								df_lib_dict[rep_type + ' mean'][barcode_index] += total_count
								df_lib_dict[rep_type + ' ' + replicate][barcode_index] = total_count
								
							df = pd.DataFrame(df_dict, columns=[
									'Genuine Barcode', 'Variants of Barcode',
									'Count', 'Total Count'
									])
							df = df.sort_values(
									by=['Total Count'], ascending=False
									)
							table_path = lib + '_' + rep_type + '_' + replicate + '.xlsx'
							table = df.to_excel(table_path, index = None)
					if replicate_count == 0:
						continue
					df_lib_dict[rep_type + ' mean'] = [
							y / replicate_count for y in df_lib_dict[rep_type + ' mean']
							]
					if rep_type != 'mapping':
						df_cont_dict[rep_type + ' mean'] = [
								d / replicate_count for d in df_cont_dict[rep_type + ' mean']
								]
					
					scatter_df = pd.DataFrame(data=df_lib_dict)
					#Make joinplot.
					if replicate_count >= 2:
						plt.figure()
						sns.set_style('ticks')
						sns.set_context('poster')
						# Create a pair grid instance
						grid = sns.PairGrid(
								data=scatter_df,
								vars =replicate_list,
								height=10
								)

						# Map the plots to the locations
						grid = grid.map_upper(
								sns.regplot, line_kws={"color": "g"},
								scatter_kws={'s': 0.1}
								)
						grid = grid.map_upper(corr)
						grid = grid.map_lower(
								sns.kdeplot, shade=True
								)
						grid = grid.map_diag(
								sns.kdeplot, linewidth=3,
								shade=True
								)
						plt.savefig(
								'pair_plot_'+ rep_type + lib + '.pdf',
								format='pdf', dpi=1000
								)
						plt.close()
						
					if total_rep_dictionary[lib][rep_type]:
						df_total_rep_dict = {
									'Genuine Barcode': [],
									'Variants of Barcode': [],
									'Count': [],
									'Total Count': []
									}
						for barcode in uniq_gold_dictionary[lib]['barcodes']:
							var_of_barcode = barcode + '\n'
							barc_count = str(
									total_rep_dictionary[lib][rep_type][barcode] / replicate_count
									) + '\n'
							total_barc_count = total_rep_dictionary[lib][rep_type][barcode] / replicate_count
							for bar_variant in uniq_gold_dictionary[lib]['barcodes'][barcode]:
								var_of_barcode += bar_variant + '\n'
								barc_count += str(
										total_rep_dictionary[lib][rep_type][bar_variant] / replicate_count
										) + '\n'
								total_barc_count += total_rep_dictionary[lib][rep_type][bar_variant] / replicate_count
								
							df_total_rep_dict['Genuine Barcode'].append(barcode)
							df_total_rep_dict['Variants of Barcode'].append(
									var_of_barcode
									)
							df_total_rep_dict['Count'].append(barc_count)
							df_total_rep_dict['Total Count'].append(
									total_barc_count
									)
							
						df_total_rep = pd.DataFrame(df_total_rep_dict, columns=[
									'Genuine Barcode', 'Variants of Barcode',
									'Count', 'Total Count'
									])
						df_total_rep = df_total_rep.sort_values(
								by=['Total Count'], ascending=False
								)
						total_rep_path = lib + '_' + rep_type + '.xlsx'
						df_total_rep.to_excel(total_rep_path, index = None)
						
				else:
					replicate_count = 0
					for file in rep_dictionary[lib][rep_type]:
						file_count += 1
						for replicate in rep_dictionary[lib][rep_type][file]:
							replicate_count += 1
							df_dict = {
									'Genuine Barcode': [],
									'Main ROI': [],
									'Variants of Barcode': [],
									'ROIs' : [],
									'Count': [],
									'Total Count': [],
									'Chimeric molecules, read count': []
									}
							if rep_type != 'mapping':
								df_cont_dict[rep_type + ' ' + replicate] = [
										'', '', '', ''
										]
							
							barcode_number = 0
							for cont_barcode in df_cont_dict['Barcode']:
								if rep_type != 'mapping' and cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
									df_cont_dict[rep_type + ' ' + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
									df_cont_dict[rep_type + ' mean'][barcode_number] += cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
								barcode_number += 1
							
							for barcode in uniq_gold_dictionary[lib]['barcodes-ROIs']:
								total_count = 0
								hybrid_count = 0
								for ROI in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['major ROI']:
									main_ROI = ROI
								variants_of_barcode = barcode + '\n'
								ROIs = main_ROI + '\n'
								if barcode in rep_dictionary[lib][rep_type][file][replicate]:
									if main_ROI in rep_dictionary[lib][rep_type][file][replicate][barcode]:
										count = str(rep_dictionary[lib][rep_type][file][replicate][barcode][main_ROI]) + '\n'
										total_count += rep_dictionary[lib][rep_type][file][replicate][barcode][main_ROI]
									else:
										count = '0\n'
								else:
									count = '0\n'
								for ROI_variant in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode][barcode]:
									if ROI_variant != main_ROI:
										variants_of_barcode += '\n'
										ROIs += ROI_variant + '\n'
										if barcode in rep_dictionary[lib][rep_type][file][replicate]:
											if ROI_variant in rep_dictionary[lib][rep_type][file][replicate][barcode]:
												count += str(rep_dictionary[lib][rep_type][file][replicate][barcode][ROI_variant]) + '\n'
												total_count += rep_dictionary[lib][rep_type][file][replicate][barcode][ROI_variant]
												if ROI_variant not in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['major ROI'][main_ROI]:
													hybrid_count += rep_dictionary[lib][rep_type][file][replicate][barcode][ROI_variant]
											else:
												count += '0\n'
										else:
											count += '0\n'
								for barc_variant in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]:
									if len(barc_variant) > 15 and barc_variant != barcode:
										variants_of_barcode += barc_variant
										for mut_variant in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode][barc_variant]:
											variants_of_barcode += '\n'
											ROIs += mut_variant + '\n'
											if barc_variant in rep_dictionary[lib][rep_type][file][replicate]:
												if mut_variant in rep_dictionary[lib][rep_type][file][replicate][barc_variant]:
													count += str(rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant])
													total_count += rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant]
													if mut_variant != main_ROI and mut_variant not in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['major ROI'][main_ROI]:
														hybrid_count += rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant]
												else:
													count += '0\n'
											else:
												count += '0\n'
								df_dict['Genuine Barcode'].append(barcode)
								df_dict['Main ROI'].append(
										main_ROI
										)
								df_dict['Variants of Barcode'].append(
										variants_of_barcode
										)
								df_dict['ROIs'].append(ROIs)
								df_dict['Count'].append(count)
								df_dict['Total Count'].append(total_count)
								df_dict['Chimeric molecules, read count'].append(hybrid_count)
								
							df = pd.DataFrame(df_dict, columns=[
									'Genuine Barcode', 'Main ROI',
									'Variants of Barcode', 'ROIs',
									'Count', 'Total Count', 'Chimeric molecules, read count'
									])
							df['Chimeric molecules, %'] = df['Chimeric molecules, read count'] / df['Total Count'] * 100
							df = df.sort_values(
									by=['Total Count'], ascending=False
									)
							table_path = lib + '_mapping_' + replicate + '.xlsx'
							table = df.to_excel(table_path, index = None)
					if replicate_count == 0:
						continue		
					if rep_type != 'mapping':
						df_cont_dict[rep_type + ' mean'] = [
								d / replicate_count for d in df_cont_dict[rep_type + ' mean']
								]
					
					df_total_map_dict = {
									'Genuine Barcode': [],
									'Main ROI': [],
									'Variants of Barcode': [],
									'ROIs' : [],
									'Count': [],
									'Total Count': [],
									'Chimeric molecules, read count': [],
									'Chimeric molecules, %': []
									}
					for barcode in uniq_gold_dictionary[lib]['barcodes-ROIs']:
						for ROI in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['major ROI']:
							main_mut = ROI
						total_map_count = uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['total reads'] / replicate_count
						hybr_count = uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]['hybrid reads'] / replicate_count
						hybrid_percent = hybr_count / total_map_count * 100
						var_of_barc = barcode + '\n'
						mut = main_mut + '\n'
						map_count = str(total_rep_dictionary[lib]['mapping'][barcode][main_mut] / replicate_count) + '\n'
						for mut_var in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode][barcode]:
							if mut_var != main_mut:
								var_of_barc += '\n'
								mut += mut_var + '\n'
								map_count += str(total_rep_dictionary[lib]['mapping'][barcode][mut_var] / replicate_count) + '\n'
						for barc_var in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode]:
							if len(barc_var) > 15 and barc_var != barcode:
								var_of_barcode += barc_var
								for mut_variant in uniq_gold_dictionary[lib]['barcodes-ROIs'][barcode][barc_var]:
									var_of_barc += '\n'
									mut += mut_variant + '\n'
									map_count += str(total_rep_dictionary[lib]['mapping'][barc_var][mut_variant] / replicate_count) + '\n'
						df_total_map_dict['Genuine Barcode'].append(barcode)
						df_total_map_dict['Main ROI'].append(main_mut)
						df_total_map_dict['Variants of Barcode'].append(
								var_of_barc
								)
						df_total_map_dict['ROIs'].append(mut)
						df_total_map_dict['Count'].append(map_count)
						df_total_map_dict['Total Count'].append(
								total_map_count
								)
						df_total_map_dict['Chimeric molecules, read count'].append(hybr_count)
						df_total_map_dict['Chimeric molecules, %'].append(
								hybrid_percent
								)
						
					df_total_map = pd.DataFrame(df_total_map_dict, columns=[
									'Genuine Barcode', 'Main ROI',
									'Variants of Barcode', 'ROIs',
									'Count', 'Total Count', 'Chimeric molecules, read count',
									'Chimeric molecules, %'
									])
					df_total_map = df_total_map.sort_values(
							by=['Total Count'], ascending=False
							)
					total_map_table_path = lib + '_mapping ' + '.xlsx'
					total_map_table = df_total_map.to_excel(
							total_map_table_path, index = None
							)
			
			columns_list = [key for key in df_lib_dict]
			df_lib = pd.DataFrame(data=df_lib_dict)
			new_columns_list = ['Barcode']
			pos_start = 1
			
			if 'ROI' in columns_list:
				new_columns_list.append('ROI')
				pos_start += 1
			
			new_columns_list += [
					k for k in columns_list if k.startswith('normalization') and not k.endswith('mean')
					]
			new_columns_list += [
					n for n in columns_list if n.startswith('expression') and not n.endswith('mean')
					]
			if 'normalization mean' in columns_list:
				new_columns_list.append('normalization mean')
			if 'expression mean' in columns_list:
				new_columns_list.append('expression mean')
			df_lib = df_lib[new_columns_list]
						
			
			cont_new_columns_list = [
					'Barcode', 'ROI'
					] + new_columns_list[pos_start:]
					
			
			if 'normalization mean' in new_columns_list and 'expression mean' in new_columns_list:
				df_lib['normalyzed expression'] = df_lib['expression mean'] / df_lib['normalization mean']
				
				cont_new_columns_list.append('normalyzed expression')
				df_cont_dict['normalyzed expression'] = ['', '', '', '']
				for n in range(4):
					if df_cont_dict['normalization mean'][n] > 0:
						df_cont_dict['normalyzed expression'][n] = df_cont_dict['expression mean'][n] / df_cont_dict['normalization mean'][n]
				if '' not in df_cont_dict['normalyzed expression']:
					cont_new_columns_list.append('normalyzed expression mean')

					df_cont_dict['normalyzed expression mean'] = [
							(df_cont_dict['normalyzed expression'][0] + df_cont_dict['normalyzed expression'][1]) / 2,
							(df_cont_dict['normalyzed expression'][0] + df_cont_dict['normalyzed expression'][1]) / 2,
							(df_cont_dict['normalyzed expression'][2] + df_cont_dict['normalyzed expression'][3]) / 2,
							(df_cont_dict['normalyzed expression'][2] + df_cont_dict['normalyzed expression'][3]) / 2
							]
					if df_cont_dict['normalyzed expression mean'][0] > 0:
						cont_new_columns_list.append('difference')
						
						if len(param.reference_BCs['forward']) > 0:
							df_lib['normalyzed to reference'] = df_lib['normalyzed expression'] / df_cont_dict['normalyzed expression mean'][0]
						
						df_cont_dict['difference'] = [
								df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
								df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
								df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
								df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0]
								]
						
				df_lib = df_lib.sort_values(
						by=['normalyzed expression'], ascending=False
						)
			
			total_table_path = lib + '.xlsx'
			total_table = df_lib.to_excel(total_table_path, index = None)
			
			if 'ROI' in new_columns_list:
				#Pairplot for ROIs that have more than two barcodes.
				lib_dict = df_lib.to_dict('index')
				ROI_dictionary = {}
				for index in lib_dict:
					new_ROI = lib_dict[index]['ROI']
					if len(param.reference_BCs['forward']) > 0:
						if new_ROI in ROI_dictionary:
							ROI_dictionary[new_ROI].append(lib_dict[index]['normalyzed to reference'])
						else:
							ROI_dictionary[new_ROI] = [lib_dict[index]['normalyzed to reference']]
					if len(param.reference_BCs['forward']) == 0:
						if new_ROI in ROI_dictionary:
							ROI_dictionary[new_ROI].append(lib_dict[index]['normalyzed expression'])
						else:
							ROI_dictionary[new_ROI] = [lib_dict[index]['normalyzed expression']]

				first_mut_pairplot_list = []
				second_mut_pairplot_list = []
				

				mut_barcode_count_list = []
				for new_ROI in ROI_dictionary:
					mut_barcode_count_list.append(
							len(ROI_dictionary[new_ROI])
							)
					if len(ROI_dictionary[new_ROI]) >= 2:
						first_random_expr = randint(0, len(ROI_dictionary[new_ROI]) - 1)
						first_mut_pairplot_list.append(ROI_dictionary[new_ROI][first_random_expr])
						second_random_expr = randint(0, len(ROI_dictionary[new_ROI]) - 1)
						while second_random_expr == first_random_expr:
							second_random_expr = randint(0, len(ROI_dictionary[new_ROI]) - 1)
						second_mut_pairplot_list.append(ROI_dictionary[new_ROI][second_random_expr])
						


				less_11 = 0
				for m in range(1, 11):
					lines_of_text.append(
							'%s ROIs with %s BCs for lib %s.\n' %(str(mut_barcode_count_list.count(m)), str(m), lib))
					less_11 += mut_barcode_count_list.count(m)
				lines_of_text.append(
						str(len(mut_barcode_count_list) - less_11) + ' ROIs with more than 10 BCs for lib ' + lib + '.\n')
				
				plt.figure(figsize=(10,10))
				sns.set_style('ticks')
				sns.set_context('poster')
				sns.displot(
						mut_barcode_count_list, kde=False,
						bins=max(mut_barcode_count_list)
						)
				sns.despine()
				plt.xlabel("Number of BCs for each ROI for lib " + lib, fontsize=12)
				plt.savefig(
						'BC_per_ROI_' + lib + '.pdf',
						format='pdf', dpi=100
						)
				plt.close()
				
				plt.figure(figsize=(10,10))
				sns.set_style('ticks')
				sns.set_context('poster')
				sns.displot(
						mut_barcode_count_list, kde=False,
						bins=max(mut_barcode_count_list)
						)
				sns.despine()
				plt.xlim(0, 10)
				plt.xlabel("Number of BCs for each ROI for " + lib, fontsize=12)
				plt.savefig(
						'BC_per_ROI_' + lib + '_limit10' + '.pdf',
						format='pdf', dpi=100
						)
				plt.close()

				pairplot_dict = {
						'first barcode': first_mut_pairplot_list,
						'second barcode': second_mut_pairplot_list
						}
				df_pairplot_dict = pd.DataFrame(
						pairplot_dict,
						columns=['first barcode', 'second barcode']
						)
				plt.figure()
				sns.set_style('ticks')
				sns.set_context('poster')
				# Create a pair grid instance
				grid = sns.PairGrid(
						data=df_pairplot_dict,
						vars = ['first barcode', 'second barcode'],
						height=10
						)

				# Map the plots to the locations
				grid = grid.map_upper(
						sns.regplot, line_kws={"color": "g"},
						scatter_kws={'s': 0.1}
						)
				grid = grid.map_upper(corr)
				grid = grid.map_lower(
						sns.kdeplot, shade=True
						)
				grid = grid.map_diag(
						sns.kdeplot, linewidth=3, shade=True
						)
				
				plt.savefig(
						'ROI_with_different_BCs_' + lib + '.pdf',
						format='pdf', dpi=1000
						)
				plt.close()
				
				#Plot mean normalysed expressions for all ROIs.
				if len(param.reference_BCs['forward']) > 0:
					if df_cont_dict['normalyzed expression mean'][0] > 0:
						mut_ser = df_lib.groupby('ROI')['normalyzed to reference'].mean()
					else:
						mut_ser = df_lib.groupby('ROI')['normalyzed expression'].mean()
				if len(param.reference_BCs['forward']) == 0:			
					mut_ser = df_lib.groupby('ROI')['normalyzed expression'].mean()
				mut_ser = mut_ser.sort_values(ascending=False)
				mut_table_path = lib + '_ROIs' + '.xlsx'
				mut_table = mut_ser.to_excel(mut_table_path)
				
				plt.figure(figsize=(10,10))
				sns.set_style('ticks')
				sns.set_context('poster')
				mut_ser_list = mut_ser.tolist()
				sns.displot(
						mut_ser, bins=round(mut_ser_list[0] * 10),
						kde=True
						)
				if len(param.reference_BCs['forward']) > 0:
					if 'normalyzed expression mean' in df_cont_dict:
						plt.axvline(1, 0,1, linewidth=2, color='r')
						plt.axvline(
								df_cont_dict['difference'][0], 0,1, linewidth=2,
								color='g'
								)
				sns.despine()
				plt.savefig(
						'density_plot_' + lib + '.pdf',
						format='pdf', dpi=1000
						)
				if param.wt_ROI[lib] in ROI_dictionary:
					lines_of_text.append(
							'In library %s %s wt ROIs.\n' %(lib, str(len(ROI_dictionary[param.wt_ROI[lib]])))
							)
					sns.displot(
							ROI_dictionary[param.wt_ROI[lib]],
							)
					sns.despine()
					plt.savefig(
						'wt_ROI_normalyzed_expression_' + lib + '.pdf',
						format='pdf', dpi=1000
						)
				plt.close()
				
				log_2_mut = [math.log2(x) for x in mut_ser_list if x > 0]
				stat, p = stats.shapiro(log_2_mut)
				if p > 0.05:
					lines_of_text.append(
							'For ROIs in lib %s log2 of normalyzed expressions looks Gaussian.\nShapiro p-value=%s.\n' % (lib, str(p))
							)
				else:
					statis, p_value = stats.normaltest(log_2_mut)
					if p_value > 0.05:
						lines_of_text.append(
								'For ROIs in lib %s log2 of normalyzed expressions looks Gaussian.\nShapiro p-value=%s, K^2 p-value=%s.\n' % (lib, str(p), str(p_value))
								)
					else:
						lines_of_text.append(
								'For ROIs in lib %s log2 of normalyzed expressions does not look Gaussian.\nShapiro p-value=%s, K^2 p-value=%s.\n' % (lib, str(p), str(p_value))
								)
					
				plt.figure(figsize=(10,10))
				sns.set_style('ticks')
				sns.set_context('poster')
				sns.displot(
						log_2_mut, bins=round(log_2_mut[0] * 20),
						kde=True
						)
				if len(param.reference_BCs['forward']) > 0:
					if 'normalyzed expression mean' in df_cont_dict:
						plt.axvline(0, 0,1, linewidth=2, color='r')
						plt.axvline(
								math.log2(df_cont_dict['difference'][0]), 0,1,
								linewidth=2,
								color='g'
								)
				sns.despine()
				plt.savefig(
						'log2_density_plot_' + lib + '.pdf',
						format='pdf', dpi=1000
						)
				if param.wt_ROI[lib] in ROI_dictionary:
					wt_mut_log2[lib] = [math.log2(y) for y in ROI_dictionary[param.wt_ROI[lib]] if y > 0]
					lines_of_text.append('In library %s %s wt ROIs above zero.\n' %(lib, str(len(wt_mut_log2[lib]))))
					statistic, p_val = stats.shapiro(wt_mut_log2[lib])
					if p_val > 0.05:
						lines_of_text.append('For wt in lib %s log2 of normalyzed expressions looks Gaussian.\np-value=%s.\n' % (lib, str(p_val)))
						log_2_wt_gaus[lib] = True
					else:
						lines_of_text.append('For wt in lib %s log2 of normalyzed expressions does not looks Gaussian.\np-value=%s.\n' % (lib, str(p_val)))
						log_2_wt_gaus[lib] = False
					
					sns.displot(
							wt_mut_log2[lib],
							)
					sns.despine()
					plt.savefig(
						'log2_wt_ROI_normalyzed_expression_' + lib + '.pdf',
						format='pdf', dpi=1000
						)
				plt.close()
				
				mut_dict[lib] = list(zip(mut_ser.index, mut_ser_list))
				mut_log2_dict[lib] = list(zip(mut_ser.index, log_2_mut))

			
			df_cont = pd.DataFrame(data=df_cont_dict)
			cont_table_path = lib + '_reference' + '.xlsx'
			if len(param.reference_BCs['forward']) > 0:
				cont_table = df_cont.to_excel(cont_table_path, index = None)
			
		else:
			for rep_type in rep_dictionary[lib]:
				if len(param.reference_BCs['forward']) > 0: 
					df_cont_dict[rep_type + ' mean'] = [0, 0, 0, 0]
					file_count = 0
					for file in rep_dictionary[lib][rep_type]:
						file_count += 1
						for replicate in rep_dictionary[lib][rep_type][file]:
							df_cont_dict[rep_type + ' ' + replicate] = [
										'', '', '', ''
										]
							barcode_number = 0
							for cont_barcode in df_cont_dict['Barcode']:
								if cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
									df_cont_dict[rep_type + ' ' + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
									df_cont_dict[rep_type + ' mean'][barcode_number] += cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
								barcode_number += 1
						
						rep_barcodes = {
								lib: {'barcodes': {}, 'barcodes-ROIs': {}}
								}
						rep_barcodes[lib]['barcodes'] = rep_dictionary[lib][rep_type][file][replicate]
						uniq_rep_dictionary = form_uniq_gold_dictionary(
								rep_barcodes
								)
						for barcode in uniq_rep_dictionary[lib]['barcodes']:
							variants_of_barcode = barcode + '\n'
							barcode_count = str(
									rep_dictionary[lib][rep_type][file][replicate][barcode]
									) + '\n'
							total_count = rep_dictionary[lib][rep_type][file][replicate][barcode]
							for bar_variant in uniq_rep_dictionary[lib]['barcodes'][barcode]:
								variants_of_barcode += bar_variant + '\n'
								barcode_count += str(rep_dictionary[lib][rep_type][file][replicate][bar_variant]) + '\n'
								total_count += rep_dictionary[lib][rep_type][file][replicate][bar_variant]
							df_dict['Genuine Barcode'].append(barcode)
							df_dict['Variants of Barcode'].append(
									variants_of_barcode
									)
							df_dict['Count'].append(barcode_count)
							df_dict['Total Count'].append(total_count)
						df = pd.DataFrame(data=df_dict)
						df = df.sort_values(
								by=['Total Count'], ascending=False
								)
						table_path = lib + '_' + rep_type + '_' + replicate + '.xlsx'
						table = df.to_excel(table_path, index = None)
			
			if len(param.reference_BCs['forward']) > 0:
				df_cont = pd.DataFrame(data=df_cont_dict)
				cont_table_path = lib + '_reference' + '.xlsx'
				cont_table = df_cont.to_excel(cont_table_path, index = None)
	
	print(lines_of_text)
	with open('run_info.txt', 'a') as info_file:
		info_file.writelines(lines_of_text)
	
	with open('ROI_dictionary.json', 'w') as f:
		json.dump(mut_dict, f)
	
	with open('ROI_log2_dictionary.json', 'w') as file:
		json.dump(mut_log2_dict, file)
	
	with open('wt_log2.json', 'w') as fil:
		json.dump(wt_mut_log2, fil)
		
	return mut_dict, mut_log2_dict, wt_mut_log2, log_2_wt_gaus


mut_dictionary, mut_log2_dictionary, wt_log2, log_2_wt_gaussian = make_tables(gold_dictionary, rep_dictionary, cont_dictionary, total_rep_dictionary)


def make_graphs(uniq_gold_dictionary):
	for lib in uniq_gold_dictionary:
		if uniq_gold_dictionary[lib]['barcodes']:
			barcode_len_list = [
					len(x) for x in uniq_gold_dictionary[lib]['barcodes']
					]
			
			plt.figure(figsize=(10,10))
			sns.set_style('ticks')
			sns.set_context('poster')
			sns.displot(
					barcode_len_list, kde=False, bins=4
					)
			sns.despine()
			plt.xlabel("BC length for " + lib, fontsize=14)
			plt.savefig(
					'BC_length_' + lib + '.pdf',
					format='pdf', dpi=1000
					)
			plt.close()


make_graphs(gold_dictionary)
				

def make_pwm(mut_dict):
	lines_of_text = []
	rounded_mismatch = round(param.ROI_mismatch)
	for lib in mut_dict:
		nucleotide_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
		total_count = 0
		for mut_tuple in mut_dict[lib]:
			for nucleotide in mut_tuple[0]:
				nucleotide_count[nucleotide] += 1
				total_count += 1
		lines_of_text.append('ROIs in lib ' + lib + ' contains\n')
		print('ROIs in lib ' + lib + ' contains ')
		for nucl in nucleotide_count:
			nucleotide_count[nucl] = nucleotide_count[nucl] / total_count
			lines_of_text.append(
					str(nucleotide_count[nucl]) + ' of ' + nucl + '.\n'
					)
			print(str(nucleotide_count[nucl]) + ' of ' + nucl)
		for i in range(
				param.ROI_length - rounded_mismatch, param.ROI_length + rounded_mismatch + 1
				):
			instances = [Seq(k[0]) for k in mut_dict[lib] if len(k[0]) == i]
			lines_of_text.append(
					str(len(instances)) + ' ROIs ' + str(i) + ' bp long for lib ' + lib + '.\n'
					)
			print(str(len(instances)) + ' ROIs ' + str(i) + ' bp long for lib ' + lib)
			if instances:
				m = motifs.create(instances)
				pwm = m.counts.normalize(pseudocounts=nucleotide_count)
				print('PWM for all ROIs' + str(i) + ' bp long in lib ' + lib)
				print(pwm)
				print(pwm.consensus)
				print(pwm.degenerate_consensus)
				if len(instances) > 20:
					m_high = motifs.create(instances[:round(len(instances) * 0.05)])
					pwm_high = m_high.counts.normalize(pseudocounts=nucleotide_count)
					print('PWM for high expressed ROIs' + str(i) + ' bp long in lib ' + lib)
					print(pwm_high)
					print(pwm_high.consensus)
					print(pwm_high.degenerate_consensus)
					m_low = motifs.create(instances[round(len(instances) * 0.95):])
					pwm_low = m_low.counts.normalize(pseudocounts=nucleotide_count)
					print('PWM for low expressed ROIs' + str(i) + ' bp long in lib ' + lib)
					print(pwm_low)
					print(pwm_low.consensus)
					print(pwm_low.degenerate_consensus)
	with open('run_info.txt', 'a') as info_file:
		info_file.writelines(lines_of_text)
			
			
make_pwm(mut_dictionary)


def find_k_mer(mut_dict, len_k_mer):
	k_mer_dict = {}
	k_mer_count = 0
	for mut_tuple in mut_dict:
		ROI = mut_tuple[0]
		for i in range(len(ROI) - len_k_mer + 1):
			k_mer = ROI[i:i + len_k_mer]
			k_mer_count += 1
			if k_mer not in k_mer_dict:
				k_mer_dict[k_mer] = 1
			else:
				k_mer_dict[k_mer] += 1
	for k_mer in k_mer_dict:
		k_mer_dict[k_mer] = k_mer_dict[k_mer] / k_mer_count * 100
	return k_mer_dict


def count_k_mer(ROI_dict):
	lines_of_text = []
	for k_mer_len in range(4, 8):
		for lib in ROI_dict:
			signif_k_mer = {
				'k_mer': [], 'Total count': [], 'High count': [], 'Low count': []
				}
			k_mer_diction = find_k_mer(ROI_dict[lib], k_mer_len)
			lines_of_text.append(
					str(len(k_mer_diction)) + ' k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
					)
			k_mer_high_diction = find_k_mer(
					ROI_dict[lib][:round(len(ROI_dict[lib]) * 0.05)], k_mer_len
					)
			lines_of_text.append(
					str(len(k_mer_high_diction)) + ' high k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
					)
			k_mer_low_diction = find_k_mer(
					ROI_dict[lib][round(len(ROI_dict[lib]) * 0.95):], k_mer_len
					)
			lines_of_text.append(
					str(len(k_mer_low_diction)) + ' low k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
					)
			print(lines_of_text[len(lines_of_text) - 3:])
			for high_k_mer in k_mer_high_diction:
				if k_mer_high_diction[high_k_mer] / (k_mer_len - 1) >= k_mer_diction[high_k_mer]:
					if high_k_mer not in k_mer_low_diction:
						signif_k_mer['k_mer'].append(high_k_mer)
						signif_k_mer['Total count'].append(
								k_mer_diction[high_k_mer]
								)
						signif_k_mer['High count'].append(
								k_mer_high_diction[high_k_mer]
								)
						signif_k_mer['Low count'].append(0)
					elif k_mer_low_diction[high_k_mer] * (k_mer_len - 1) <= k_mer_diction[high_k_mer]:
						signif_k_mer['k_mer'].append(high_k_mer)
						signif_k_mer['Total count'].append(
								k_mer_diction[high_k_mer]
								)
						signif_k_mer['High count'].append(
								k_mer_high_diction[high_k_mer]
								)
						signif_k_mer['Low count'].append(
								k_mer_low_diction[high_k_mer]
								)
			for low_k_mer in k_mer_low_diction:
				if k_mer_low_diction[low_k_mer] / (k_mer_len - 1) >= k_mer_diction[low_k_mer]:
					if low_k_mer not in k_mer_high_diction:
						signif_k_mer['k_mer'].append(low_k_mer)
						signif_k_mer['Total count'].append(
								k_mer_diction[low_k_mer]
								)
						signif_k_mer['High count'].append(0)
						signif_k_mer['Low count'].append(
								k_mer_low_diction[low_k_mer]
								)
					elif k_mer_high_diction[low_k_mer] * (k_mer_len - 1) <= k_mer_diction[low_k_mer]:
						signif_k_mer['k_mer'].append(low_k_mer)
						signif_k_mer['Total count'].append(
								k_mer_diction[low_k_mer]
								)
						signif_k_mer['High count'].append(
								k_mer_high_diction[low_k_mer]
								)
						signif_k_mer['Low count'].append(
								k_mer_low_diction[low_k_mer]
								)
			df = pd.DataFrame(data=signif_k_mer)
			df['High count / Total count'] = df['High count'] / df['Total count']
			df['Low count / Total count'] = df['Low count'] / df['Total count']
			df = df.sort_values(
								by=['High count / Total count'], ascending=False
								)
			table_path = lib + '_k-mer_' + str(k_mer_len) + '.xlsx'
			table = df.to_excel(table_path, index = None)


count_k_mer(mut_dictionary)


def find_log2_k_mer(ROI_lst, len_k_mer):
	k_mer_dict = {}
	for mut_tuple in ROI_lst:
		ROI = mut_tuple[0]
		for i in range(len(ROI) - len_k_mer + 1):
			k_mer = ROI[i:i + len_k_mer]
			if k_mer not in k_mer_dict:
				k_mer_dict[k_mer] = [mut_tuple[1]]
			else:
				k_mer_dict[k_mer].append(mut_tuple[1])
	return k_mer_dict


def log2_k_mer_stat_test(mut_log2_dict, wt_mut_log2, log_2_wt_gaus):
	for k_mer_len in range(4, 7):
		for lib in mut_log2_dict:
			signif_k_mer = {'k_mer': [], 'Mean expression': []}
			k_mer_diction = find_log2_k_mer(mut_log2_dict[lib], k_mer_len)
			for k_mer in k_mer_diction:
				if len(k_mer_diction[k_mer]) >= 20:
					stat, p = stats.shapiro(k_mer_diction[k_mer])
					if p > 0.05 and log_2_wt_gaus[lib]:
						statis, p_value = stats.ttest_ind(
								wt_mut_log2[lib], k_mer_diction[k_mer]
								)
					else:
						statis, p_value = stats.mannwhitneyu(
								wt_mut_log2[lib], k_mer_diction[k_mer]
								)
					if p_value <= 0.05:
						signif_k_mer['k_mer'].append(k_mer)
						k_mer_mean = sum(k_mer_diction[k_mer]) / float(len(k_mer_diction[k_mer]))
						signif_k_mer['Mean expression'].append(k_mer_mean)
			df = pd.DataFrame(data=signif_k_mer)
			df = df.sort_values(by=['Mean expression'], ascending=False)
			table_path = lib + '_stat_k-mer_' + str(k_mer_len) + '.xlsx'
			table = df.to_excel(table_path, index = None)


log2_k_mer_stat_test(mut_log2_dictionary, wt_log2, log_2_wt_gaussian)
