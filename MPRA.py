# -*- coding: utf-8 -*-
"""
MPFA Anya
"""


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


#Sequences of all indexes, constant parts and controls.
barcode_length = 18
control_barcode_length = 20
mutation_length = 8
index_dict = {
        'A1': 'TTCGGAGT', 'A2': 'ACTCATTT', 'A3': 'GGGATCCG', 'A4': 'TCAAGCAA',
        'A5': 'CAAGATAA', 'A6': 'GGACAACG', 'A7': 'AGCGAGCT', 'A8': 'CTGCACGT',
        'A9': 'GCACTAGT', 'A10': 'AGTCGCCG', 'A11': 'TAAACATC',
        'A12': 'ACAATTCG', 'A13': 'TACTTGTC', 'A14': 'GTACCGTT',
        'A15': 'GCCACATA', 'A16': 'CCTATGGT', 'A17': 'AACGTCGC',
        'A18': 'AGGCAGCA', 'A19': 'AGCTTTCT', 'A20': 'GGTATGTT',
        'A21': 'GAGGGACC', 'A22': 'TAGCTCTA', 'A23': 'TAATTGCG',
        'A24': 'GAAATGGG', 'A25': 'TCGAGACT', 'A26': 'CAGAGAGG',
        'A27': 'ATTAGTCA', 'R701': 'ATCACG'
        }
mapping_indexes = [
        'A1', 'A2',  'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'R701'
        ]
norm_or_expr_indexes = [
        'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19',
        'A20', 'A21', 'A22', 'A23', 'A24', 'A25', 'A26', 'A27'
        ]
library_list = [
        '17-24', '21-28', '25-32', '29-36', '33-40', '37-44', '41-48', '45-52',
        '49-56', 'H'
        ]
replicate_type_list = ['normalization', 'expression', 'mapping']

expr_or_norm_const_part1 = 'CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTC'
expr_or_norm_const_part2 = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGC'
expr_or_norm_control_const_part2 = 'CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCG'

map_const_part1 = 'GACACTCGAGGATCGAG'
map_const_part2 = {
        '17-24':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGT',
        '21-28':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATG',
        '25-32':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTC',
        '29-36':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATT',
        '33-40':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTAC',
        '37-44':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCAT',
        '41-48':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCATGATT',
        '45-52':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCATGATTATCT',
        '49-56':
            'GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTCAATTTTACGCATGATTATCTTTAA',
        'hybrid_test':
            'GAATTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTCATGCGTC'
            }
wt_mutation = {
        '17-24': 'CATGCGTC', '21-28': 'CGTCAATT', '25-32': 'AATTTTAC',
        '29-36': 'TTACGCAT', '33-40': 'GCATGATT', '37-44': 'GATTATCT',
        '41-48': 'ATCTTTAA', '45-52': 'TTAACGTA', '49-56': 'CGTACGTC'
        }
map_const_part3 = {
        '17-24': 'AATTTTACGCATGATTATCTTTAACGTACGTCACAAT',
        '21-28': 'TTACGCATGATTATCTTTAACGTACGTCACAAT',
        '25-32': 'GCATGATTATCTTTAACGTACGTCACAAT',
        '29-36': 'GATTATCTTTAACGTACGTCACAAT',
        '33-40': 'ATCTTTAACGTACGTCACAAT', '37-44': 'TTAACGTACGTCACAAT',
        '41-48': 'CGTACGTCACAAT', '45-52': 'CGTCACAAT', '49-56': 'ACAAT',
        'control': 'GATTATCTTTAACGTACGTCACA'
        }

control_barcodes = {
        'forward': {
                'bc_1': 'TTCCAAGTGCAGGTTAGGCG', 'bc_2': 'TGTGTACGGCTTGCTCTCAA',
                'bc_3': 'GAGCCCGGATCCACTCCAAG', 'bc_4': 'TGTCACGTCAGCTAACCCAC'
                },
        'reversed':{
                'bc_1': 'CGCCTAACCTGCACTTGGAA', 'bc_2': 'TTGAGAGCAAGCCGTACACA',
                'bc_3': 'CTTGGAGTGGATCCGGGCTC', 'bc_4': 'GTGGGTTAGCTGACGTGACA'
                }
             }

control_mutations = {'wt': 'TTACGCAT', 'del_c': 'TTAGCATG'}

hybrid_test_barcodes = {
        'bc15': 'ACTAGGCAAGGCACCAGT', 'bc21': 'GGATATAGAGCGAGTTAA'
        }
hybrid_test_mutations = {'mut15': 'CGACTTAT', 'mut21': 'GCCATATA'}


#Default mismatches levels, lengths of last constant parts, minimum read counts.
index_mismatche_8bp = 1.1
index_mismatche_R701 = 0
expr_norm_const1_mismatche = 4.4
expr_norm_const2_length = 20
expr_norm_const2_mismatche = 2.2

map_const1_mismatche = 2.2
map_const2_mismatche = {
        '17-24': 6.6, '21-28': 7.7, '25-32': 7.7, '29-36': 8.8, '33-40': 8.8,
        '37-44': 8.8, '41-48': 9.9, '45-52': 9.9, '49-56': 10.9
        }
map_const3_length = {
        '17-24': 10, '21-28': 10, '25-32': 10, '29-36': 10, '33-40': 10,
        '37-44': 10, '41-48': 10, '45-52': 9, '49-56': 5
        }
map_const3_mismatche = {
        '17-24': 1.1, '21-28': 1.1, '25-32': 1.1, '29-36': 1.1, '33-40': 1.1,
        '37-44': 1.1, '41-48': 1.1, '45-52': 1.1, '49-56': 0
        }
barcode_mismatche = 2.2
mutation_mismatche = 1.1

minimum_read_count = {'normalization': 3, 'expression': 0}
major_mutation_portion = 0.9


#Some functions to collect experiment information.
#Ask user how many files or libraries will be analyzed.
def input_count(name, additional):
    input_number = input(
            'How many %s%s would you like to analyze? Press number: '
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


#Check if the path to file or library name exist.
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
            print('Wrong file format. File name should ends with .fastq.gz.\nPlease, try again.')
            return check_name(n, name_2)


    elif name_2 == 'library name or H for hybrid test':
        if next_name.upper() in library_list:
            return next_name.upper()
        else:
            print('No such library: %s. The library name should be in this format: 29-36.\nPlease, try again.'
                    % next_name
                    )
            return check_name(n, name_2)


#Validate index.
def check_index(
        rep_type, index_name_list, file_dictionary, name_file, libr_dictionary,
        libr_name
        ):
    valid_indexes = []
    for index_name in index_name_list:
        if index_name in index_dict:
            if (index_name in mapping_indexes and rep_type == 'mapping') or (index_name in norm_or_expr_indexes and rep_type == 'expression') or (index_name in norm_or_expr_indexes and rep_type == 'normalization'):
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


#Continue or restart step.
def yes_or_no():
    answer = input(
            'Is everything correct? Press Y to continue, press N to restart previous step. '
            ).upper()
    if answer == 'Y' or answer == 'N':
        return answer
    else:
        print('You entered an incorrect value: %s. Try again.' % answer)
        return yes_or_no()


#Form dictionary of input files.
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


#Add path to output file, make a text file.
def add_out_path():
    out_path = input(
            'Enter absolute path to your output file (for example D:\\test-folder): '
            )
    new_results = out_path + '\\results'
    run_info = new_results + '\\run_info.txt'
    try:
        os.makedirs(new_results)
        with open(run_info, 'w') as file:
            file.write('Hello, user!\n')
        return new_results, run_info
    except OSError:
        print('Creation of the directory %s failed' % out_path)
        return add_out_path()


#Check presence of normalization, expression or mapping replicates. 
def have_any(rep_name, lib_name, file_name):
    answ = input(
            'Do you have any %s replicates for library %s in file %s? Press Y or N. '
            %(rep_name, lib_name, file_name)
            )
    if answ.upper() == 'Y' or answ.upper() == 'N':
        return answ.upper()
    else:
        print('You entered an incorrect value: %s. Try again.' % answ)
        return have_any(rep_name, lib_name, file_name)


#Form dictionary of libraries in format: {
#    '29-36': {'normalization': {1.fastq.gz: [[A11, A13], [A15]]},
#              'expression': {1.fastq.gz: [[A12], [A14]]},
#              'mapping': {1.fastq.gz: [[A1]]},
#             }
#     'H' : {1.fastq.gz: [A2, A3, A4]}
#    }
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

                for replicate_type in replicate_type_list:
                    response = have_any(replicate_type, library_name, file)
                    if response == 'Y':
                        library_dict[library_name][replicate_type][file] = []
                        print('Enter indexes for %s replicates of library %s in file %s.\nIf there are more then one index per replicate, print indexes in one line separated by spaces. Finally press N.'
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
                print('Enter indexes for hybrid test in file %s.\nIf there are more then one index per replicate, print indexes in one line separated by spaces. Finally press N.'
                        % file
                        )
                if library_name not in library_dict:
                    library_dict['H'] = {file: []}
                else:
                    library_dict['H'][file] = []

                ind = input().upper().strip()
                while ind != 'N':
                    if ind in mapping_indexes:
                        if ind not in file_diction[file]:
                            file_diction[file][ind] = {'read count': 0}
                            library_dict['H'][file].append(ind)
                            print('Index %s was successfully added.' % ind)
                        else:
                            print('Index %s was already used in file %s.'
                                  %(ind, file)
                                  )
                    else:
                        print('Index %s is not appropriate for hybrid test.'
                              % ind
                              )
                    ind = input().upper().strip()
    
    for lib_name in library_dict:
        if library_dict[library_name] == {
                'normalization': {}, 'expression': {}, 'mapping': {}
                }:
            raise Exception('Error! Indexes for library %s was not added!'
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
                print('Hybrid test indexes {0} in file {1}.'.format(
                        library_dict[library][file], file
                        ))
                    
    resp = yes_or_no()
    if resp == 'Y':
        return file_diction, library_dict
    else:
        for file in file_diction:
            file_diction[file] = {}
        return add_lib(name_plur, file_diction, name_sing)


#Asks the user if he wants to change the level of mismatches.
def change_or_not():
    ans = input('Would you like to change something? Press Y or N: ').upper()
    if ans == 'Y' or ans == 'N':
        return ans
    else:
        print('You entered an incorrect value: %s. Try again.' % ans)
        return change_or_not()


#Check if the current_value lower or equal the max_value.
def mismatche_max(input_message, max_value, libr_name=''):
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
        return mismatche_max(input_message, max_value, libr_name='')


#Chenges values in dictionary minimum_read_count
def change_read_count(rep):
    try:
        min_read_count = int(input('Enter minimum read count for ' + rep + ': '))
        if rep == 'normalization' and min_read_count <= 0:
            print('Minimum read count for normalization can\'t be 0. Please, try again.')
            change_read_count(rep)
        return min_read_count
    except ValueError:
        print('You entered incorrect value. Try again.')
        change_read_count(rep)


#Change the level of mismatches.
def change_mismatches(libr_dict):
    global index_mismatche_8bp, index_mismatche_R701, barcode_mismatche, mutation_mismatche, expr_norm_const1_mismatche, expr_norm_const2_mismatche, expr_norm_const2_length, map_const1_mismatche, map_const2_mismatche, map_const3_length, map_const3_mismatche, minimum_read_count
    
    index_mismatche_8bp = mismatche_max(
            'mismatche value for indexes A1-A27', 4
            )
    index_mismatche_R701 =  mismatche_max('mismatche value for index R701', 3)
    barcode_mismatche = mismatche_max(
            'mismatche value for barcode', barcode_length // 2 - 1
            )
    mutation_mismatche = mismatche_max(
            'mismatche value for mutation', mutation_length // 2 - 1
            )
    expr_norm_const1_mismatche = mismatche_max(
            'mismatche value for normalization or expression constant part 1',
            len(expr_or_norm_const_part1) // 2
            )
    expr_norm_const2_length = mismatche_max(
            'length of normalization or expression constant part 2',
            len(expr_or_norm_control_const_part2)
            )
    expr_norm_const2_mismatche = mismatche_max(
            'mismatche value for normalization or expression constant part 2',
            expr_norm_const2_length // 2
            )
    map_const1_mismatche = mismatche_max(
            'mismatche value for mapping constant part 1',
            len(map_const_part1) // 2
            )
        
    for libr in libr_dict:
        if libr == 'H':
            map_const2_mismatche['25-32'] = mismatche_max(
                    'mismatche value for hybrid test constant part 2',
                    len(map_const_part2['25-32']) // 2
                    )
            map_const3_length['25-32'] = mismatche_max(
                    'length of constant part 3 for hybrid test',
                    len(map_const_part3['25-32'])
                    )
            map_const3_mismatche['25-32'] = mismatche_max(
                    'mismatche level for hybrid test constant part 3',
                    map_const3_length['25-32'] // 2
                    )
        else:
            if  libr_dict[libr]['mapping']:
                map_const2_mismatche[libr] = mismatche_max(
                        'mismatche value for mapping constant part 2 in library ',
                        len(map_const_part2[libr]) // 2, libr
                        )
                
                map_const3_length[libr] = mismatche_max(
                        'length of constant part 3 in library ',
                        len(map_const_part3[libr]), libr
                        )
                map_const3_mismatche[libr] = mismatche_max(
                        'mismatche value for mapping constant part 3 in library ',
                        map_const3_length[libr] // 2, libr
                        )                   

    map_const2_mismatche['29-36'] = mismatche_max(
            'mismatche level for mutation control constant part 2',
            len(map_const_part2['29-36']) // 2
            )
    map_const3_length['29-36'] = mismatche_max(
            'length of constant part 3 for mutation control',
            len(map_const_part3['29-36'])
            )
    map_const3_mismatche['29-36'] = mismatche_max(
            'mismatche level for mutation control constant part 3',
            map_const3_length['29-36'] // 2
            )
    
    for rep_type in minimum_read_count:
        minimum_read_count[rep_type] = change_read_count(rep_type)
    
    report = yes_or_no()
    if report == 'N':
        change_mismatches(libr_dict)


#Display and change mismathes levels
def mismathes_count(libr_dict, run_inf_text):
    print('Let\'s check acceptable mismatche levels.\n{0} mismathe in indexes A1-A27, {1} mismathes in index R701.'.format(
            str(index_mismatche_8bp), str(index_mismatche_R701)
            ))
    print('{0} mismathes in barcode, {1} in mutation.'.format(
            str(barcode_mismatche), str(mutation_mismatche)
            ))
    print('{0} mismathes in normalization or expression constant part 1;\n{1} mismathes in normalization or expression constant part 2 {2} bp long.'.format(
            str(expr_norm_const1_mismatche), str(expr_norm_const2_mismatche),
            str(expr_norm_const2_length)
            ))
    print('{0} mismathes in mapping or hybrid test constant part 1'.format(
            str(map_const1_mismatche)
            ))

    for lib in libr_dict:
        if lib == 'H':
            print('For hybrid test {0} mismathes for mapping constant part 2;\n{1} mismathes for hybrid test constant part 3 {2} bp long.'.format(
                    str(map_const2_mismatche['25-32']),
                    str(map_const3_mismatche['25-32']),
                    str(map_const3_length['25-32'])
                    ))
        elif libr_dict[lib]['mapping']:
            print('For library {0} {1} mismathes in mapping constant part 2;\n{2} mismathes for mapping constant part 3 {3} bp long.'.format(
                    lib, str(map_const2_mismatche[lib]),
                    str(map_const3_mismatche[lib]), str(map_const3_length[lib])
                    ))

    print('For mapping control {0} mismathes for mapping constant part 2;\n{1} mismathes for mapping control constant part 3 {2} bp long.'.format(
            str(map_const2_mismatche['29-36']),
            str(map_const3_mismatche['29-36']), str(map_const3_length['29-36'])
            ))

    for rep_type in minimum_read_count:
        print('Minimum read count for ' + rep_type + ' is ' +
              str(minimum_read_count[rep_type]) + '.')

    respons = change_or_not()

    if respons == 'Y':
        change_mismatches(libr_dict)
        
    lines_of_text = [
            '{0} mismathe in indexes A1-A27, {1} mismathes in index R701.\n'.format(str(index_mismatche_8bp), str(index_mismatche_R701)),
            '{0} mismathes in barcode, {1} in mutation.\n'.format(str(barcode_mismatche), str(mutation_mismatche)),
            '{0} mismathes in normalization or expression constant part 1;\n{1} mismathes in normalization or expression constant part 2 {2} bp long.\n'.format(str(expr_norm_const1_mismatche), str(expr_norm_const2_mismatche), str(expr_norm_const2_length)),
            '{0} mismathes in mapping or hybrid test constant part 1\n'.format(str(map_const1_mismatche)),
            'For mapping control {0} mismathes for mapping constant part 2;\n{1} mismathes for mapping control constant part 3 {2} bp long.\n'.format(str(map_const2_mismatche['29-36']), str(map_const3_mismatche['29-36']), str(map_const3_length['29-36']))
            ]
    for lib in libr_dict:
        if lib == 'H':
            lines_of_text.append(
                    'For hybrid test {0} mismathes for mapping constant part 2;\n{1} mismathes for hybrid test constant part 3 {2} bp long.\n'.format(str(map_const2_mismatche['25-32']), str(map_const3_mismatche['25-32']), str(map_const3_length['25-32']))
                    )
        elif libr_dict[lib]['mapping']:
            lines_of_text.append(
                    'For library {0} {1} mismathes in mapping constant part 2;\n{2} mismathes for mapping constant part 3 {3} bp long.\n'.format(lib, str(map_const2_mismatche[lib]), str(map_const3_mismatche[lib]), str(map_const3_length[lib]))
                    )
    with open(run_inf_text, 'a') as info_file:
        info_file.writelines(lines_of_text)


#Greeting user and collect experiment information.
def collect_info():
    
    print('Hello, my dear friend!')
    
    fastq_file_dict = add_file('FASTQ files', '', 'absolute path to file')
    
    output_path, output_file = add_out_path()
    
    fastq_file_dict, lib_dict = add_lib(
            'libraries', fastq_file_dict, 'library name or H for hybrid test'
            )

    mismathes_count(lib_dict, output_file)
    
    with open(output_path + '\\' + 'fastq_file_dictionary.json', 'w') as f:
        json.dump(fastq_file_dict, f)
        
    with open(output_path + '\\' + 'lib_diction.json', 'w') as file:
        json.dump(lib_dict, file)

    return fastq_file_dict, lib_dict, output_path, output_file


fastq_file_dictionary, lib_diction, output_folder, run_info_text = collect_info()


#Align two sequences
def count_score(query, ref):
    return pairwise2.align.globalxs(query, ref, -0.1, -0.1, score_only=True)


#Find the best alignment
def aligner(read, reference, mismatche):
    insertion = 0
    deletion = 0
    max_score = count_score(read[:len(reference)], reference)
    best_align_length = len(reference)
    if mismatche > 0:
        while insertion < mismatche:
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
            while deletion < mismatche:
                deletion += 1
                deletion_score = count_score(
                        read[:len(reference) - deletion], reference
                        )
                if deletion_score >= max_score:
                    max_score = deletion_score
                    best_align_length = len(reference) - deletion
                else:
                    break
    minimum_score = len(reference) - mismatche + insertion
    return max_score, best_align_length, minimum_score


#Check if read quality upper than 20.
def check_quality(seq, start, end):
    quality_list = seq.letter_annotations['phred_quality']
    for quality in range(start, end):
        if quality_list[quality] < 10:
            return False
    return True


#Find barcodes
def find_barcode(
        rec, start, const1, const1_mismatche, const2, const2_mismatche, mode
        ):
    score1, const1_length, min_score1 = aligner(
            rec.seq[start:], const1, const1_mismatche
            )
    if score1 >= min_score1:
        for control_barcode in control_barcodes[mode]:
            score_conbc, conbc_length, min_score_conbc = aligner(
                    rec.seq[(start + const1_length):],
                    control_barcodes[mode][control_barcode], barcode_mismatche
                    )
            if score_conbc >= min_score_conbc:
                break
        else:
            best_score2, const2_length, best_min_score2 = aligner(
                rec.seq[(start + const1_length + barcode_length):],
                const2, const2_mismatche
                )
            best_barcode_length = barcode_length
            
            insertion = 0
            deletion = 0
            while True:
                insertion += 1
                ins_score2, const2_length, min_score2 = aligner(
                        rec.seq[(start + const1_length + barcode_length + insertion):],
                        const2, const2_mismatche
                        )
                if ins_score2 > best_score2:
                    best_score2 = ins_score2
                    best_min_score2 = min_score2
                    best_barcode_length = barcode_length + insertion
                else:
                    break
                    
            if best_barcode_length == barcode_length:
                while True:
                    deletion += 1
                    del_score2, const2_length, min_score2 = aligner(
                            rec.seq[(start + const1_length + barcode_length - deletion):],
                            const2, const2_mismatche
                            )
                    if del_score2 > best_score2:
                        best_score2 = del_score2
                        best_min_score2 = min_score2
                        best_barcode_length = barcode_length - deletion
                    else:
                        break
            if best_barcode_length <= barcode_length + barcode_mismatche and best_barcode_length >= barcode_length - barcode_mismatche:
                if best_score2 >= best_min_score2:
                    if check_quality(rec, (start + const1_length), (start + const1_length + best_barcode_length)):
                        return rec.seq[start + const1_length:start + const1_length + best_barcode_length]
                    return 'low quality'
        if check_quality(rec, (start + const1_length), (start + const1_length + conbc_length)):
            score2, const2_length, min_score2 = aligner(
                    rec.seq[(start + const1_length + conbc_length):],
                    const2, const2_mismatche
                    )
            if score2 >= min_score2:
                return control_barcode
        return 'low quality'
    return 'undetermined'
    

#Find barcodes and mutations for hybrid test    
def hybrid_test_analyses(rec, start):
    seq = rec.seq
    score1, const1_length, min_score1 = aligner(
            seq[start:], map_const_part1, map_const1_mismatche
            )
        
    if score1 >= min_score1:
        for hybrid_test_barcode in hybrid_test_barcodes:
            score_hbc, hbc_length, min_score_hbc = aligner(
                    seq[start + const1_length:],
                    hybrid_test_barcodes[hybrid_test_barcode],
                    barcode_mismatche
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
                        map_const_part2['hybrid_test'],
                        map_const2_mismatche['25-32']
                        )
            else:
                score2, const2_length, min_score2 = aligner(
                        seq[start + const1_length + hbc_length:],
                        map_const_part2['25-32'],
                        map_const2_mismatche['25-32']
                        )
            if score2 >= min_score2:
                for hybrid_test_mutation in hybrid_test_mutations:
                    score_hmut, hmut_length, min_score_hmut = aligner(
                            seq[start + const1_length + hbc_length + const2_length:],
                            hybrid_test_mutations[hybrid_test_mutation],
                            mutation_mismatche
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
                            map_const_part3['25-32'][:map_const3_length['25-32']],
                            map_const3_mismatche['25-32']
                            )
                    if score3 >= min_score3:
                        return hybrid_test_barcode, hybrid_test_mutation
                return hybrid_test_barcode, 'low quality'
            return hybrid_test_barcode, 'undetermined'
        return 'low quality', 'undetermined'
    return 'undetermined', 'undetermined'


#Find barcodes and mutations for mapping
def find_barcode_mutation(rec, ind, start, lib_dictionar, file):
    seq = rec.seq
    for lib in lib_dictionar:
        if lib == 'H':
            if file in lib_dictionar['H']:
                if ind in lib_dictionar['H'][file]:
                    barcode, mutation = hybrid_test_analyses(rec, start)
                    return barcode, mutation
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
            seq[start:], map_const_part1, map_const1_mismatche
            )
    if score1 >= min_score1:
        for control_barcode in control_barcodes['forward']:
            score_conbc, conbc_length, min_score_conbc = aligner(
                    seq[(start + const1_length):],
                    control_barcodes['forward'][control_barcode],
                    barcode_mismatche
                    )
            if score_conbc >= min_score_conbc:
                cont_barcode = control_barcode
                break
        else:
            best_score2, best_const2_length, best_min_score2 = aligner(
                    seq[start + const1_length + barcode_length:],
                    map_const_part2[librar], map_const2_mismatche[librar]
                    )
            best_barcode_length = barcode_length
            insertion = 0
            deletion = 0
            while True:
                insertion += 1
                ins_score2, const2_length, min_score2 = aligner(
                        seq[start + const1_length + barcode_length + insertion:],
                        map_const_part2[librar],
                        map_const2_mismatche[librar]
                        )
                if ins_score2 > best_score2:
                    best_score2 = ins_score2
                    best_const2_length = const2_length
                    best_min_score2 = min_score2
                    best_barcode_length = barcode_length + insertion
                else:
                    break
            if best_barcode_length == barcode_length:
                while True:
                    deletion += 1
                    del_score2, const2_length, min_score2 = aligner(
                            seq[(start + const1_length + barcode_length - deletion):],
                            map_const_part2[librar],
                            map_const2_mismatche[librar]
                            )
                    if del_score2 > best_score2:
                        best_score2 = del_score2
                        best_const2_length = const2_length
                        best_min_score2 = min_score2
                        best_barcode_length = barcode_length - deletion
                    else:
                        break
            if best_barcode_length <= barcode_length + barcode_mismatche and best_barcode_length >= barcode_length - barcode_mismatche:
                if best_score2 >= best_min_score2:
                    if check_quality(
                            rec, (start + const1_length),
                            (start + const1_length + best_barcode_length)
                            ):
                        best_score3, const3_length, best_min_score3 = aligner(
                                seq[start + const1_length + best_barcode_length + best_const2_length + mutation_length:],
                                map_const_part3[librar][:map_const3_length[librar]],
                                map_const3_mismatche[librar]
                                )
                        best_mutation_length = mutation_length
                        insertion = 0
                        deletion = 0
                        while True:
                            insertion += 1
                            ins_score3, const3_length, min_score3 = aligner(
                                    seq[start + const1_length + best_barcode_length + best_const2_length + mutation_length + insertion:],
                                    map_const_part3[librar][:map_const3_length[librar]],
                                    map_const3_mismatche[librar]
                                    )
                            if ins_score3 > best_score3:
                                best_score3 = ins_score3
                                best_min_score3 = min_score3
                                best_mutation_length = mutation_length + insertion
                            else:
                                break
                    
                        if best_mutation_length == mutation_length:
                            while True:
                                deletion += 1
                                del_score3, const3_length, min_score3 = aligner(
                                        seq[start + const1_length + best_barcode_length + best_const2_length + mutation_length - deletion:],
                                        map_const_part3[librar][:map_const3_length[librar]],
                                        map_const3_mismatche[librar]
                                        )
                                if del_score3 > best_score3:
                                    best_score3 = del_score3
                                    best_min_score3 = min_score3
                                    best_mutation_length = mutation_length - deletion
                                else:
                                    break
                        if best_mutation_length >= mutation_length - mutation_mismatche and best_mutation_length <= mutation_length + mutation_mismatche:
                            if best_score3 >= best_min_score3:
                                if check_quality(
                                        rec,
                                        (start + const1_length + best_barcode_length + best_const2_length),
                                        (start + const1_length + best_barcode_length + best_const2_length + best_mutation_length)
                                        ):
                                    return seq[
                                            start + const1_length:start + const1_length + best_barcode_length
                                            ], seq[start + const1_length + best_barcode_length + best_const2_length:start + const1_length + best_barcode_length + best_const2_length + best_mutation_length]
                                return 'low quality mutation', 'low quality mutation'
                        return 'undetermined', 'undetermined'
                    return 'low quality barcode', 'low quality barcode'
            return 'undetermined', 'undetermined'
        if check_quality(
                rec, (start + const1_length),
                (start + const1_length + conbc_length)
                ):
            score2, const2_length, min_score2 = aligner(
                    seq[(start + const1_length + conbc_length):],
                    map_const_part2['29-36'],
                    map_const2_mismatche['29-36']
                    )
            if score2 >= min_score2:
                if cont_barcode == 'bc_1' or cont_barcode == 'bc_2':
                    score_cmut, cmut_length, min_score_cmut = aligner(
                            seq[start + const1_length + conbc_length + const2_length:],
                            control_mutations['wt'],
                            mutation_mismatche
                            )
                    cont_mutation = 'wt'
                else:
                    score_cmut, cmut_length, min_score_cmut = aligner(
                            seq[start + const1_length + conbc_length + const2_length:],
                            control_mutations['del_c'],
                            mutation_mismatche
                            )
                    cont_mutation = 'del_c'
                if score_cmut >= min_score_cmut:
                    if check_quality(
                            rec,
                            (start + const1_length + conbc_length + const2_length),
                            (start + const1_length + conbc_length + const2_length + cmut_length)):
                        score3, const3_length, min_score3 = aligner(
                                seq[start + const1_length + conbc_length + const2_length + cmut_length:],
                                map_const_part3['29-36'][:map_const3_length['29-36']],
                                map_const3_mismatche['29-36']
                                )
                        if score3 >= min_score3:
                            return cont_barcode, cont_mutation
                    return 'low quality mutation', 'low quality mutation'
            return 'undetermined', 'undetermined'
        return 'low quality barcode', 'low quality barcode'                    
    return 'undetermined', 'undetermined'
        

#Form dictionary of barcodes and mutations in format:
#{'1.fastq.gz':
#              {A11: {'read count': 1}, {'eff read count': 1}, {'bc1': ['seq1']}},
#              {A13: {'read count': 1}, {'eff read count': 1}, {'bc2': ['seq2']}},
#              {A15: {'read count': 0}} {'eff read count': 0},,
#              {A12: {'read count': 1}, {'eff read count': 1}, {'bc3': ['seq3']}},
#              {A14: {'read count': 1}, {'eff read count': 1}, {'bc4': ['seq4']}},
#              {A1: {'read count': 1}, {'eff read count': 1}, {'bc4': {'mut4': ['seq8']}}},
#              {A2: {'read count': 1}, {'eff read count': 1}, {'bc1': {'mut1': ['seq5']}}},
#              {A3: {'read count': 1}, {'eff read count': 1}, {'bc2': {'mut2': ['seq6']}}},
#              {A4: {'read count': 1}, {'eff read count': 1}, {'bc3': {'mut3': ['seq7']}}}
#}
def fastq_parsing(ind_dictionary, lib_dictionary):
    start_time = datetime.now()
    print('Analyses was started ' + str(start_time))
    for inp_file in ind_dictionary:
        undetermined_indexes = 0
        low_quality_index = 0
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
                                record.seq, index_dict['R701'],
                                index_mismatche_R701
                                )
                    else:
                        score, index_length, min_score = aligner(
                                record.seq, index_dict[ind],
                                index_mismatche_8bp
                                )
                    if score >= min_score:
                        index = ind
                        break
                else:
                    undetermined_indexes += 1
                    continue
                if check_quality(record, 0, index_length):
                    ind_dictionary[inp_file][index]['read count'] += 1
                    if index in norm_or_expr_indexes:
                        barcode = find_barcode(
                                record, index_length,
                                expr_or_norm_const_part1,
                                expr_norm_const1_mismatche,
                                expr_or_norm_const_part2[:expr_norm_const2_length],
                                expr_norm_const2_mismatche, 'reversed'
                                )
                        if type(barcode) is not str:
                            barcode = str(barcode.reverse_complement())
                        if barcode not in ind_dictionary[inp_file][index]:
                            ind_dictionary[inp_file][index][barcode] = [record_sequence]
                        else:
                            ind_dictionary[inp_file][index][barcode].append(record_sequence)
                    else:
                        barcode, mutation = find_barcode_mutation(
                                record, index, index_length,
                                lib_dictionary, inp_file
                                )
                        barcode = str(barcode)
                        mutation = str(mutation)
                        if barcode not in ind_dictionary[inp_file][index]:
                            ind_dictionary[inp_file][index][barcode] = {mutation: [record_sequence]}
                        else:
                            if mutation not in ind_dictionary[inp_file][index][barcode]:
                                ind_dictionary[inp_file][index][barcode][mutation] = [record_sequence]
                            else:
                                ind_dictionary[inp_file][index][barcode][mutation].append(record_sequence)
                else:
                    low_quality_index += 1
        
        lines_of_text = [
                'Total read count in file %s is %s.\n' %(inp_file, str(total_read_count)),
                'In file ' + inp_file + ' was found:\n'
                ]
        for index in ind_dictionary[inp_file]:
            if index in mapping_indexes:
                try:
                    ind_dictionary[inp_file][index]['low quality count'] = len(
                            ind_dictionary[inp_file][index]['low quality barcode']['low quality barcode']
                            ) + len(ind_dictionary[inp_file][index]['low quality mutation']['low quality mutation'])
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
            lines_of_text.append('%s total reads, %s effective reads, %s low quality reads and %s identical barcodes for index %s.\n'
                  %(str(ind_dictionary[inp_file][index]['read count']),
                    str(ind_dictionary[inp_file][index]['eff read count']),
                    str(ind_dictionary[inp_file][index]['low quality count']),
                    str(len(ind_dictionary[inp_file][index]) - delete), index))
        lines_of_text += [
                str(undetermined_indexes) + ' reads without any index.\n',
                str(low_quality_index) + ' reads with low quality index.\n'
                ]
        
        print(*lines_of_text)
        
        with open(run_info_text, 'a') as info_file:
            info_file.writelines(lines_of_text)
    
    print('Analyses was finished ' + str(datetime.now()))
    with open(output_folder + '\\' + 'fastq_file_dictionary.json', 'w') as f:
        json.dump(ind_dictionary, f)
    return ind_dictionary    


fastq_file_dictionary = fastq_parsing(fastq_file_dictionary, lib_diction)

#print(lib_diction)
#n = 0
#for read in fastq_file_dictionary['C:/Users/user/1_S1_L001_R1_001.fastq.gz']['R701']['bc15']['undetermined']:
#    print(read)
#    n += 1
#    if n > 50:
#        break
#for barcode in fastq_file_dictionary['C:\\Users\\user\\1_S1_L001_R1_001_2018_10_05.fastq.gz']['A1']:
#    print(barcode)
#print(fastq_file_dictionary['C:\\Users\\user\\1_S1_L001_R1_001_2018_11_20.fastq.gz']['A1']['CCCGGAGAGGGGGGGAAG']['CATCTTAT'])

#print(len(fastq_file_dictionary['C:/Users/user/1_S1_L001_R1_001.fastq.gz']['A4']['bc21']['undetermined']))


#Form dictionary of barcodes-mutations, found in at least two replicates.
#gold_mapping = {
#                'barcode1': {'mut1': mut1_count, 'mut2': mut2_count},
#                'barcode2': {'mut3': mut3_count, 'mut4': mut4_count}
#               }
#mapping_count = {'file name': {'A1': {
#                                      'barcode1': {
#                                                   'mut1': mut1_count1,
#                                                   'mut2': mut2_count1
#                                                  },
#                                      'barcode2': {
#                                                   'mut3': mut3_count1,
#                                                   'mut4': mut4_count1
#                                                   }
#                                     },
#                               'A2': {
#                                      'barcode1': {
#                                                   'mut1': mut1_count2,
#                                                   'mut2': mut2_count2
#                                                   },
#                                      'barcode2': {
#                                                   'mut3': mut3_count2,
#                                                   'mut4': mut4_count2
#                                                  }
#                                     }}}
def form_gold_mapping(ind_dictionary, lib_dictionary, lib):
    gold_mapping = {}
    #This includes barcodes that are not found in other replicates.
    single_mapping = {}
    #This is temporal dictionary. Includes reads, foubd only in this replicate.
    rep_mapping = {}
    #This dectionary includes all barcodes, mutations and mutation counts for each replicate.
    mapping_count = {}
    control_barcodes_dict = {}
    if len(lib_dictionary[lib]['mapping']) == 1:
        for file in lib_dictionary[lib]['mapping']:
            mapping_count[file] = {}
            control_barcodes_dict[file] = {}
            if len(lib_dictionary[lib]['mapping'][file]) == 1:
                for index_list in lib_dictionary[lib]['mapping'][file]:
                    control_barcodes_dict[file][', '.join(index_list)] = {}
                    for index in index_list:
                        for barcode in ind_dictionary[file][index]:
                            #Check that it is not 'read count' or control barcode.
                            if len(barcode) > 15 and barcode not in [
                                    'low quality count',
                                    'low quality barcode',
                                    'low quality mutation'
                                    ]:
                                for mutation in ind_dictionary[file][index][barcode]:
                                    mutation_count = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                    if barcode in gold_mapping:
                                        if mutation in gold_mapping[barcode]:
                                            gold_mapping[barcode][mutation] += mutation_count
                                        else:
                                            gold_mapping[barcode][mutation] = mutation_count
                                    else:
                                        gold_mapping[barcode] = {
                                                mutation: mutation_count
                                                }
                            elif barcode in control_barcodes['forward']:
                                if barcode == 'bc_1' or barcode == 'bc_2':
                                    mutation = 'wt'
                                else:
                                    mutation = 'del_c'
                                if barcode in control_barcodes_dict[file][', '.join(index_list)]:
                                    control_barcodes_dict[file][', '.join(index_list)][barcode] += len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                else:
                                    control_barcodes_dict[file][', '.join(index_list)][barcode] = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                    mapping_count[file][', '.join(index_list)] = gold_mapping
            else:    
                for i in range(len(lib_dictionary[lib]['mapping'][file])):
                    rep_name = ', '.join(
                            lib_dictionary[lib]['mapping'][file][i]
                            )
                    mapping_count[file][rep_name] = {}
                    control_barcodes_dict[file][rep_name] = {}
                    for index in lib_dictionary[lib]['mapping'][file][i]:
                        for barcode in ind_dictionary[file][index]:
                            #Check that it is not 'read count' or control barcode.
                            if len(barcode) > 15 and barcode not in [
                                    'low quality count',
                                    'low quality barcode',
                                    'low quality mutation'
                                    ]:
                                for mutation in ind_dictionary[file][index][barcode]:
                                    mutation_count = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                    if barcode in mapping_count[file][rep_name]:
                                        if mutation in mapping_count[file][rep_name][barcode]:
                                            mapping_count[file][rep_name][barcode][mutation] += mutation_count
                                        else:
                                            mapping_count[file][rep_name][barcode][mutation] = mutation_count
                                    else:
                                        mapping_count[file][rep_name][barcode] = {
                                                mutation: mutation_count
                                                }
                                    if barcode in gold_mapping:
                                        if mutation in gold_mapping[barcode]:
                                            gold_mapping[barcode][mutation] += mutation_count
                                        else:
                                            for n in range(i + 1, len(
                                                    lib_dictionary[lib]['mapping'][file]
                                                    )):
                                                for ind in lib_dictionary[lib]['mapping'][file][n]:
                                                    if barcode in ind_dictionary[file][ind]:
                                                        if mutation in ind_dictionary[file][ind][barcode]:
                                                            gold_mapping[barcode][mutation] = mutation_count
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
                                                    if mutation in ind_dictionary[file][ind][barcode]:
                                                        gold_mapping[barcode] = {
                                                                mutation: mutation_count
                                                                }
                                                        break
                                            else:
                                                continue
                                            break
                            elif barcode in control_barcodes['forward']:
                                if barcode == 'bc_1' or barcode == 'bc_2':
                                    mutation = 'wt'
                                else:
                                    mutation = 'del_c'
                                if barcode in control_barcodes_dict[file][rep_name]:
                                    control_barcodes_dict[file][rep_name][barcode] += len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                else:
                                    control_barcodes_dict[file][rep_name][barcode] = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
    else:
        for file in lib_dictionary[lib]['mapping']:
            mapping_count[file] = {}
            control_barcodes_dict[file] = {}
            if len(lib_dictionary[lib]['mapping'][file]) == 1:
                for index_list in lib_dictionary[lib]['mapping'][file]:
#This is temporal dictionary. Includes reads, found only in this replicate.
                    rep_mapping = {}
                    rep_name = ', '.join(index_list)
                    mapping_count[file][rep_name] = {}
                    control_barcodes_dict[file][rep_name] = {}
                    for index in index_list:
                        for barcode in ind_dictionary[file][index]:
                            #Check that it is not 'read count' or control barcode.
                            if len(barcode) > 15 and barcode not in [
                                    'low quality count',
                                    'low quality barcode',
                                    'low quality mutation'
                                    ]:
                                for mutation in ind_dictionary[file][index][barcode]:
                                    mutation_count = len(
                                                ind_dictionary[file][index][barcode][mutation]
                                                ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                    if barcode in mapping_count[file][rep_name]:
                                        if mutation in mapping_count[file][rep_name][barcode]:
                                            mapping_count[file][rep_name][barcode][mutation] += mutation_count
                                        else:
                                            mapping_count[file][rep_name][barcode][mutation] = mutation_count
                                    else:
                                        mapping_count[file][rep_name][barcode] = {
                                                mutation: mutation_count
                                                }
                                    if barcode in gold_mapping:
                                        if mutation in gold_mapping[barcode]:
                                            gold_mapping[barcode][mutation] += mutation_count
                                        else:
                                            if barcode in single_mapping:
                                                if mutation in single_mapping[barcode]:
                                                    gold_mapping[barcode][mutation] = mutation_count
                                                    single_mapping[barcode].pop(
                                                            mutation
                                                            )
                                                else:
                                                    if barcode in rep_mapping:
                                                        if mutation in rep_mapping:
                                                            rep_mapping[barcode][mutation] += mutation_count
                                                        else:
                                                            rep_mapping[barcode][mutation] = mutation_count
                                                    else:
                                                        rep_mapping[barcode] = {
                                                                mutation: mutation_count
                                                                }
                                            else:
                                                if barcode in rep_mapping:
                                                    if mutation in rep_mapping:
                                                        rep_mapping[barcode][mutation] += mutation_count
                                                    else:
                                                        rep_mapping[barcode][mutation] = mutation_count
                                                else:
                                                    rep_mapping[barcode] = {
                                                            mutation: mutation_count
                                                            }
                                    else:
                                        if barcode in single_mapping:
                                            if mutation in single_mapping[barcode]:
                                                gold_mapping[barcode] = {
                                                        mutation: mutation_count
                                                        }
                                                single_mapping[barcode].pop(
                                                        mutation
                                                        )
                                            else:
                                                if barcode in rep_mapping:
                                                    if mutation in rep_mapping:
                                                        rep_mapping[barcode][mutation] += mutation_count
                                                    else:
                                                        rep_mapping[barcode][mutation] = mutation_count
                                                else:
                                                    rep_mapping[barcode] = {
                                                            mutation: mutation_count
                                                            }
                                        else:
                                            if barcode in rep_mapping:
                                                if mutation in rep_mapping:
                                                    rep_mapping[barcode][mutation] += mutation_count
                                                else:
                                                    rep_mapping[barcode][mutation] = mutation_count
                                            else:
                                                rep_mapping[barcode] = {
                                                        mutation: mutation_count
                                                        }
                            elif barcode in control_barcodes['forward']:
                                if barcode == 'bc_1' or barcode == 'bc_2':
                                    mutation = 'wt'
                                else:
                                    mutation = 'del_c'
                                if barcode in control_barcodes_dict[file][rep_name]:
                                    control_barcodes_dict[file][rep_name][barcode] += len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                else:
                                    control_barcodes_dict[file][rep_name][barcode] = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                    for barcode in rep_mapping:
                        for mutation in rep_mapping[barcode]:
                            mutation_count = rep_mapping[barcode][mutation]
                            if barcode in single_mapping:
                                single_mapping[barcode][mutation] = mutation_count
                            else:
                                single_mapping[barcode] = {
                                        mutation: mutation_count
                                        }
            else:
                for i in range(len(lib_dictionary[lib]['mapping'][file])):
#This is temporal dictionary. Includes reads, found only in this replicate.
                    rep_mapping = {}
                    rep_name = ', '.join(
                            lib_dictionary[lib]['mapping'][file][i]
                            )
                    mapping_count[file][rep_name] = {}
                    control_barcodes_dict[file][rep_name] = {}
                    for index in lib_dictionary[lib]['mapping'][file][i]:
                        for barcode in ind_dictionary[file][index]:
                            #Check that it is not 'read count' or control barcode.
                            if len(barcode) > 15 and barcode not in [
                                    'low quality count',
                                    'low quality barcode',
                                    'low quality mutation'
                                    ]:
                                for mutation in ind_dictionary[file][index][barcode]:
                                    mutation_count = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                    if barcode in mapping_count[file][rep_name]:
                                        if mutation in mapping_count[file][rep_name][barcode]:
                                            mapping_count[file][rep_name][barcode][mutation] += mutation_count
                                        else:
                                            mapping_count[file][rep_name][barcode][mutation] = mutation_count
                                    else:
                                        mapping_count[file][rep_name][barcode] = {
                                                mutation: mutation_count
                                                }
                                    if barcode in gold_mapping:
                                        if mutation in gold_mapping[barcode]:
                                            gold_mapping[barcode][mutation] += mutation_count
                                        else:
                                            for n in range(i + 1, len(
                                                    lib_dictionary[lib]['mapping'][file]
                                                    )):
                                                for ind in lib_dictionary[lib]['mapping'][file][n]:
                                                    if barcode in ind_dictionary[file][ind]:
                                                        if mutation in ind_dictionary[file][ind][barcode]:
                                                            gold_mapping[barcode][mutation] = mutation_count
                                                            break
                                                else:
                                                    continue
                                                break
                                            else:
                                                if barcode in rep_mapping:
                                                    if mutation in rep_mapping:
                                                        rep_mapping[barcode][mutation] += mutation_count
                                                    else:
                                                        rep_mapping[barcode][mutation] = mutation_count
                                                else:
                                                    rep_mapping[barcode] = {
                                                            mutation: mutation_count
                                                            }
                                    else:
                                        for n in range(i + 1, len(
                                                lib_dictionary[lib]['mapping'][file]
                                                )):
                                            for ind in lib_dictionary[lib]['mapping'][file][n]:
                                                if barcode in ind_dictionary[file][ind]:
                                                    if mutation in ind_dictionary[file][ind][barcode]:
                                                        gold_mapping[barcode] = {
                                                                mutation: mutation_count
                                                                }
                                                        break
                                            else:
                                                continue
                                            break
                                        else:
                                            if barcode in rep_mapping:
                                                if mutation in rep_mapping:
                                                    rep_mapping[barcode][mutation] += mutation_count
                                                else:
                                                    rep_mapping[barcode][mutation] = mutation_count
                                            else:
                                                rep_mapping[barcode] = {
                                                        mutation: mutation_count
                                                        }
                            elif barcode in control_barcodes['forward']:
                                if barcode == 'bc_1' or barcode == 'bc_2':
                                    mutation = 'wt'
                                else:
                                    mutation = 'del_c'
                                if barcode in control_barcodes_dict[file][rep_name]:
                                    control_barcodes_dict[file][rep_name][barcode] += len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                                else:
                                    control_barcodes_dict[file][rep_name][barcode] = len(
                                            ind_dictionary[file][index][barcode][mutation]
                                            ) * 1000000 / ind_dictionary[file][index]['eff read count']
                    for barcode in rep_mapping:
                        for mutation in rep_mapping[barcode]:
                            mutation_count = rep_mapping[barcode][mutation]
                            if barcode in single_mapping:
                                if mutation in single_mapping:
                                    if barcode in gold_mapping:
                                        gold_mapping[barcode][mutation] = mutation_count
                                    else:
                                        gold_mapping[barcode] = {
                                                mutation: mutation_count
                                                }
                                else:
                                    single_mapping[barcode][mutation] = mutation_count
                            else:
                                single_mapping[barcode] = {
                                        mutation: mutation_count
                                        }
    if not gold_mapping:
        raise Exception('Error! There was not found any valid barcode-mutation pair for mapping of library %s'
                        % lib)

    return gold_mapping, mapping_count, control_barcodes_dict


#Form dictionary of barcodes, found in all replicates and
#dictionary of read count per million effective reads for each replicate.
#gold_barcode_list = {
#                     'barcode1': (barcode1_count1 + barcode1_count2),
#                     'barcode2': (barcode2_count1 + barcode2_count2)
#                    }
#barcode_count_dict = {'file name': {'A10': {
#                                            'barcode1': barcode1_count1,
#                                            'barcode2': barcode2_count1},
#                                    'A11': {
#                                            'barcode1': barcode1_count2,
#                                            'barcode2': barcode2_count2
#                                           }}}
def form_gold_norm_or_expr(ind_dictionary, lib_dictionary, lib, rep_type):
    gold_barcode_dict = {}
    #Includes all barcodes of first replicate.
    first_barcode_dict = {}
    first_file_dict = {}
    #This dectionary includes all barcodes and their counts for each replicate.
    barcode_count_dict = {}
    control_barcodes_dict = {}
    if len(lib_dictionary[lib][rep_type]) == 1:
        for file in lib_dictionary[lib][rep_type]:
            barcode_count_dict[file] = {}
            control_barcodes_dict[file] = {}
            if len(lib_dictionary[lib][rep_type][file]) == 1:
                for index_list in lib_dictionary[lib][rep_type][file]:
                    control_barcodes_dict[file][', '.join(index_list)] = {}
                    rep_eff_read_count = 0
                    for index in index_list:
                        rep_eff_read_count += ind_dictionary[file][index]['eff read count']
                        for barcode in ind_dictionary[file][index]:
                            #Check that it is not 'read count' or control barcode.
                            if len(barcode) > 15 and barcode != 'low quality count':
                                barcode_count = len(
                                        ind_dictionary[file][index][barcode]
                                        )
                                if barcode in gold_barcode_dict:
                                    first_barcode_dict[barcode] += barcode_count
                                else:
                                    first_barcode_dict[barcode] = barcode_count
                            elif barcode in control_barcodes['forward']:
                                if barcode in control_barcodes_dict[file][', '.join(index_list)]:
                                    control_barcodes_dict[file][', '.join(index_list)][barcode] += len(
                                            ind_dictionary[file][index][barcode]
                                            )
                                else:
                                    control_barcodes_dict[file][', '.join(index_list)][barcode] = len(
                                            ind_dictionary[file][index][barcode]
                                            )
                    for control_barcode in control_barcodes_dict[file][', '.join(index_list)]:
                        control_barcodes_dict[file][', '.join(index_list)][control_barcode] = control_barcodes_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
                for barcode in first_barcode_dict:
                    #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
                    if first_barcode_dict[barcode] >= minimum_read_count[rep_type]:
                        gold_barcode_dict[barcode] = first_barcode_dict[barcode] / rep_eff_read_count * 1000000
                barcode_count_dict[file][', '.join(index_list)] = gold_barcode_dict
            else:
                first_rep_eff_read_count = 0
                control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])] = {}
                for index in lib_dictionary[lib][rep_type][file][0]:
                    first_rep_eff_read_count += ind_dictionary[file][index]['eff read count']
                    for barcode in ind_dictionary[file][index]:
                        #Check that it is not 'read count' or control barcode.
                        if len(barcode) > 15 and barcode != 'low quality count':
                            barcode_count = len(
                                    ind_dictionary[file][index][barcode]
                                    )
                            if barcode in first_barcode_dict:
                                first_barcode_dict[barcode] += barcode_count
                            else:
                                first_barcode_dict[barcode] = barcode_count
                        elif barcode in control_barcodes['forward']:
                            if barcode in control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])]:
                                control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][barcode] += len(
                                        ind_dictionary[file][index][barcode]
                                        )
                            else:
                                control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][barcode] = len(
                                        ind_dictionary[file][index][barcode]
                                        )
                for control_barcode in control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])]:
                    control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][control_barcode] = control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])][control_barcode] / first_rep_eff_read_count * 1000000
                for barcode in first_barcode_dict:
                    first_barcode_dict[barcode] = first_barcode_dict[barcode] / first_rep_eff_read_count * 1000000
                barcode_count_dict[file][', '.join(lib_dictionary[lib][rep_type][file][0])] = first_barcode_dict.copy()
                
                for n in range(1, len(lib_dictionary[lib][rep_type][file])):
                    barcode_count_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])] = {}
                    control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])] = {}
                    rep_eff_read_count = 0
                    for ind in lib_dictionary[lib][rep_type][file][n]:
                        rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
                        for control_barcode in control_barcodes['forward']:
                            if control_barcode in ind_dictionary[file][ind]:
                                if control_barcode in control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])]:
                                    control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] += len(
                                            ind_dictionary[file][ind][control_barcode]
                                            )
                                else:
                                    control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] = len(
                                            ind_dictionary[file][ind][control_barcode]
                                            )
                    for control_barcode in control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])]:
                        control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] = control_barcodes_dict[file][', '.join(lib_dictionary[lib][rep_type][file][n])][control_barcode] / rep_eff_read_count * 1000000
                
                for barcode in first_barcode_dict:
                    #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
                    if round(first_barcode_dict[barcode] / 1000000 * first_rep_eff_read_count) >= minimum_read_count[rep_type]:
                        for n in range(1, len(lib_dictionary[lib][rep_type][file])):
                            rep_eff_read_count = 0
                            rep_barcode_count = 0
                            for ind in lib_dictionary[lib][rep_type][file][n]:
                                rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
                                if barcode in ind_dictionary[file][ind]:
                                    rep_barcode_count += len(
                                            ind_dictionary[file][ind][barcode]
                                            )
                            #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
                            if rep_barcode_count < minimum_read_count[rep_type]:
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
            control_barcodes_dict[file] = {}
        first_rep_eff_read_count = 0
        for frst_file in lib_dictionary[lib][rep_type]:
            first_file = frst_file
            break
        
        control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])] = {}
        for index in lib_dictionary[lib][rep_type][first_file][0]:
            first_rep_eff_read_count += ind_dictionary[first_file][index]['eff read count']
            for barcode in ind_dictionary[first_file][index]:
            #Check that it is not 'read count' or control barcode.
                if len(barcode) > 15 and barcode != 'low quality count':
                    barcode_count = len(
                            ind_dictionary[first_file][index][barcode]
                            )
                    if barcode in first_barcode_dict:
                        first_barcode_dict[barcode] += barcode_count
                    else:
                        first_barcode_dict[barcode] = barcode_count
                elif barcode in control_barcodes['forward']:
                    if barcode in control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])]:
                        control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][barcode] += len(
                                ind_dictionary[first_file][index][barcode]
                                )
                    else:
                        control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][barcode] = len(
                                ind_dictionary[first_file][index][barcode]
                                )
        for control_barcode in control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])]:
            control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][control_barcode] = control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])][control_barcode] / first_rep_eff_read_count * 1000000
        
        for barcode in first_barcode_dict:
            first_barcode_dict[barcode] = first_barcode_dict[barcode] / first_rep_eff_read_count * 1000000
        barcode_count_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][0])] = first_barcode_dict.copy()
            
        for n in range(1, len(lib_dictionary[lib][rep_type][first_file])):
            barcode_count_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = {}
            control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = {}
            rep_eff_read_count = 0
            for ind in lib_dictionary[lib][rep_type][first_file][n]:
                rep_eff_read_count += ind_dictionary[first_file][ind]['eff read count']
                for control_barcode in control_barcodes['forward']:
                    if control_barcode in ind_dictionary[first_file][ind]:
                        if control_barcode in control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])]:
                            control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] += len(
                                    ind_dictionary[first_file][ind][control_barcode]
                                    )
            for control_barcode in control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])]:
                control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] = control_barcodes_dict[first_file][', '.join(lib_dictionary[lib][rep_type][first_file][n])] / rep_eff_read_count * 1000000
        
        for barcode in first_barcode_dict:
            #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
            if round(first_barcode_dict[barcode] * first_rep_eff_read_count / 1000000) >= minimum_read_count[rep_type]:
                for n in range(1, len(lib_dictionary[lib][rep_type][first_file])):
                    rep_eff_read_count = 0
                    rep_barcode_count = 0
                    for ind in lib_dictionary[lib][rep_type][first_file][n]:
                        rep_eff_read_count += ind_dictionary[first_file][ind]['eff read count']
                        if barcode in ind_dictionary[first_file][ind]:
                            rep_barcode_count += len(
                                    ind_dictionary[first_file][ind][barcode]
                                    )
                    #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
                    if rep_barcode_count < minimum_read_count[rep_type]:
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
                    control_barcodes_dict[file][', '.join(index_list)] = {}
                    rep_eff_read_count = 0
                    for ind in index_list:
                        rep_eff_read_count += ind_dictionary[file][ind]['eff read count']
                        for control_barcode in control_barcodes['forward']:
                            if control_barcode in ind_dictionary[file][ind]:
                                if control_barcode in control_barcodes_dict[file][', '.join(index_list)]:
                                    control_barcodes_dict[file][', '.join(index_list)][control_barcode] += len(
                                            ind_dictionary[file][ind][control_barcode]
                                            )
                                else:
                                    control_barcodes_dict[file][', '.join(index_list)][control_barcode] = len(
                                            ind_dictionary[file][ind][control_barcode]
                                            )
                    for control_barcode in control_barcodes_dict[file][', '.join(index_list)]:
                        control_barcodes_dict[file][', '.join(index_list)][control_barcode] = control_barcodes_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
                        
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
                        #Check that number of reads and number of reads per million effective reads are both greater then minimum_read_count.
                        if rep_barcode_count < minimum_read_count[rep_type]:
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
        raise Exception('Error! There was not found any barcode for normalization of library %s'
                        % lib)
    
    return gold_barcode_dict, barcode_count_dict, control_barcodes_dict


#This function works only if minimum read count for expression = 0.
def form_rep_expr(ind_dictionary, lib_dictionary, lib):
    barcode_count_dict = {}
    control_barcodes_dict = {}
    total_expression = {}
    total_eff_read_count = 0
    for file in lib_dictionary[lib]['expression']:
        barcode_count_dict[file] = {}
        control_barcodes_dict[file] = {}
        for index_list in lib_dictionary[lib]['expression'][file]:
            barcode_count_dict[file][', '.join(index_list)] = {}
            control_barcodes_dict[file][', '.join(index_list)] = {}
            rep_eff_read_count = 0
            for index in index_list:
                rep_eff_read_count += ind_dictionary[file][index]['eff read count']
                total_eff_read_count += ind_dictionary[file][index]['eff read count']
                for barcode in ind_dictionary[file][index]:
                    #Check that it is not 'read count' or control barcode.
                    if len(barcode) > 15 and barcode != 'low quality count':
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
                    elif barcode in control_barcodes['forward']:
                        if barcode in control_barcodes_dict[file][', '.join(index_list)]:
                            control_barcodes_dict[file][', '.join(index_list)][barcode] += len(
                                    ind_dictionary[file][index][barcode]
                                    )
                        else:
                            control_barcodes_dict[file][', '.join(index_list)][barcode] = len(
                                    ind_dictionary[file][index][barcode]
                                    )
                for barcode in barcode_count_dict[file][', '.join(index_list)]:
                    barcode_count_dict[file][', '.join(index_list)][barcode] = barcode_count_dict[file][', '.join(index_list)][barcode] / rep_eff_read_count * 1000000
                for control_barcode in control_barcodes_dict[file][', '.join(index_list)]:
                    control_barcodes_dict[file][', '.join(index_list)][control_barcode] = control_barcodes_dict[file][', '.join(index_list)][control_barcode] / rep_eff_read_count * 1000000
    for barcode in total_expression:
        total_expression[barcode] = total_expression[barcode] / total_eff_read_count * 1000000
    return barcode_count_dict, control_barcodes_dict, total_expression


#Return True if query is a mutant variant of ref.
def align_mut_variants(ref, query, mode):
    if mode == 'barcode':
        mismatche = barcode_mismatche
        length = barcode_length
    else:
        mismatche = mutation_mismatche
        length = mutation_length
    align_score = pairwise2.align.globalms(ref, query, 1, -0.1, -0.1, -0.1, score_only=True)
    if align_score >= length - mismatche:
        return True
    return False


def not_empty(barc_dict, rep_types, library):
    if not barc_dict:
        raise Exception('Error! There was not found identical barcodes in %s of library %s.'
                        %(rep_types, library))


#Form gold dictionary of identical barcodes
#(only works if minimum_read_count > 0).
#gold_dictionary = {
#                 '29-36': {
#                           'barcodes': {'barcode1': barcode1_count1, 'barcode2': barcode1_count2},
#                           'barcodes-mutations': {
#                                                   'barcode1': {'mut1': mut1_count, 'mut2': mut2_count},
#                                                   'barcode2': {'mut3': mut3_count, 'mut4': mut4_count}
#                                                   }
#                          }
#                 '25-32': {
#                           'barcodes': {'barcode3': barcode1_count3, 'barcode4': barcode1_count4},
#                           'barcodes-mutations': {}
#                           }
#                   }
#rep_dictionary = {'25-32':{
#                  'normalization': {'file name': {
#                                    'A10': {
#                                            'barcode1': barcode1_count1,
#                                            'barcode2': barcode2_count1},
#                                    'A11': {
#                                            'barcode1': barcode1_count2,
#                                            'barcode2': barcode2_count2
#                                           }
#                                                  }},
#                   'expression': {},
#                   'mapping': {'file name': {
#                                 'A1': {
#                                      'barcode1': {
#                                                   'mut1': mut1_count1,
#                                                   'mut2': mut2_count1
#                                                  },
#                                      'barcode2': {
#                                                   'mut3': mut3_count1,
#                                                   'mut4': mut4_count1
#                                                   }
#                                     },
#                               'A2': {
#                                      'barcode1': {
#                                                   'mut1': mut1_count2,
#                                                   'mut2': mut2_count2
#                                                   },
#                                      'barcode2': {
#                                                   'mut3': mut3_count2,
#                                                   'mut4': mut4_count2
#                                                  }
#                                     }
#                                          }
#                              }}}
def form_gold_dictionary(ind_dictionary, lib_dictionary):
    gold_dictionary = {}
    rep_dictionary = {}
    cont_dictionary = {}
    total_rep_dictionary = {}
    for lib in lib_dictionary:
        if lib != 'H':
            gold_dictionary[lib] = {'barcodes': {}, 'barcodes-mutations': {}}
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
                        str(len(gold_mapping)) + ' identical barcodes for mapping of library ' + lib + '.\n'
                        )
                print(lines_of_text[len(lines_of_text) - 1])
                    
            if lib_dictionary[lib]['normalization']:
                gold_normalization, rep_dictionary[lib]['normalization'], cont_dictionary[lib]['normalization'] = form_gold_norm_or_expr(
                        ind_dictionary, lib_dictionary, lib,
                        'normalization'
                        )
                lines_of_text.append(
                        str(len(gold_normalization)) + ' identical barcodes for normalization of library ' + lib + '.\n'
                        )
                print(lines_of_text[len(lines_of_text) - 1])
                    
            if lib_dictionary[lib]['expression']:
                if minimum_read_count['expression'] > 0:
                    gold_expression, rep_dictionary[lib]['expression'], cont_dictionary[lib]['expression'] = form_gold_norm_or_expr(
                            ind_dictionary, lib_dictionary, lib, 'expression'
                            )
                    lines_of_text.append(
                            str(len(gold_expression)) + ' identical barcodes for expression of library ' + lib + '.\n'
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
                        for mutation in gold_mapping[barcode]:
                            barcode_count += gold_mapping[barcode][mutation]
                        gold_dictionary[lib]['barcodes'][barcode] = barcode_count
                        gold_dictionary[lib]['barcodes-mutations'][barcode] = gold_mapping[barcode]
                not_empty(
                        gold_dictionary[lib]['barcodes'],
                        'mapping, normalization and expression',
                        lib
                        )
            elif gold_mapping and gold_normalization:
                for barcode in gold_mapping:
                    if barcode in gold_normalization:
                        barcode_count = gold_normalization[barcode]
                        for mutation in gold_mapping[barcode]:
                            barcode_count += gold_mapping[barcode][mutation]
                        if barcode in total_expr:
                            barcode_count += total_expr[barcode]
                        gold_dictionary[lib]['barcodes'][barcode] = barcode_count
                        gold_dictionary[lib]['barcodes-mutations'][barcode] = gold_mapping[barcode]
                not_empty(
                        gold_dictionary[lib]['barcodes'],
                        'mapping and normalization',
                        lib
                        )
            elif gold_mapping and gold_expression:
                for barcode in gold_mapping:
                    if barcode in gold_expression:
                        barcode_count = gold_expression[barcode]
                        for mutation in gold_mapping[barcode]:
                            barcode_count += gold_mapping[barcode][mutation]
                        gold_dictionary[lib]['barcodes'][barcode] = barcode_count
                        gold_dictionary[lib]['barcodes-mutations'][barcode] = gold_mapping[barcode]
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
                    for mutation in gold_mapping[barcode]:
                        barcode_count += gold_mapping[barcode][mutation]
                    if barcode in total_expr:
                        barcode_count += total_expr[barcode]
                    gold_dictionary[lib]['barcodes'][barcode] = barcode_count
                    gold_dictionary[lib]['barcodes-mutations'][barcode] = gold_mapping[barcode]
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
                        '%s identical barcodes for library %s.\n' %(str(len(gold_dictionary[lib]['barcodes'])), lib)
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
                    df_dict = {'Barcode': [], 'Mutation': [], 'Count': []}
                    for barcode in ind_dictionary[file][index]:
                        if barcode not in ['read count', 'low quality count', 'eff read count']:
                            for mutation in ind_dictionary[file][index][barcode]:
                                df_dict['Barcode'].append(barcode)
                                df_dict['Mutation'].append(mutation)
                                df_dict['Count'].append(len(ind_dictionary[file][index][barcode][mutation]))
                    df = pd.DataFrame(data=df_dict)
                    table_path = output_folder + '\\' + 'Hybrid_test' + index + str(file_count) + '.xlsx'
                    table = df.to_excel(table_path, index = None)
            
    with open(output_folder + '\\gold_dictionary.json', 'w') as f:
        json.dump(gold_dictionary, f)
        
    with open(output_folder + '\\rep_dictionary.json', 'w') as fi:
        json.dump(rep_dictionary, fi)
        
    with open(output_folder + '\\cont_dictionary.json', 'w') as fil:
        json.dump(cont_dictionary, fil)
        
    with open(output_folder + '\\total_rep_dictionary.json', 'w') as file:
        json.dump(total_rep_dictionary, file)
            
    return gold_dictionary, rep_dictionary, cont_dictionary, total_rep_dictionary


gold_dict, rep_dict, cont_dict, total_rep_dict = form_gold_dictionary(fastq_file_dictionary, lib_diction)


#Form gold dictionary of unique barcodes
#uniq_gold_dictionary = {
#                 '29-36': {
#                           'barcodes': {'barcode1': ['barcode2']},
#                           'barcodes-mutations': {
#                                            'barcode1':{
#                                                   'barcode1': ['mut1', 'mut2'],
#                                                   'barcode2': ['mut3', 'mut4'],
#                                                   'major mutation' : {
#                                                              'mut1': ['mut3']
#                                                                      }
#                                                   'total reads': n,
#                                                   'hybrid reads': k
#                                                   }
#                                                 }
#                           }
#                 '25-32': {
#                           'barcodes': {'barcode3': [], 'barcode4': []},
#                           'barcodes-mutations': {}
#                           }
#                         }
def form_uniq_gold_dictionary(gold_dictionary):
    print(datetime.now())
    uniq_gold_dictionary = {}
    hash_length = int(barcode_length // (barcode_mismatche / 1.1 + 1))
    for lib in gold_dictionary:
        if gold_dictionary[lib]['barcodes']:
            uniq_gold_dictionary[lib] = {
                    'barcodes': {}, 'barcodes-mutations': {}
                    }
            hash_dict = {}
            uniq_barcodes_list = []
            mut_variant_dict = {}
            
#Form hash dictionary:
#hash_dict = {
#             'ATGCAG': [
#                        ('TGACTGGATCGAATGCAG', 5),
#                        ('ATGCAGAAGGCTCTTGAT', 3),
#                        ('ATGCAGAAGGCTCTTGAC', 1),
#                       ],
#              'TGCAGA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'GCAGAA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'CAGAAG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'AGAAGG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'GAAGGC': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'AAGGCT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'AGGCTC': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'GGCTCT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'GCTCTT': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'CTCTTG': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'TCTTGA': [('ATGCAGAAGGCTCTTGAT', 3), ('ATGCAGAAGGCTCTTGAC', 1)],
#              'CTTGAT': [('ATGCAGAAGGCTCTTGAT', 3)],
#              'CTTGAC': [('ATGCAGAAGGCTCTTGAC', 1)],
#              'TGACTG': [('TGACTGGATCGAATGCAG', 5)],
#              'GACTGG': [('TGACTGGATCGAATGCAG', 5)],
#              'ACTGGA': [('TGACTGGATCGAATGCAG', 5)],
#              'CTGGAT': [('TGACTGGATCGAATGCAG', 5)],
#              'TGGATC': [('TGACTGGATCGAATGCAG', 5)],
#              'GGATCG': [('TGACTGGATCGAATGCAG', 5)],
#              'GATCGA': [('TGACTGGATCGAATGCAG', 5)],
#              'ATCGAA': [('TGACTGGATCGAATGCAG', 5)],
#              'TCGAAT': [('TGACTGGATCGAATGCAG', 5)],
#              'CGAATG': [('TGACTGGATCGAATGCAG', 5)],
#              'GAATGC': [('TGACTGGATCGAATGCAG', 5)],
#              'AATGCA': [('TGACTGGATCGAATGCAG', 5)],
#              'ATGCAG': [('TGACTGGATCGAATGCAG', 5)]
#             }
            for barcode in gold_dictionary[lib]['barcodes']:
                for i in range(len(barcode) - hash_length):
                    new_hash = barcode[i:i + hash_length]
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
#                      ['TGACTGGATCGAATGCAG', 5],
#                      ['ATGCAGAAGGCTCTTGAT', 3, 'ATGCAGAAGGCTCTTGAC'],
#                      ['ATGCAGAAGGCTCTTGAC', 1]
#                     ]
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
#                                  'barcodes': {
#                                             'TGACTGGATCGAATGCAG': [],
#                                             'ATGCAGAAGGCTCTTGAT':
#                                                        ['ATGCAGAAGGCTCTTGAC']
#                                              },
#                                  'barcodes-mutations': {}
#                                  }
#                       }
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
#If so, check mutations.
#If such unique barcodes have similar mutations, they will be combined.
            for mut_variant in mut_variant_dict:
                if len(mut_variant_dict[mut_variant]) > 1:
                    if gold_dictionary[lib]['barcodes-mutations']:
                        same_barcodes = []
                        for barcode in mut_variant_dict[mut_variant]:
                            for mutation in gold_dictionary[lib]['barcodes-mutations'][mut_variant]:
                                for mut in gold_dictionary[lib]['barcodes-mutations'][barcode]:
                                    if align_mut_variants(
                                            mutation, mut, 'mutation'
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

#Form dictionary 'barcodes-mutations'.
#Find major mutation for each unique barcode.
            barcodes_to_remove = []
            if gold_dictionary[lib]['barcodes-mutations']:
                for uniq_barcode in uniq_gold_dictionary[lib]['barcodes']:
                    uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode] = {
                            uniq_barcode: []
                            }
                    mutation_lists = []
                    for mutation in gold_dictionary[lib]['barcodes-mutations'][uniq_barcode]:
                        mutation_lists.append([
                                mutation,
                                gold_dictionary[lib]['barcodes-mutations'][uniq_barcode][mutation]
                                ])
                        uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode][uniq_barcode].append(
                                mutation
                                )
                    for barcode_mut_var in uniq_gold_dictionary[lib]['barcodes'][uniq_barcode]:
                        uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode][barcode_mut_var] = []
                        for mut in gold_dictionary[lib]['barcodes-mutations'][barcode_mut_var]:
                            for mutation_list in mutation_lists:
                                if mut == mutation_list[0]:
                                    mutation_list[1] += gold_dictionary[lib]['barcodes-mutations'][barcode_mut_var][mut]
                                    break
                            else:
                                mutation_lists.append([
                                        mut,
                                        gold_dictionary[lib]['barcodes-mutations'][barcode_mut_var][mut]
                                        ])
                    mutation_lists.sort(key=lambda x: (x[1], x[0]), reverse=True)
                    total_mutation_count = mutation_lists[len(mutation_lists) - 1][1]
                    for i in range(len(mutation_lists) - 1):
                        total_mutation_count += mutation_lists[i][1]
                        if mutation_lists[i][0] != 'mut var':
                            for n in range(i + 1, len(mutation_lists)):
                                if align_mut_variants(
                                        mutation_lists[i][0],
                                        mutation_lists[n][0],
                                        'mutation'):
                                    mutation_lists[i][1] += mutation_lists[n][1]
                                    mutation_lists[i].append(mutation_lists[n][0])
                                    mutation_lists[n][0] = 'mut var'
                    mutation_lists.sort(key=lambda x: (x[1], x[0]), reverse=True)
                    major_mutation = mutation_lists[0][0]
                    if major_mutation in gold_dictionary[lib]['barcodes-mutations'][uniq_barcode] and mutation_lists[0][1] >= total_mutation_count * major_mutation_portion:
                        uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode]['major mutation'] = {
                                major_mutation: []
                                }
                        if len(mutation_lists[0]) > 2:
                            for k in range(2, len(mutation_lists[0])):
                                uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode]['major mutation'][major_mutation].append(
                                        mutation_lists[0][k]
                                        )
                        uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode]['total reads'] = total_mutation_count
                        uniq_gold_dictionary[lib]['barcodes-mutations'][uniq_barcode]['hybrid reads'] = round(total_mutation_count - mutation_lists[0][1], 6)
                    else:
                        uniq_gold_dictionary[lib]['barcodes-mutations'].pop(
                                uniq_barcode
                                )
                        barcodes_to_remove.append(uniq_barcode)
                for barcode_to_remove in barcodes_to_remove:
                    uniq_gold_dictionary[lib]['barcodes'].pop(
                            barcode_to_remove
                            )
            lines_of_text = [
                    str(len(uniq_gold_dictionary[lib]['barcodes'])) + ' unique barcodes for library ' + lib + '.\n',
                    str(len(barcodes_to_remove)) + ' barcodes were removed because of too high hybrids percent.\n'
                    ]
            print(*lines_of_text)
            with open(run_info_text, 'a') as info_file:
                info_file.writelines(lines_of_text)
            print(datetime.now())
            
    if not os.path.isfile(output_folder + '\\' + 'uniq_gold_dictionary.json'):
        with open(output_folder + '\\uniq_gold_dictionary.json', 'w') as f:
            json.dump(uniq_gold_dictionary, f)
        
    return uniq_gold_dictionary
                            

gold_dict = form_uniq_gold_dictionary(gold_dict)


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
                'Barcode': ['bc_1', 'bc_2', 'bc_3', 'bc_4'],
                'Mutation': ['wt', 'wt', 'del_c', 'del_c']
                }
        
        if uniq_gold_dictionary[lib]['barcodes']:
            if uniq_gold_dictionary[lib]['barcodes-mutations']:
                df_lib_dict = {'Barcode': [], 'Mutation': []}
                for barcode in uniq_gold_dictionary[lib]['barcodes-mutations']:
                    df_lib_dict['Barcode'].append(barcode)
                    for mutation in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['major mutation']:
                        df_lib_dict['Mutation'].append(mutation)
            else:
                df_lib_dict = {'Barcode': []}
                for barcode in uniq_gold_dictionary[lib]['barcodes']:
                    df_lib_dict['Barcode'].append(barcode)
            unique_barcode_count = len(uniq_gold_dictionary[lib]['barcodes'])
            
            for rep_type in rep_dictionary[lib]:
                file_count = 0
                if rep_dictionary[lib][rep_type]:
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
                            replicate_count += 1
                            replicate_list.append(
                                    rep_type + str(file_count) + replicate
                                    )
                            df_dict = {
                                    'Unique Barcode': [],
                                    'Variants of Barcode': [],
                                    'Count': [],
                                    'Total Count': []
                                    }
                            df_lib_dict[rep_type + str(file_count) + replicate] = [
                                    '' for a in range(unique_barcode_count)
                                    ]
                            
                            df_cont_dict[rep_type + str(file_count) + replicate] = [
                                    '', '', '', ''
                                    ]
                            barcode_number = 0
                            for cont_barcode in df_cont_dict['Barcode']:
                                if cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
                                    df_cont_dict[rep_type + str(file_count) + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
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
                                        
                                df_dict['Unique Barcode'].append(barcode)
                                df_dict['Variants of Barcode'].append(
                                        variants_of_barcode
                                        )
                                df_dict['Count'].append(barcode_count)
                                df_dict['Total Count'].append(total_count)
                                df_lib_dict[rep_type + ' mean'][barcode_index] += total_count
                                df_lib_dict[rep_type + str(file_count) + replicate][barcode_index] = total_count
                                
                            df = pd.DataFrame(df_dict, columns=[
                                    'Unique Barcode', 'Variants of Barcode',
                                    'Count', 'Total Count'
                                    ])
                            df = df.sort_values(
                                    by=['Total Count'], ascending=False
                                    )
                            table_path = output_folder + '\\' + lib + rep_type + replicate + str(file_count) + '.xlsx'
                            table = df.to_excel(table_path, index = None)
                    if replicate_count == 0:
                        continue
                    df_lib_dict[rep_type + ' mean'] = [
                            y / replicate_count for y in df_lib_dict[rep_type + ' mean']
                            ]
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
                                sns.kdeplot, n_levels=100, shade=True,
                                shade_lowest=False
                                )
                        grid = grid.map_diag(
                                sns.kdeplot, linewidth=3,
                                shade=True
                                )
                        plt.savefig(
                                output_folder + '\\pair_plot_'+ rep_type + lib + '.eps',
                                format='eps', dpi=1000
                                )
                        plt.close()
                        
                    if total_rep_dictionary[lib][rep_type]:
                        df_total_rep_dict = {
                                    'Unique Barcode': [],
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
                                
                            df_total_rep_dict['Unique Barcode'].append(barcode)
                            df_total_rep_dict['Variants of Barcode'].append(
                                    var_of_barcode
                                    )
                            df_total_rep_dict['Count'].append(barc_count)
                            df_total_rep_dict['Total Count'].append(
                                    total_barc_count
                                    )
                            
                        df_total_rep = pd.DataFrame(df_total_rep_dict, columns=[
                                    'Unique Barcode', 'Variants of Barcode',
                                    'Count', 'Total Count'
                                    ])
                        df_total_rep = df_total_rep.sort_values(
                                by=['Total Count'], ascending=False
                                )
                        total_rep_path = output_folder + '\\' + lib + rep_type + '.xlsx'
                        df_total_rep.to_excel(total_rep_path, index = None)
                        
                else:
                    replicate_count = 0
                    for file in rep_dictionary[lib][rep_type]:
                        file_count += 1
                        for replicate in rep_dictionary[lib][rep_type][file]:
                            replicate_count += 1
                            df_dict = {
                                    'Unique Barcode': [],
                                    'Major Mutation': [],
                                    'Variants of Barcode': [],
                                    'Mutations' : [],
                                    'Count': [],
                                    'Total Count': [],
                                    'Hybrid Reads': []
                                    }
                            df_cont_dict[rep_type + str(file_count) + replicate] = [
                                    '', '', '', ''
                                    ]
                            
                            barcode_number = 0
                            for cont_barcode in df_cont_dict['Barcode']:
                                if cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
                                    df_cont_dict[rep_type + str(file_count) + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
                                    df_cont_dict[rep_type + ' mean'][barcode_number] += cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
                                barcode_number += 1
                            
                            for barcode in uniq_gold_dictionary[lib]['barcodes-mutations']:
                                total_count = 0
                                hybrid_count = 0
                                for mutation in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['major mutation']:
                                    major_mutation = mutation
                                variants_of_barcode = barcode + '\n'
                                mutations = major_mutation + '\n'
                                if barcode in rep_dictionary[lib][rep_type][file][replicate]:
                                    if major_mutation in rep_dictionary[lib][rep_type][file][replicate][barcode]:
                                        count = str(rep_dictionary[lib][rep_type][file][replicate][barcode][major_mutation]) + '\n'
                                        total_count += rep_dictionary[lib][rep_type][file][replicate][barcode][major_mutation]
                                    else:
                                        count = '0\n'
                                else:
                                    count = '0\n'
                                for mutation_variant in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode][barcode]:
                                    if mutation_variant != major_mutation:
                                        variants_of_barcode += '\n'
                                        mutations += mutation_variant + '\n'
                                        if barcode in rep_dictionary[lib][rep_type][file][replicate]:
                                            if mutation_variant in rep_dictionary[lib][rep_type][file][replicate][barcode]:
                                                count += str(rep_dictionary[lib][rep_type][file][replicate][barcode][mutation_variant]) + '\n'
                                                total_count += rep_dictionary[lib][rep_type][file][replicate][barcode][mutation_variant]
                                                if mutation_variant not in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['major mutation'][major_mutation]:
                                                    hybrid_count += rep_dictionary[lib][rep_type][file][replicate][barcode][mutation_variant]
                                            else:
                                                count += '0\n'
                                        else:
                                            count += '0\n'
                                for barc_variant in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]:
                                    if len(barc_variant) > 15 and barc_variant != barcode:
                                        variants_of_barcode += barc_variant
                                        for mut_variant in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode][barc_variant]:
                                            variants_of_barcode += '\n'
                                            mutations += mut_variant + '\n'
                                            if barc_variant in rep_dictionary[lib][rep_type][file][replicate]:
                                                if mut_variant in rep_dictionary[lib][rep_type][file][replicate][barc_variant]:
                                                    count += str(rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant])
                                                    total_count += rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant]
                                                    if mut_variant != major_mutation and mut_variant not in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['major mutation'][major_mutation]:
                                                        hybrid_count += rep_dictionary[lib][rep_type][file][replicate][barc_variant][mut_variant]
                                                else:
                                                    count += '0\n'
                                            else:
                                                count += '0\n'
                                df_dict['Unique Barcode'].append(barcode)
                                df_dict['Major Mutation'].append(
                                        major_mutation
                                        )
                                df_dict['Variants of Barcode'].append(
                                        variants_of_barcode
                                        )
                                df_dict['Mutations'].append(mutations)
                                df_dict['Count'].append(count)
                                df_dict['Total Count'].append(total_count)
                                df_dict['Hybrid Reads'].append(hybrid_count)
                                
                            df = pd.DataFrame(df_dict, columns=[
                                    'Unique Barcode', 'Major Mutation',
                                    'Variants of Barcode', 'Mutations',
                                    'Count', 'Total Count', 'Hybrid Reads'
                                    ])
                            df['Hybrids Percent'] = df['Hybrid Reads'] / df['Total Count'] * 100
                            df = df.sort_values(
                                    by=['Total Count'], ascending=False
                                    )
                            table_path = output_folder + '\\' + lib + '_mapping_' + replicate + str(file_count) + '.xlsx'
                            table = df.to_excel(table_path, index = None)
                    if replicate_count == 0:
                        continue        
                    df_cont_dict[rep_type + ' mean'] = [
                            d / replicate_count for d in df_cont_dict[rep_type + ' mean']
                            ]
                    
                    df_total_map_dict = {
                                    'Unique Barcode': [],
                                    'Major Mutation': [],
                                    'Variants of Barcode': [],
                                    'Mutations' : [],
                                    'Count': [],
                                    'Total Count': [],
                                    'Hybrid Reads': [],
                                    'Hybrids Percent': []
                                    }
                    for barcode in uniq_gold_dictionary[lib]['barcodes-mutations']:
                        for mutation in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['major mutation']:
                            major_mut = mutation
                        total_map_count = uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['total reads'] / replicate_count
                        hybr_count = uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]['hybrid reads'] / replicate_count
                        hybrid_percent = hybr_count / total_map_count * 100
                        var_of_barc = barcode + '\n'
                        mut = major_mut + '\n'
                        map_count = str(total_rep_dictionary[lib]['mapping'][barcode][major_mut] / replicate_count) + '\n'
                        for mut_var in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode][barcode]:
                            if mut_var != major_mut:
                                var_of_barc += '\n'
                                mut += mut_var + '\n'
                                map_count += str(total_rep_dictionary[lib]['mapping'][barcode][mut_var] / replicate_count) + '\n'
                        for barc_var in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode]:
                            if len(barc_var) > 15 and barc_var != barcode:
                                var_of_barcode += barc_var
                                for mut_variant in uniq_gold_dictionary[lib]['barcodes-mutations'][barcode][barc_var]:
                                    var_of_barc += '\n'
                                    mut += mut_variant + '\n'
                                    map_count += str(total_rep_dictionary[lib]['mapping'][barc_var][mut_variant] / replicate_count) + '\n'
                        df_total_map_dict['Unique Barcode'].append(barcode)
                        df_total_map_dict['Major Mutation'].append(major_mut)
                        df_total_map_dict['Variants of Barcode'].append(
                                var_of_barc
                                )
                        df_total_map_dict['Mutations'].append(mut)
                        df_total_map_dict['Count'].append(map_count)
                        df_total_map_dict['Total Count'].append(
                                total_map_count
                                )
                        df_total_map_dict['Hybrid Reads'].append(hybr_count)
                        df_total_map_dict['Hybrids Percent'].append(
                                hybrid_percent
                                )
                        
                    df_total_map = pd.DataFrame(df_total_map_dict, columns=[
                                    'Unique Barcode', 'Major Mutation',
                                    'Variants of Barcode', 'Mutations',
                                    'Count', 'Total Count', 'Hybrid Reads',
                                    'Hybrids Percent'
                                    ])
                    df_total_map = df_total_map.sort_values(
                            by=['Total Count'], ascending=False
                            )
                    total_map_table_path = output_folder + '\\' + lib + '_mapping ' + '.xlsx'
                    total_map_table = df_total_map.to_excel(
                            total_map_table_path, index = None
                            )
            
            columns_list = [key for key in df_lib_dict]
            df_lib = pd.DataFrame(data=df_lib_dict)
            new_columns_list = ['Barcode']
            pos_start = 1
            
            if 'Mutation' in columns_list:
                new_columns_list.append('Mutation')
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
                    'Barcode', 'Mutation'
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
                        df_lib['normalyzed to control'] = df_lib['normalyzed expression'] / df_cont_dict['normalyzed expression mean'][0]
                        
                        df_cont_dict['difference'] = [
                                df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
                                df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
                                df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0],
                                df_cont_dict['normalyzed expression mean'][2] / df_cont_dict['normalyzed expression mean'][0]
                                ]
                        
                df_lib = df_lib.sort_values(
                        by=['normalyzed expression'], ascending=False
                        )
            
            total_table_path = output_folder + '\\' + lib + '.xlsx'
            total_table = df_lib.to_excel(total_table_path, index = None)
            
            if 'Mutation' in new_columns_list:
                #Pairplot for mutations that have more than two barcodes.
                lib_dict = df_lib.to_dict('index')
                mutation_dictionary = {}
                for index in lib_dict:
                    new_mutation = lib_dict[index]['Mutation']
                    if new_mutation in mutation_dictionary:
                        mutation_dictionary[new_mutation].append(lib_dict[index]['normalyzed to control'])
                    else:
                        mutation_dictionary[new_mutation] = [lib_dict[index]['normalyzed to control']]
                first_mut_pairplot_list = []
                second_mut_pairplot_list = []
                
                mut_barcode_count_list = []
                for new_mutation in mutation_dictionary:
                    mut_barcode_count_list.append(
                            len(mutation_dictionary[new_mutation])
                            )
                    if len(mutation_dictionary[new_mutation]) >= 2:
                        first_random_expr = randint(0, len(mutation_dictionary[new_mutation]) - 1)
                        first_mut_pairplot_list.append(mutation_dictionary[new_mutation][first_random_expr])
                        second_random_expr = randint(0, len(mutation_dictionary[new_mutation]) - 1)
                        while second_random_expr == first_random_expr:
                            second_random_expr = randint(0, len(mutation_dictionary[new_mutation]) - 1)
                        second_mut_pairplot_list.append(mutation_dictionary[new_mutation][second_random_expr])
                
                less_11 = 0
                for m in range(1, 11):
                    lines_of_text.append(
                            '%s mutations with %s barcodes for lib %s.\n' %(str(mut_barcode_count_list.count(m)), str(m), lib))
                    less_11 += mut_barcode_count_list.count(m)
                lines_of_text.append(
                        str(len(mut_barcode_count_list) - less_11) + ' mutations with more than 10 barcodes for lib ' + lib + '.\n')
                
                plt.figure(figsize=(10,10))
                sns.set_style('ticks')
                sns.set_context('poster')
                sns.distplot(
                        mut_barcode_count_list, kde=False,
                        bins=max(mut_barcode_count_list),
                        hist_kws=dict(linewidth=0),
                        axlabel='Number of barcodes for each mutation for lib ' + lib
                        )
                sns.despine()
                plt.savefig(
                        output_folder + '\\barc_per_mut_' + lib + '.eps',
                        format='eps', dpi=1000
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
                        sns.kdeplot, n_levels=30, shade=True,
                        shade_lowest=False
                        )
                grid = grid.map_diag(
                        sns.kdeplot, linewidth=3, shade=True
                        )
                
                plt.savefig(
                        output_folder + '\\mutation_with_different_barcodes_' + lib + '.eps',
                        format='eps', dpi=1000
                        )
                plt.close()
                
                #Plot mean normalysed expressions for all mutations.
                if df_cont_dict['normalyzed expression mean'][0] > 0:
                    mut_ser = df_lib.groupby('Mutation')['normalyzed to control'].mean()
                else:
                    mut_ser = df_lib.groupby('Mutation')['normalyzed expression'].mean()
                mut_ser = mut_ser.sort_values(ascending=False)
                mut_table_path = output_folder + '\\' + lib + '_mutations' + '.xlsx'
                mut_table = mut_ser.to_excel(mut_table_path)
                
                plt.figure(figsize=(10,10))
                sns.set_style('ticks')
                sns.set_context('poster')
                mut_ser_list = mut_ser.tolist()
                sns.distplot(
                        mut_ser, bins=round(mut_ser_list[0] * 10),
                        kde_kws={'color': 'g', 'linewidth': 6, 'label': 'total KDE'},
                        hist_kws=dict(linewidth=0, alpha=0.1)
                        )
                if 'normalyzed expression mean' in df_cont_dict:
                    plt.axvline(1, 0,1, linewidth=4, color='r')
                    plt.axvline(
                            df_cont_dict['difference'][0], 0,1, linewidth=4,
                            color='r'
                            )
                sns.despine()
                plt.savefig(
                        output_folder + '\\kdeplot_' + lib + '.eps',
                        format='eps', dpi=1000
                        )
                if wt_mutation[lib] in mutation_dictionary:
                    lines_of_text.append(
                            'In library %s %s wt mutations.\n' %(lib, str(len(mutation_dictionary[wt_mutation[lib]])))
                            )
                    sns.distplot(
                            mutation_dictionary[wt_mutation[lib]],
                            kde_kws={'color': 'k', 'linewidth': 3, 'label': 'wt KDE'},
                            hist_kws=dict(linewidth=0, alpha=0.1)
                            )
                    sns.despine()
                    plt.savefig(
                        output_folder + '\\kdeplot_and_wt_' + lib + '.eps',
                        format='eps', dpi=1000
                        )
                plt.close()
                
                log_2_mut = [math.log2(x) for x in mut_ser_list if x > 0]
                stat, p = stats.shapiro(log_2_mut)
                if p > 0.05:
                    lines_of_text.append(
                            'For mutations in lib %s log2 of normalyzed expressions looks Gaussian.\nShapiro p-value=%s.\n' % (lib, str(p))
                            )
                else:
                    statis, p_value = stats.normaltest(log_2_mut)
                    if p_value > 0.05:
                        lines_of_text.append(
                                'For mutations in lib %s log2 of normalyzed expressions looks Gaussian.\nShapiro p-value=%s, K^2 p-value=%s.\n' % (lib, str(p), str(p_value))
                                )
                    else:
                        lines_of_text.append(
                                'For mutations in lib %s log2 of normalyzed expressions does not look Gaussian.\nShapiro p-value=%s, K^2 p-value=%s.\n' % (lib, str(p), str(p_value))
                                )
                    
                plt.figure(figsize=(10,10))
                sns.set_style('ticks')
                sns.set_context('poster')
                sns.distplot(
                        log_2_mut, bins=round(log_2_mut[0] * 20),
                        kde_kws={'color': 'g', 'linewidth': 6, 'label': 'total KDE'},
                        hist_kws=dict(linewidth=0, alpha=0.1)
                        )
                if 'normalyzed expression mean' in df_cont_dict:
                    plt.axvline(0, 0,1, linewidth=4, color='r')
                    plt.axvline(
                            math.log2(df_cont_dict['difference'][0]), 0,1,
                            linewidth=4,
                            color='r'
                            )
                sns.despine()
                plt.savefig(
                        output_folder + '\\log2_kdeplot_' + lib + '.eps',
                        format='eps', dpi=1000
                        )
                if wt_mutation[lib] in mutation_dictionary:
                    wt_mut_log2[lib] = [math.log2(y) for y in mutation_dictionary[wt_mutation[lib]] if y > 0]
                    lines_of_text.append('In library %s %s wt mutations above zero.\n' %(lib, str(len(wt_mut_log2[lib]))))
                    statistic, p_val = stats.shapiro(wt_mut_log2[lib])
                    if p_val > 0.05:
                        lines_of_text.append('For wt in lib %s log2 of normalyzed expressions looks Gaussian.\np-value=%s.\n' % (lib, str(p_val)))
                        log_2_wt_gaus[lib] = True
                    else:
                        lines_of_text.append('For wt in lib %s log2 of normalyzed expressions does not looks Gaussian.\np-value=%s.\n' % (lib, str(p_val)))
                        log_2_wt_gaus[lib] = False
                    
                    sns.distplot(
                            wt_mut_log2[lib],
                            kde_kws={'color': 'k', 'linewidth': 3, 'label': 'wt KDE'},
                            hist_kws=dict(linewidth=0, alpha=0.1)
                            )
                    sns.despine()
                    plt.savefig(
                        output_folder + '\\log2_kdeplot_and_wt_' + lib + '.eps',
                        format='eps', dpi=1000
                        )
                plt.close()
                
                mut_dict[lib] = list(zip(mut_ser.index, mut_ser_list))
                mut_log2_dict[lib] = list(zip(mut_ser.index, log_2_mut))

            
            df_cont = pd.DataFrame(data=df_cont_dict)
            cont_columns_list = df_cont.columns.tolist()
            cont_new_columns_list += [
                    l for l in cont_columns_list if l.startswith('mapping') and not l.endswith('mean')
                    ]
            
            if 'mapping mean' in cont_columns_list:
                cont_new_columns_list.append('mapping mean')
            df_cont = df_cont[cont_new_columns_list]
            cont_table_path = output_folder + '\\' + lib + '_control' + '.xlsx'
            cont_table = df_cont.to_excel(cont_table_path, index = None)
                
        else:
            for rep_type in rep_dictionary[lib]:
                df_cont_dict[rep_type + ' mean'] = [0, 0, 0, 0]
                file_count = 0
                for file in rep_dictionary[lib][rep_type]:
                    file_count += 1
                    for replicate in rep_dictionary[lib][rep_type][file]:
                        df_cont_dict[rep_type + str(file_count) + replicate] = [
                                    '', '', '', ''
                                    ]
                        barcode_number = 0
                        for cont_barcode in df_cont_dict['Barcode']:
                            if cont_barcode in cont_dictionary[lib][rep_type][file][replicate]:
                                df_cont_dict[rep_type + str(file_count) + replicate][barcode_number] = cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
                                df_cont_dict[rep_type + ' mean'][barcode_number] += cont_dictionary[lib][rep_type][file][replicate][cont_barcode]
                            barcode_number += 1
                        
                        rep_barcodes = {
                                lib: {'barcodes': {}, 'barcodes-mutations': {}}
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
                            df_dict['Unique Barcode'].append(barcode)
                            df_dict['Variants of Barcode'].append(
                                    variants_of_barcode
                                    )
                            df_dict['Count'].append(barcode_count)
                            df_dict['Total Count'].append(total_count)
                        df = pd.DataFrame(data=df_dict)
                        df = df.sort_values(
                                by=['Total Count'], ascending=False
                                )
                        table_path = output_folder + '\\' + lib + rep_type + replicate + str(file_count) + '.xlsx'
                        table = df.to_excel(table_path, index = None)
            
            df_cont = pd.DataFrame(data=df_cont_dict)
            cont_table_path = output_folder + '\\' + lib + '_control' + '.xlsx'
            cont_table = df_cont.to_excel(cont_table_path, index = None)
    
    print(*lines_of_text)
    with open(run_info_text, 'a') as info_file:
        info_file.writelines(lines_of_text)
    
    with open(output_folder + '\\mut_dictionary.json', 'w') as f:
        json.dump(mut_dict, f)
    
    with open(output_folder + '\\mut_log2_dictionary.json', 'w') as file:
        json.dump(mut_log2_dict, file)
    
    with open(output_folder + '\\wt_log2.json', 'w') as fil:
        json.dump(wt_mut_log2, fil)
        
    return mut_dict, mut_log2_dict, wt_mut_log2, log_2_wt_gaus


mut_dictionary, mut_log2_dictionary, wt_log2, log_2_wt_gaussian = make_tables(gold_dict, rep_dict, cont_dict, total_rep_dict)


def make_graphs(uniq_gold_dictionary):
    for lib in uniq_gold_dictionary:
        if uniq_gold_dictionary[lib]['barcodes']:
            barcode_len_list = [
                    len(x) for x in uniq_gold_dictionary[lib]['barcodes']
                    ]
            
            plt.figure(figsize=(10,10))
            sns.set_style('ticks')
            sns.set_context('poster')
            sns.distplot(
                    barcode_len_list, kde=False, bins=4,
                    hist_kws=dict(linewidth=0),
                    axlabel='Barcode length for library ' + lib
                    )
            sns.despine()
            plt.savefig(
                    output_folder + '\\barcode_len_' + lib + '.eps',
                    format='eps', dpi=1000
                    )
            plt.close()


make_graphs(gold_dict)
                

def make_pwm(mut_dict):
    lines_of_text = []
    rounded_mismatche = round(mutation_mismatche)
    for lib in mut_dict:
        nucleotide_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        total_count = 0
        for mut_tuple in mut_dict[lib]:
            for nucleotide in mut_tuple[0]:
                nucleotide_count[nucleotide] += 1
                total_count += 1
        lines_of_text.append('Mutations in lib ' + lib + ' contains\n')
        print('Mutations in lib ' + lib + ' contains ')
        for nucl in nucleotide_count:
            nucleotide_count[nucl] = nucleotide_count[nucl] / total_count
            lines_of_text.append(
                    str(nucleotide_count[nucl]) + ' of ' + nucl + '.\n'
                    )
            print(str(nucleotide_count[nucl]) + ' of ' + nucl)
        for i in range(
                mutation_length - rounded_mismatche, mutation_length + rounded_mismatche + 1
                ):
            instances = [Seq(k[0]) for k in mut_dict[lib] if len(k[0]) == i]
            lines_of_text.append(
                    str(len(instances)) + ' mutations ' + str(i) + ' nucleotides long for lib ' + lib + '.\n'
                    )
            print(str(len(instances)) + ' mutations ' + str(i) + ' nucleotides long for lib ' + lib)
            if instances:
                m = motifs.create(instances)
                pwm = m.counts.normalize(pseudocounts=nucleotide_count)
                print('PWM for all mutations' + str(i) + ' bp long in lib ' + lib)
                print(pwm)
                print(pwm.consensus)
                print(pwm.degenerate_consensus)
                if len(instances) > 20:
                    m_high = motifs.create(instances[:round(len(instances) * 0.05)])
                    pwm_high = m_high.counts.normalize(pseudocounts=nucleotide_count)
                    print('PWM for high expressed mutations' + str(i) + ' bp long in lib ' + lib)
                    print(pwm_high)
                    print(pwm_high.consensus)
                    print(pwm_high.degenerate_consensus)
                    m_low = motifs.create(instances[round(len(instances) * 0.95):])
                    pwm_low = m_low.counts.normalize(pseudocounts=nucleotide_count)
                    print('PWM for low expressed mutations' + str(i) + ' bp long in lib ' + lib)
                    print(pwm_low)
                    print(pwm_low.consensus)
                    print(pwm_low.degenerate_consensus)
    with open(run_info_text, 'a') as info_file:
        info_file.writelines(lines_of_text)
            
            
make_pwm(mut_dictionary)


def find_k_mer(mut_dict, len_k_mer):
    k_mer_dict = {}
    k_mer_count = 0
    for mut_tuple in mut_dict:
        mutation = mut_tuple[0]
        for i in range(len(mutation) - len_k_mer + 1):
            k_mer = mutation[i:i + len_k_mer]
            k_mer_count += 1
            if k_mer not in k_mer_dict:
                k_mer_dict[k_mer] = 1
            else:
                k_mer_dict[k_mer] += 1
    for k_mer in k_mer_dict:
        k_mer_dict[k_mer] = k_mer_dict[k_mer] / k_mer_count * 100
    return k_mer_dict


def count_k_mer(mutation_dict):
    lines_of_text = []
    for k_mer_len in range(4, 8):
        for lib in mutation_dict:
            signif_k_mer = {
                'k_mer': [], 'Total count': [], 'High count': [], 'Low count': []
                }
            k_mer_diction = find_k_mer(mutation_dict[lib], k_mer_len)
            lines_of_text.append(
                    str(len(k_mer_diction)) + ' k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
                    )
            k_mer_high_diction = find_k_mer(
                    mutation_dict[lib][:round(len(mutation_dict[lib]) * 0.05)], k_mer_len
                    )
            lines_of_text.append(
                    str(len(k_mer_high_diction)) + ' high k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
                    )
            k_mer_low_diction = find_k_mer(
                    mutation_dict[lib][round(len(mutation_dict[lib]) * 0.95):], k_mer_len
                    )
            lines_of_text.append(
                    str(len(k_mer_low_diction)) + ' low k-mers ' + str(k_mer_len) +  ' bp long were found in library ' + lib + '.\n'
                    )
            print(*lines_of_text[len(lines_of_text) - 3:])
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
            table_path = output_folder + '\\' + lib + 'k_mers' + str(k_mer_len) + '.xlsx'
            table = df.to_excel(table_path, index = None)


count_k_mer(mut_dictionary)


def find_log2_k_mer(mutation_lst, len_k_mer):
    k_mer_dict = {}
    for mut_tuple in mutation_lst:
        mutation = mut_tuple[0]
        for i in range(len(mutation) - len_k_mer + 1):
            k_mer = mutation[i:i + len_k_mer]
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
            table_path = output_folder + '\\' + lib + 'stats_k_mers' + str(k_mer_len) + '.xlsx'
            table = df.to_excel(table_path, index = None)


log2_k_mer_stat_test(mut_log2_dictionary, wt_log2, log_2_wt_gaussian)
