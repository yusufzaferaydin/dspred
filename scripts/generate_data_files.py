__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "26 June 2018"

usage = """
	Usage: python generate_obs_files.py <psiblast_pssm_filename> <hhblits_pssm_filename> <obs_filename_psiblast_dbn> <obs_filename_hhmake_dbn>
"""

def read_psiblast_pssm_from_file(pssm_filename, normalization_type):

        pssm_file = open( pssm_filename, 'r')

        pssm = []

        #Skip the first three lines of the PSSM file
        line = pssm_file.readline()
        line = pssm_file.readline()
        line = pssm_file.readline()

        done = 0
        number_of_aas = 0
        number_of_columns = 20

        while not done:

                line = pssm_file.readline()

                if ( len(line) < 10 ):

                        done = 1
                        continue

                token_list = line.split()
                aa_index = int(token_list[0])
                aa_type = token_list[1]

                pssm_row = []

                for i in range(0, number_of_columns, 1):

                        pssm_value = int(token_list[i+2])

                        if (normalization_type == 1): #sigmoidal transformation

                                transformed_pssm = 1 / (1 + exp(-1*pssm_value))

                        elif (normalization_type == 2): #linear mapping

                                if (pssm_value > 7):

                                        transformed_pssm = 1.0

                                elif (pssm_value < -7):

                                        transformed_pssm = 0.0
                                else:
                                        transformed_pssm = pssm_value / 14.0 + 0.5

                        pssm_row += [ transformed_pssm ]

                pssm.extend( [pssm_row] )
                number_of_aas += 1

        pssm_file.close()

        return pssm

def give_me_aa_number_hmmer( aa ):

        aa_labels = 'ACDEFGHIKLMNPQRSTVWY'
        aa_number = aa_labels.find(aa)

        if ( aa_number == -1 ):

                return 0
        else:
                return aa_number

def read_hhsearch_pssm_from_file(pssm_filename, sigmoid_transformation, linear_transformation, seq_length):

        pssm_file = open( pssm_filename, 'r')

        pssm = []

        done_skip_header = 0
        while not done_skip_header:

                line = pssm_file.readline()
                if ("NULL" in line):
                        line = line.rstrip()
                        token_list = line.split()
                        background_profile = token_list[1:21]
                        line = pssm_file.readline()
                        line = pssm_file.readline()
                        line = pssm_file.readline()
                        done_skip_header = 1

        aa_labels = 'ARNDCQEGHILKMFPSTWYV'
        done = 0
        while not done:

                line = pssm_file.readline().rstrip()

                if (len(line) == 0):
                        done = 1
                        continue

                token_list = line.split()

                if (token_list[0] == '//'):
                        done = 1
                        continue

                if (token_list[0] == '!'):
                        line = pssm_file.readline().rstrip()
                        line = pssm_file.readline().rstrip()
                        continue

                #if (token_list[0] == 'X'):
                #        line = pssm_file.readline().rstrip()
                #        line = pssm_file.readline().rstrip()
                #       continue

                aa_index = int(token_list[1])
                if (aa_index > seq_length):
                        done = 1
                        continue

                profile_values = token_list[2:22]
                pssm_row = []

                for i in range(20):

                        aa_index_hmmer = give_me_aa_number_hmmer(aa_labels[i])

                        if (str(profile_values[aa_index_hmmer]) == '*'):
                                background_value = float(background_profile[aa_index_hmmer])/1000.0
                                pssm_value = -0.00001
                        else:
                                pssm_value = float(profile_values[aa_index_hmmer])/1000.0
                                background_value = float(background_profile[aa_index_hmmer])/1000.0
                                #pssm_value = log(exp(-1*pssm_value) / exp(-1*background_value)) #we convert to -log(p) to probability first
                                pssm_value = -1*pssm_value + background_value

                        pssm_row += [ pssm_value ]

                pssm.extend( [pssm_row] )

                #Skip two lines (for other states in HMM-profile)
                line = pssm_file.readline().rstrip()
                line = pssm_file.readline().rstrip()

        #Find the minimum and maximum values
        n_rows = len(pssm)
        n_cols = len(pssm[0])

        min_value = 1e9
        max_value = -1e9
        for i in range(n_rows):
                for j in range(n_cols):
                        pssm_value = pssm[i][j]
                        if (pssm_value == -0.00001):
                                continue
                        if (pssm_value > max_value):
                                max_value = pssm_value
                        if (pssm_value < min_value):
                                min_value = pssm_value

        if (sigmoid_transformation == 1):
                for i in range(n_rows):
                        for j in range(20):
                                pssm_value = pssm[i][j]
                                if (pssm_value == -0.00001):
                                        transformed_pssm = 0.0
                                else:
                                        if (-1*pssm_value > 500):
                                                transformed_pssm = 0.0
                                        else:
                                                transformed_pssm = 1 / (1 + exp(-1*pssm_value))

                                pssm[i][j] = transformed_pssm

        elif (linear_transformation == 1):
                for i in range(n_rows):
                        for j in range(20):
                                pssm_value = pssm[i][j]
                                if (pssm_value == -0.00001):
                                        transformed_pssm = -10
                                else:
                                        transformed_pssm = (pssm_value - min_value) / (max_value - min_value)
                                pssm[i][j] = transformed_pssm
        pssm_file.close()
        return pssm

def read_hhblits_pssm_from_file(pssm_filename, normalization_type, seq_length):

        if (normalization_type == 1):
                pssm = read_hhsearch_pssm_from_file(pssm_filename, 1, 0, seq_length)

        elif (normalization_type == 2):
                pssm = read_hhsearch_pssm_from_file(pssm_filename, 0, 1, seq_length)

        return pssm


def generate_obs_file_test(pssm, seq_length, pssm_window_size, n_dim, obs_filename):

	obs_file = open(obs_filename, 'w')
	
	for i in range(pssm_window_size):
                for j in range(n_dim):
                        obs_file.write('%.7f ' % 0.0)
                obs_file.write('\n')

        for i in range(seq_length):
		for j in range(n_dim):
			obs_file.write('%.7f ' % pssm[i][j])	
		obs_file.write('\n')

	obs_file.close()

def reverse_pssm( pssm ):

        n_rows = len( pssm )
        n_iterations = n_rows / 2

        for i in range(n_iterations):

                temp = pssm[i]
                pssm[i] = pssm[n_rows-1-i]
                pssm[n_rows-1-i] = temp

        return pssm

import sys
from math import *

if __name__ == "__main__":

	if (len(sys.argv) != 8):
		sys.stderr.write(usage)
		sys.exit(1)

	psiblast_pssm_filename = sys.argv[1]
	hhblits_pssm_filename = sys.argv[2]
	pssm_window_size = int(sys.argv[3])
	obs_filename_psiblast_dbn_nc = sys.argv[4]
        obs_filename_psiblast_dbn_cn = sys.argv[5]
	obs_filename_hhmake_dbn_nc = sys.argv[6]
        obs_filename_hhmake_dbn_cn = sys.argv[7]

	#check generate_obs_file function in generate_data.py	

	n_dim = 20
	
	normalization_type = 1 #sigmoidal transformation
        psiblast_pssm = read_psiblast_pssm_from_file(psiblast_pssm_filename, normalization_type)
	seq_length = len(psiblast_pssm)

        hhmake_pssm = read_hhblits_pssm_from_file(hhblits_pssm_filename, normalization_type, seq_length)

	generate_obs_file_test(psiblast_pssm, seq_length, pssm_window_size, n_dim, obs_filename_psiblast_dbn_nc)
	psiblast_pssm_reverse = reverse_pssm(psiblast_pssm)	
        generate_obs_file_test(psiblast_pssm_reverse, seq_length, pssm_window_size, n_dim, obs_filename_psiblast_dbn_cn)

        generate_obs_file_test(hhmake_pssm, seq_length, pssm_window_size, n_dim, obs_filename_hhmake_dbn_nc)
        hhmake_pssm_reverse = reverse_pssm(hhmake_pssm)
        generate_obs_file_test(hhmake_pssm_reverse, seq_length, pssm_window_size, n_dim, obs_filename_hhmake_dbn_cn)

        #if ((prediction_task == 1) or (prediction_task == 2)):
        #	label_to_append = 2
	#else:
        #        label_to_append = 0

	

