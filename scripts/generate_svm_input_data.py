__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "1 July 2018"

usage = """
	Usage: python generate_svm_input_data.py <dat_filename_psiblast> <dat_filename_hhmake> <average_postp_filename_psiblast> <average_postp_filename_hhmake> <average_distribution_3_filename> <pssm_window_size_two_sided> <svm_input_data_filename> 
"""

import sys

if __name__ == "__main__":

	if (len(sys.argv) != 8):
		sys.stderr.write(usage)
		sys.exit(1)

	dat_filename_psiblast = sys.argv[1]
	dat_filename_hhmake = sys.argv[2]
	average_postp_filename_psiblast = sys.argv[3]
	average_postp_filename_hhmake = sys.argv[4]
	average_distribution_3_filename = sys.argv[5]
	pssm_window_size_two_sided = int(sys.argv[6])
	svm_input_data_filename = sys.argv[7]

	dat_file_psiblast = open(dat_filename_psiblast, 'r')
	dat_file_psiblast_inside = dat_file_psiblast.readlines()
	dat_file_psiblast.close()
	n_lines_psiblast = len(dat_file_psiblast_inside)
	token_list = dat_file_psiblast_inside[0].rstrip().split()
	n_dim_psiblast = len(token_list)

	dat_file_hhmake = open(dat_filename_hhmake, 'r')
	dat_file_hhmake_inside = dat_file_hhmake.readlines()
	dat_file_hhmake.close()
	n_lines_hhmake = len(dat_file_hhmake_inside)
	token_list = dat_file_hhmake_inside[0].rstrip().split()
	n_dim_hhmake = len(token_list)

	average_distribution_3_file = open(average_distribution_3_filename, 'r')
	average_distribution_3_file_inside = average_distribution_3_file.readlines()
	average_distribution_3_file.close()
	n_lines_distribution_3 = len(average_distribution_3_file_inside)
	token_list = average_distribution_3_file_inside[0].rstrip().split()
	n_dim_distribution_3 = len(token_list)

	average_postp_file_psiblast = open(average_postp_filename_psiblast, 'r')
	average_postp_file_psiblast_inside = average_postp_file_psiblast.readlines()
	average_postp_file_psiblast.close()
	n_lines_distribution_1 = len(average_postp_file_psiblast_inside)
	token_list = average_postp_file_psiblast_inside[0].rstrip().split()
	n_dim_distribution_1 = len(token_list)

        average_postp_file_hhmake = open(average_postp_filename_hhmake, 'r')
        average_postp_file_hhmake_inside = average_postp_file_hhmake.readlines()
        average_postp_file_hhmake.close()
        n_lines_distribution_2 = len(average_postp_file_hhmake_inside)
        token_list = average_postp_file_hhmake_inside[0].rstrip().split()
        n_dim_distribution_2 = len(token_list)

	n_dim_total = n_dim_psiblast + n_dim_hhmake + n_dim_distribution_1 + n_dim_distribution_2 + n_dim_distribution_3

	if ((n_lines_psiblast == n_lines_hhmake) and (n_lines_distribution_1 == n_lines_distribution_2) and (n_lines_distribution_2 == n_lines_distribution_3)):
		n_lines_equal_flag = 1
	else:
		n_lines_equal_flag = 0

	if (n_lines_equal_flag == 0):
		print 'Number of lines are not all the same for the psiblast, hhmake, dist 1, dist 2 and dist 3 files. Exiting...'
		sys.exit(1)

	seq_length = n_lines_psiblast
	half_window_size = (pssm_window_size_two_sided - 1) / 2

	svm_input_data_file = open(svm_input_data_filename, 'w')

	for m in range(seq_length):

		#token_list_psiblast_row_m = dat_file_psiblast_inside[m].rstrip().split()
		central_position_label = 0
		svm_input_data_file.write('%d ' % central_position_label)
		feature_index = 1

		for k in range(pssm_window_size_two_sided):

			pssm_row_selected = m + k - half_window_size

			if (pssm_row_selected < 0) or (pssm_row_selected > seq_length - 1):
				for jj in range(n_dim_total):
					svm_input_data_file.write('%d:%f ' % (feature_index, 0.0))
                                        feature_index += 1
			else:
				token_list_psiblast_row = dat_file_psiblast_inside[pssm_row_selected].rstrip().split()
                                token_list_hhmake_row = dat_file_hhmake_inside[pssm_row_selected].rstrip().split()
                                token_list_distribution_3_row = average_distribution_3_file_inside[pssm_row_selected].rstrip().split()
                                token_list_distribution_1_row = average_postp_file_psiblast_inside[pssm_row_selected].rstrip().split()
                                token_list_distribution_2_row = average_postp_file_hhmake_inside[pssm_row_selected].rstrip().split()
				
				for jj in range(n_dim_psiblast):
                                	svm_input_data_file.write('%d:%f ' % (feature_index, float(token_list_psiblast_row[jj])))
					feature_index += 1

                                for jj in range(n_dim_hhmake):
                                        svm_input_data_file.write('%d:%f ' % (feature_index, float(token_list_hhmake_row[jj])))
                                        feature_index += 1

                                for jj in range(n_dim_distribution_3):
                                        svm_input_data_file.write('%d:%f ' % (feature_index, float(token_list_distribution_3_row[jj])))
                                        feature_index += 1

                                for jj in range(n_dim_distribution_1):
                                        svm_input_data_file.write('%d:%f ' % (feature_index, float(token_list_distribution_1_row[jj])))                                   
                                        feature_index += 1

                                for jj in range(n_dim_distribution_2):
                                        svm_input_data_file.write('%d:%f ' % (feature_index, float(token_list_distribution_2_row[jj])))                                   
                                        feature_index += 1

		svm_input_data_file.write('\n')

	svm_input_data_file.close()

