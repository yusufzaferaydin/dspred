__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "29 June 2018"

usage = """
	Usage: python compute_average_postp_file.py <postp_filename_nc> <postp_filename_cn> <average_postp_filename>
"""

import sys

if __name__ == "__main__":

	if (len(sys.argv) != 4):
		sys.stderr.write(usage)
		sys.exit(1)

	postp_filename_nc = sys.argv[1]
        postp_filename_cn = sys.argv[2]
	average_postp_filename = sys.argv[3]

	postp_file_nc = open(postp_filename_nc, 'r')
	postp_file_nc_inside = postp_file_nc.readlines()
	postp_file_nc.close()
	n_lines = len(postp_file_nc_inside)

        postp_file_cn = open(postp_filename_cn, 'r')
        postp_file_cn_inside = postp_file_cn.readlines()
        postp_file_cn.close()
        n_lines_2 = len(postp_file_cn_inside)
	
	if (n_lines != n_lines_2):
		print 'Number of amino acids is not the same in postp files of NC and CN models. Exiting...'
		sys.exit(1)

	average_postp_file = open(average_postp_filename, 'w')

	for i in range(n_lines):

		prediction_score_aa_nc = postp_file_nc_inside[i].rstrip()
		token_list_nc = prediction_score_aa_nc.split()
		prediction_score_aa_cn = postp_file_cn_inside[i].rstrip()
                token_list_cn = prediction_score_aa_cn.split()
		
		n_classes_nc = len(token_list_nc)
                n_classes_cn = len(token_list_cn)

		if (n_classes_nc != n_classes_cn):
			print 'Number of classes are not the same in postp files of NC and CN models. Exiting...'
			sys.exit(1)

		for j in range(n_classes_nc):
			average_prediction_score = (float(token_list_nc[j]) + float(token_list_cn[j])) / 2.0
			average_postp_file.write('%.5f ' % average_prediction_score)
		average_postp_file.write('\n')			

	average_postp_file.close()

