__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "02 July 2018"

usage = """
	Usage: python reformat_svm_prediction_output.py <prediction_task> <libsvm_output_filename> <dspred_output_filename> 
"""

ss_labels_3_state = 'HEL'
ss_labels_8_state = 'HGIEBST ' #this is the DSSP convention
torsion_labels_5_state = 'ABEGO'
torsion_labels_7_state = 'LAMBEGO'
sa_labels_2_state = 'eb' #meaning exposed or buried

import sys

if __name__ == "__main__":

        if (len(sys.argv) != 4):
                sys.stderr.write(usage)
                sys.exit(1)

	prediction_task = int(sys.argv[1])
	libsvm_output_filename = sys.argv[2]
	dspred_output_filename = sys.argv[3]

	libsvm_output_file = open(libsvm_output_filename, 'r')
	libsvm_output_file_inside = libsvm_output_file.readlines()
	libsvm_output_file.close()
	n_lines = len(libsvm_output_file_inside)
	n_aas = n_lines-1
	label_order_line = libsvm_output_file_inside[0].rstrip()
        token_list = label_order_line.split()
	n_classes = len(token_list)-1	
        labels_libsvm_order = [0]*n_classes
        for i in range(n_classes):
                labels_libsvm_order[i] = int(token_list[i+1])

        prediction_outputs_aa = [0.0]*n_classes

	dspred_output_file = open(dspred_output_filename, 'w')

	for i in range(n_aas):

		prediction_line = libsvm_output_file_inside[i+1].rstrip()
		token_list = prediction_line.split()

		predicted_label = int(token_list[0])

		if (prediction_task == 1):
			predicted_label_str = ss_labels_3_state[predicted_label]

		elif (prediction_task == 2):
			predicted_label_str = ss_labels_8_state[predicted_label]

		elif (prediction_task == 3):
			predicted_label_str = torsion_labels_5_state[predicted_label]

		elif (prediction_task == 4):
			predicted_label_str = torsion_labels_7_state[predicted_label]

		elif (prediction_task == 5):
			predicted_label_str = sa_labels_2_state[predicted_label]

                for k in range(n_classes):
                	prediction_outputs_aa[labels_libsvm_order[k]] = float(token_list[-1*n_classes+k])		

		dspred_output_file.write('%d\t%s\t' % ((i+1), predicted_label_str))

                for k in range(n_classes):
			dspred_output_file.write('%.5f\t' % prediction_outputs_aa[k])
		dspred_output_file.write('\n')

	dspred_output_file.close()

