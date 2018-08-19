__author__ = "Zafer Aydin"
__date__ = "11/05/2009"
__project__ = "Protein Structure Prediction"

usage = """
	Usage: python compute_average_postp_files.py <number_of_state_classes> <dataset_dir> <fastas all list filename> <obs files list filename> <postp files list filename forward model> <postp files list filename reverse_model> <pssm window size>
"""

import sys

ss_alphabet_3 = "HEL"
ss_alphabet_5 = "HEebL"
ss_alphabet_8 = "HBEGITSLLL"   #Label 9 and 10 are also L in SS labels
torsion_alphabet_5 = "ORGAB"
joint_alphabet_7 = "ABCDEFG"
joint_alphabet_11 = "ABCDEFGHIJK"

def generate_true_ss_file(prediction_task, data_dir, obs_files_list_filename, pssm_window_size, true_ss_filename):

	obs_files_list_file = open( obs_files_list_filename, 'r' )
	true_ss_file = open( true_ss_filename, 'w' )

	for line in obs_files_list_file:

		obs_filename = line.rstrip()
		obs_file = open( obs_filename, 'r' )
		
		true_ss_labels_str = ''
		line_counter = 0

		for obs_line in obs_file:
			
			if (line_counter < pssm_window_size):
				line_counter += 1
				continue
			
			token_list = obs_line.split()	
			#ss_label_int = int( token_list[-1] )
			ss_label_int = int( token_list[-3] )

			if (prediction_task == 1):

				true_ss_labels_str += ss_alphabet_3[ss_label_int]  
		
			elif (prediction_task == 2):

				true_ss_labels_str += ss_alphabet_5[ss_label_int]

			elif (prediction_task == 3):
			
				true_ss_labels_str += torsion_alphabet_5[ss_label_int]
	
		obs_file.close()

		protein_id = obs_filename.replace(data_dir, '')
		protein_id = protein_id.replace('.obs', '')

		true_ss_file.write( '>' + protein_id + '\n')
		true_ss_file.write(true_ss_labels_str + '\n');

	true_ss_file.close()
	obs_files_list_file.close()


def generate_true_and_predicted_ss_file( prediction_task, number_of_state_classes, data_dir, pssm_window_size, fasta_filenames_list_filename, obs_files_list_filename, postp_files_list_filename_forward, postp_files_list_filename_reverse ):

        obs_files_list_file = open( obs_files_list_filename, 'r' )
	postp_files_list_file_forward = open( postp_files_list_filename_forward, 'r' )
	postp_files_list_file_reverse = open( postp_files_list_filename_reverse, 'r' )
	fasta_filenames_file = open( fasta_filenames_list_filename, 'r')

	done  = 0

	while not done:
		
		line_forward = postp_files_list_file_forward.readline()
		line_reverse = postp_files_list_file_reverse.readline()
	
		if (( len( line_forward ) == 0 ) or ( len( line_reverse ) == 0 )):
		#if ( len( line_forward ) == 0 ):

			done = 1
			continue

                line_fasta_filename = fasta_filenames_file.readline().rstrip()
                fasta_filename = data_dir + line_fasta_filename
                fasta_file = open(fasta_filename, 'r')
                fasta_id = fasta_file.readline()
                aa_seq = fasta_file.readline().rstrip()
                fasta_file.close()

		postp_filename_forward = line_forward.rstrip()
		postp_filename_reverse = line_reverse.rstrip()

		postp_file_forward = open( postp_filename_forward, 'r' )
		postp_file_reverse = open( postp_filename_reverse, 'r' )

                if ('nc.aapw.test.postp.sigmoid' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.aapw.test.postp.sigmoid', 'ss3')

                elif ('cn.aapw.test.postp.sigmoid' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.aapw.test.postp.sigmoid', 'ss3')

                elif ('nc.dcw.test.postp.sigmoid' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.dcw.test.postp.sigmoid', 'ss3')

                elif ('cn.dcw.test.postp.sigmoid' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.dcw.test.postp.sigmoid', 'ss3')

                elif ('nc.test.postp.sigmoid' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.test.postp.sigmoid', 'ss3')

                elif ('cn.test.postp.sigmoid' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.test.postp.sigmoid', 'ss3')

                elif ('nc.aapw.test.postp.linear' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.aapw.test.postp.linear', 'ss3')

                elif ('cn.aapw.test.postp.linear' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.aapw.test.postp.linear', 'ss3')

                elif ('nc.dcw.test.postp.linear' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.dcw.test.postp.linear', 'ss3')

                elif ('cn.dcw.test.postp.linear' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.dcw.test.postp.linear', 'ss3')

                elif ('nc.test.postp.linear' in postp_filename_forward):

                        postp_file_out_filename = postp_filename_forward.replace('nc.test.postp.linear', 'ss3')

                elif ('cn.test.postp.linear' in postp_filename_reverse):

                        postp_file_out_filename = postp_filename_reverse.replace('cn.test.postp.linear', 'ss3')

		#postp_file_out_filename = postp_filename_forward.replace('.test.postp.sigmoid', '.out')
		#postp_file_out_filename = postp_file_out_filename.replace('.test.postp.linear', '.out')
		postp_file_out = open( postp_file_out_filename, 'w')
		postp_file_out.write('%s\n\n' % '# DBNPRED_v1.0 3-State Secondary Structure Prediction H:Helix, E:strand, L:loop')

		pred_ss_labels_str = ''

		line_counter = 1

		postp_file_forward_inside =  postp_file_forward.readlines()
                postp_file_reverse_inside =  postp_file_reverse.readlines()

		postp_file_forward_inside_length = len(postp_file_forward_inside)
		postp_file_reverse_inside_length = len(postp_file_reverse_inside)

		if (postp_file_forward_inside_length != postp_file_reverse_inside_length):
			
			print postp_file_forward_inside_length
			print postp_file_reverse_inside_length
			print 'Error the length of the postp files for the forward and reverse models are not the same'
			print line_fasta_filename
			sys.exit(-1)


		#Read the TRUE labels 
		obs_filename = obs_files_list_file.readline().rstrip()
                obs_file = open( obs_filename, 'r' )

                true_labels = ''
                line_counter = 0

                for obs_line in obs_file:

                        if (line_counter < pssm_window_size):
                                line_counter += 1
                                continue

                        token_list = obs_line.split()
                        #ss_label_int = int( token_list[-1] )
                        #ss_label_int = int( token_list[-3] )

                        #print ss_label_int

			#if (prediction_task == 1):

                        #        true_labels += ss_alphabet_3[ss_label_int]

                        #elif (prediction_task == 2):

                        #        true_labels += ss_alphabet_5[ss_label_int]

                        #elif (prediction_task == 3):

                        #        true_labels += torsion_alphabet_5[ss_label_int]

                        #elif (prediction_task == 4):

                        #        true_labels += joint_alphabet_7[ss_label_int]

                        #elif (prediction_task == 5):

                        #        true_labels += joint_alphabet_11[ss_label_int]

                obs_file.close()

		for i in range(0, postp_file_forward_inside_length, 1):

			postp_line_forward = postp_file_forward_inside[i]
			postp_line_reverse = postp_file_reverse_inside[-1-i]

			token_list_forward = postp_line_forward.split()
			token_list_reverse = postp_line_reverse.split()
			
			#print '%s' % token_list_forward
			#print '%s' % token_list_reverse

			#print len( token_list_forward )
			#print number_of_state_classes

			if ( len( token_list_forward ) == number_of_state_classes ):

				max_value = -1000
				max_index = 0			

				postp_values = []
	
				for j in range(0, len( token_list_forward ), 1):
				
					#Postp value is the average of postp values from the forward and reverse models
					postp_value = ( float( token_list_forward[j] ) + float( token_list_reverse[j] ) ) / 2.0
					#postp_value = float(token_list_forward[j])
					#postp_value = float(token_list_reverse[j])

					postp_values += [ postp_value ]

					#print '%f' % float(token_list_forward[j])
					#print '%f' % float(token_list_reverse[j])
					#print '%f' % postp_value

					if ( postp_value > max_value):
						
						max_value = postp_value
						max_index = j


				#print postp_values
				
				#print '%s' % ssalphabet[max_index]

				if (prediction_task == 1):
			 	
					pred_label = ss_alphabet_3[max_index]
					pred_ss_labels_str += pred_label
					
				elif (prediction_task == 2):

					pred_label = ss_alphabet_5[max_index]
	                                pred_ss_labels_str += pred_label
	
				elif (prediction_task == 3):

					pred_label = torsion_alphabet_5[max_index]
					pred_ss_labels_str += pred_label

                                elif (prediction_task == 4):

                                        pred_label = joint_alphabet_7[max_index]
                                        pred_ss_labels_str += pred_label

                                elif (prediction_task == 5):

                                        pred_label = joint_alphabet_11[max_index]
                                        pred_ss_labels_str += pred_label
			
				aa_type = aa_seq[i]
				#true_label = true_labels[i]

				if (pred_label == 'R'):
					pred_label = 'E'

				kk=i+1
                                postp_file_out.write("%d\t%s\t%s\t" % (kk, aa_type, pred_label))

				for j in range(len(token_list_forward)):

					postp_file_out.write("%f\t" % postp_values[j])

				postp_file_out.write("\n")

				#print max_index	
				#print pred_ss_labels_str

			#print pred_ss_labels_str
			line_counter += 1

		postp_file_forward.close()
		postp_file_reverse.close()
		postp_file_out.close()

		#print '%s' % pred_ss_labels_str

		#protein_id = postp_filename_forward.replace(data_dir, '')
		#protein_id = protein_id.replace('.w_test.postp.sigmoid', '')
		#protein_id = protein_id.replace('.test.postp.sigmoid', '')
		#protein_id = protein_id.replace('.test.postp.linear', '')
		#protein_id = protein_id.replace('.test.postp', '')

                #pred_ss_file.write( '>' + protein_id + '\n' )
                #pred_ss_file.write( pred_ss_labels_str + '\n' );

	obs_files_list_file.close()
	postp_files_list_file_forward.close()
	postp_files_list_file_reverse.close()
	fasta_filenames_file.close()


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# main program is here

import sys
import math

if __name__ == "__main__":
    
	if ( len(sys.argv) != 9 ):
		
		sys.stderr.write(usage)
		sys.exit(-1)
	
	prediction_task = int(sys.argv[1])
	number_of_state_classes = int(sys.argv[2])
	data_dir = sys.argv[3]
	fasta_filenames_list_filename = sys.argv[4]
	obs_files_list_filename = sys.argv[5]
	postp_files_list_filename_forward = sys.argv[6]
	postp_files_list_filename_reverse = sys.argv[7]
	pssm_window_size = int(sys.argv[8])

	#generate_true_ss_file( number_of_state_classes, data_dir, pssm_window_size )
	generate_true_and_predicted_ss_file( prediction_task, number_of_state_classes, data_dir, pssm_window_size, fasta_filenames_list_filename, obs_files_list_filename, postp_files_list_filename_forward, postp_files_list_filename_reverse ) 
	


 


