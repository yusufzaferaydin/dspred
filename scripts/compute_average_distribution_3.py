__author__ = "Zafer Aydin"
__date__ = "01 July 2018"
__project__ = "Protein Structure Prediction"

usage = """
	python compute_average_distribution_3.py <prediction_task> <n_states> <psiblast_dbn_filename> <hhmake_dbn_filename> <sp_1_filename> <psiblast_weight> <psiblast_weight> <hhmake_weight> <sp_1_weight> <average_dbn_filename>
"""

alphabet_ss_3_state = "HEL"
alphabet_ss_8_state = "HGIEBSTL"
alphabet_torsion_5_state = "ABEGO"
alphabet_torsion_7_state = "LAMBEGO"
alphabet_sa_2_state = "eb"

import sys
import os

def give_me_ss_label(ss_index):

	ss_label = (ss_index == 0)*'H' + (ss_index == 1)*'E' + (ss_index == 2)*'L'
	return ss_label

def give_me_torsion_label(torsion_index):

	torsion_label = (torsion_index == 0)*'A' + (torsion_index == 1)*'B' + (torsion_index == 2)*'E' + (torsion_index == 3)*'G' + (torsion_index == 4)*'O'
	return torsion_label

if __name__ == "__main__":

	if (len(sys.argv) != 12):
		sys.stderr.write(usage)
		sys.exit(1)

	prediction_task = int(sys.argv[1])
	n_states = int(sys.argv[2])
	psiblast_dbn_filename = sys.argv[3]
        hhmake_dbn_filename = sys.argv[4]
        sp_1_filename = sys.argv[5]
        sp_2_filename = sys.argv[6]
	psiblast_weight = float(sys.argv[7])
	hhmake_weight = float(sys.argv[8])
        sp_1_weight = float(sys.argv[9])
        sp_2_weight = float(sys.argv[10])
	average_dbn_filename = sys.argv[11]

	psiblast_dbn_file = open(psiblast_dbn_filename, 'r')
	psiblast_dbn_file_inside = psiblast_dbn_file.readlines()
	psiblast_dbn_file.close()
	n_lines_psiblast = len(psiblast_dbn_file_inside)

	hhmake_dbn_file = open(hhmake_dbn_filename, 'r')
        hhmake_dbn_file_inside = hhmake_dbn_file.readlines()
        hhmake_dbn_file.close()
	n_lines_hhmake = len(hhmake_dbn_file_inside)

        sp_1_file = open(sp_1_filename, 'r')
        sp_1_file_inside = sp_1_file.readlines()
        sp_1_file.close()
        n_lines_sp_1 = len(sp_1_file_inside)

	if (os.path.exists(sp_2_filename)):
	        sp_2_file = open(sp_2_filename, 'r')
        	sp_2_file_inside = sp_2_file.readlines()
	        sp_2_file.close()
        	n_lines_sp_2 = len(sp_2_file_inside)
		sp_2_file_exists_flag = 1
	else:
		sp_2_weight = 0.0
		sp_2_file_exists_flag = 0

	average_dbn_file = open(average_dbn_filename, 'w')
	
	if (n_lines_psiblast != n_lines_hhmake):
		print 'Error: The DBN out files from psiblast and hhmake profiles have different number of amino acids'
		print psiblast_dbn_filename
		print hhmake_dbn_filename
		sys.exit(1)

        if (n_lines_psiblast != n_lines_sp_1):
                print 'Warning: The DBN out files from psiblast and sp_1 profiles have different number of amino acids'
                print psiblast_dbn_filename
                print sp_1_filename
                sys.exit(1)

        if (sp_2_file_exists_flag == 1) and (n_lines_psiblast != n_lines_sp_2):
                print 'Warning: The DBN out files from psiblast and sp_2 profiles have different number of amino acids'
                print psiblast_dbn_filename
                print sp_2_filename
                sys.exit(1)

	#print psiblast_weight
	#print hhmake_weight
	#print sp_1_weight

        psiblast_weight_orig = psiblast_weight
        hhmake_weight_orig = hhmake_weight
        sp_1_weight_orig = sp_1_weight
        sp_2_weight_orig = sp_2_weight

	n_aas = 0
	n_aas_no_hit = 0

	for i in range(n_lines_psiblast):
		
		line = psiblast_dbn_file_inside[i].rstrip()
		token_list = line.split()
		backward_index = -1*n_states
		posteriors_psiblast = token_list[backward_index:]

                line = hhmake_dbn_file_inside[i].rstrip()
                token_list = line.split()
                posteriors_hhmake = token_list[backward_index:]

                psiblast_weight = psiblast_weight_orig
                hhmake_weight = hhmake_weight_orig
                sp_1_weight = sp_1_weight_orig
                sp_2_weight = sp_2_weight_orig

		if (n_lines_sp_1 > 0):
	                line = sp_1_file_inside[i].rstrip()
        	        token_list = line.split()
                	posteriors_sp_1 = token_list[backward_index:]

                if (n_lines_sp_2 > 0):
                        line = sp_2_file_inside[i].rstrip()
                        token_list_2 = line.split()
                        posteriors_sp_2 = token_list_2[backward_index:]

		posteriors_average = [0]*n_states
		max_posterior = 0.0
		max_label_index = 0

                if (n_lines_sp_1 > 0):

                        no_hit_aa_flag = 0
                        for k in range(n_states-1):
                                if (posteriors_sp_1[k] == posteriors_sp_1[k+1]):
                                        no_hit_aa_flag = 1
                                else:
                                        no_hit_aa_flag = 0
                                        break
                else:
                        no_hit_aa_flag = 1

                if (n_lines_sp_2 > 0):

                        no_hit_aa_flag_2 = 0
                        for k in range(n_states-1):
                                if (posteriors_sp_2[k] == posteriors_sp_2[k+1]):
                                        no_hit_aa_flag_2 = 1
                                else:
                                        no_hit_aa_flag_2 = 0
                                        break
                else:
                        no_hit_aa_flag_2 = 1

                if (no_hit_aa_flag == 1):
                        sp_1_weight = 0.0

                if (no_hit_aa_flag_2 == 1):
                        sp_2_weight = 0.0

                if (sp_1_weight == 0.0) and (sp_2_weight == 0.0):
                        psiblast_weight = 0.5
                        hhmake_weight = 0.5
			n_aas_no_hit += 1

                total_weight = psiblast_weight + hhmake_weight + sp_1_weight + sp_2_weight

		for j in range(n_states):

                        if (total_weight != 0.0):

				posteriors_average[j] = (psiblast_weight*float(posteriors_psiblast[j]) + hhmake_weight*float(posteriors_hhmake[j]) + sp_1_weight*float(posteriors_sp_1[j]) + sp_2_weight*float(posteriors_sp_2[j])) / total_weight
			else:
                                posteriors_average[j] = 0.0
			

                        if (posteriors_average[j] < 0.0):
                                print "Posterior average is negative for %s" % psiblast_dbn_filename
                                print "Posterior average is negative for %s" % hhmake_dbn_filename
                                print "Posterior average is negative for %s" % struct_1_dbn_filename
                                print "Posterior average is negative for %s" % struct_2_dbn_filename
                                print psiblast_weight
                                print hhmake_weight
                                print sp_1_weight
                                print sp_2_weight
                                print posteriors_psiblast[j]
                                print posteriors_hhmake[j]
                                print posteriors_sp_1[j]
                                print posteriors_sp_2[j]
                                sys.exit(1)

                        if ( round(posteriors_average[j],5) > round(max_posterior,5) ):
                                max_posterior = posteriors_average[j]
                                max_label_index = j

                n_aas += 1

		if (prediction_task == 1):
			pred_label = alphabet_ss_3_state[max_label_index]

		elif (prediction_task == 2):
                        pred_label = alphabet_ss_8_state[max_label_index]

		elif (prediction_task == 3):
                        pred_label = alphabet_torsion_5_state[max_label_index]

                elif (prediction_task == 4):
                        pred_label = alphabet_torsion_7_state[max_label_index]

                elif (prediction_task == 5):
                        pred_label = alphabet_sa_2_state[max_label_index]

		#average_dbn_file.write('%d\t%s\t%s\t' % (position_index, pred_label))
		
		for j in range(n_states):
			average_dbn_file.write('%f\t' % posteriors_average[j]) 
		
		average_dbn_file.write('\n')

        average_dbn_file.close()
	
	percentage_aas_skipped = 100*float(n_aas_no_hit) / n_aas
	#print 'percentage_aas_skipped = %f' % percentage_aas_skipped

