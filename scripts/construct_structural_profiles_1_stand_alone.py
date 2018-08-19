__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "08/10/2015"

usage = """

	Usage: python construct_structural_profiles_1.py prediction_task $hmm_alignments_dir $hmm_alignments_list_filename $get_pdb_script_filename $compute_torsion_angles_script_filename $rosetta_database_dir $hhblits_score_percentage_threshold $power_parameter $blosum_matrix_filename $output_filename 
"""

import sys
import os
import math
import string

alphabet_ss = "HEL"
alphabet_torsion_5 = "ABEGO"
alphabet_torsion_7 = "LAMBEGO"

#the following max accessbility is obtained from https://en.wikipedia.org/wiki/Relative_accessible_surface_area
max_accessibility = {'A': 106, 'B': 160, 'C': 135, 'D': 163, 'E': 194, 'F': 197, 'G': 84, 'H': 184, 'I': 169, 'K': 205, 'L': 164, 'M': 188, 'N': 157, 'P': 136, 'Q': 198, 'R': 248, 'S': 130, 'T': 142, 'V': 142, 'W': 227, 'X': 180, 'Y': 222, 'Z': 196}

#three letter codes are available at https://en.wikipedia.org/?title=Amino_acid
three_letter_aa_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J', 'Xaa': 'X', 'XAA': 'X', 'UNK': 'X', 'Unk': 'X', 'SEC': 'U', 'Pyl': 'O'}

def give_me_aa_number( aa ):

        aa_labels = 'ARNDCQEGHILKMFPSTWYV'
        aa_number = aa_labels.find(aa)

        if ( aa_number == -1 ):

                return 0
        else:
                return aa_number

def extract_aa_sequence(torsion_filename):

        aa_sequence = ''

        torsion_file = open(torsion_filename, 'r')
        torsion_file_inside = torsion_file.readlines()
        torsion_file.close()
        n_lines = len(torsion_file_inside)

        phrase_found_flag = 0
        for i in range(n_lines):
                line = torsion_file_inside[i].rstrip()
                if ('idx res ss      phi      psi    omega' in line):
                        torsion_labels_starting_index = i+1
                        phrase_found_flag = 1
                        break

        if (phrase_found_flag == 0):
                return aa_sequence

        #print torsion_labels_starting_index

        done = 0
        while not done:
                line = torsion_file_inside[torsion_labels_starting_index].rstrip()
                token_list = line.split()
                position_index = token_list[0]
                if (position_index != str(1)):
                        torsion_labels_starting_index += 1
                else:
                        done = 1

                if (torsion_labels_starting_index >= n_lines):
                        done = 1

        for i in range(torsion_labels_starting_index, n_lines, 1):

                line = torsion_file_inside[i].rstrip()
                token_list = line.split()
                aa_type = token_list[1]
                aa_sequence += aa_type

        return aa_sequence

def extract_torsion_label_sequence(torsion_filename):

        torsion_label_sequence = ''

        torsion_file = open(torsion_filename, 'r')
        torsion_file_inside = torsion_file.readlines()
        torsion_file.close()
        n_lines = len(torsion_file_inside)

        phrase_found_flag = 0
        for i in range(n_lines):
                line = torsion_file_inside[i].rstrip()
                if ('idx res ss      phi      psi    omega' in line):
                        torsion_labels_starting_index = i+1
                        phrase_found_flag = 1
                        break

        if (phrase_found_flag == 0):
                return torsion_label_sequence

        #print torsion_labels_starting_index

        done = 0
        while not done:
                line = torsion_file_inside[torsion_labels_starting_index].rstrip()
                token_list = line.split()
                position_index = token_list[0]
                if (position_index != str(1)):
                        torsion_labels_starting_index += 1
                else:
                        done = 1

                if (torsion_labels_starting_index >= n_lines):
                        done = 1

        for i in range(torsion_labels_starting_index, n_lines, 1):

                line = torsion_file_inside[i].rstrip()
                token_list = line.split()
                torsion_label_ABEGO = token_list[10]
                ss_label = token_list[2]

                if ((torsion_label_ABEGO == 'A') and (ss_label == 'L')):
                        torsion_label_LAMBEGO = 'L'
                elif ((torsion_label_ABEGO == 'A') and (ss_label != 'L')):
                        torsion_label_LAMBEGO = 'A'
                elif ((torsion_label_ABEGO == 'B') and (ss_label == 'L')):
                        torsion_label_LAMBEGO = 'M'
                elif ((torsion_label_ABEGO == 'B') and (ss_label != 'L')):
                        torsion_label_LAMBEGO = 'B'
                else:
                        torsion_label_LAMBEGO = torsion_label_ABEGO

                torsion_label_sequence += torsion_label_LAMBEGO

        return torsion_label_sequence

def retrieve_phi_psi_omega(torsion_angles_filename):

        torsion_angles_file = open(torsion_angles_filename, 'r')
        torsion_angles_file_inside = torsion_angles_file.readlines()
        torsion_angles_file.close()
        n_lines = len(torsion_angles_file_inside)

	#print torsion_angles_filename
	#print torsion_angles_file_inside

        for i in range(n_lines):
                line = torsion_angles_file_inside[i].rstrip()
                if ('NBRES' in line):
                        token_list = line.split()
                        seq_length = int(token_list[1])
                        break

        phi_psi_omegas = [0]*seq_length
        for i in range(seq_length):
                phi_psi_omegas[i] = [0]*3

        aa_seq_pdb = ''
        line_shift = n_lines - seq_length
        for i in range(seq_length):
                line = torsion_angles_file_inside[i+line_shift].rstrip()
                token_list = line.split()
                aa_type = token_list[0]
                aa_seq_pdb += three_letter_aa_code[aa_type]
                phi = float(token_list[2])
                psi = float(token_list[3])
                omega = float(token_list[4])

                if (phi == 500) or (phi == 9999):
                        phi = 0.0

                if (psi == 500) or (psi == 9999):
                        psi = 0.0

                if (omega == 500) or (omega == 9999):
                        omega = 0.0

                phi_psi_omegas[i][0] = phi
                phi_psi_omegas[i][1] = psi
                phi_psi_omegas[i][2] = omega

        return (phi_psi_omegas, aa_seq_pdb)

def reduce_from_8_to_3_state_EHL_mapping( ss_seq ):

        str_length = len(ss_seq)

        ss_seq_red = ''

        for i in range(0, str_length, 1):

                if ( (ss_seq[i] == 'H') or (ss_seq[i] == 'G') or (ss_seq[i] == 'I') ):

                        ss_seq_red += 'H'

                elif ( (ss_seq[i] == 'E') or (ss_seq[i] == 'B') or (ss_seq[i] == 'e') or (ss_seq[i] == 'b')):

                        ss_seq_red += 'E'

                elif ( (ss_seq[i] == 'S') or (ss_seq[i] == 'T') or (ss_seq[i] == ' ') or (ss_seq[i] == 'L')):

                        ss_seq_red += 'L'

        return ss_seq_red

def convert_sa_to_2_state(sa_seq, aa_seq):

        token_list = sa_seq.split()
        seq_length = len(token_list)

        sa_seq_2 = '' #this one has solvent accessibility labels

        for i in range(seq_length):

                aa = aa_seq[i]
                max_acc = float(max_accessibility[aa])
                absolute_acc = float(token_list[i])
                relative_acc = 100*absolute_acc / max_acc

                if (relative_acc < 25):
                        sa_seq_2 += 'b'
                else:
                        sa_seq_2 += 'e'

        return sa_seq_2

def convert_torsion_angles_to_5_state(phi_psi_omegas, aa_seq, ss_seq_3_state, aa_seq_2, exclamation_present_flag, dssp_seq_in_pdb_seq_flag):

        seq_length = len(phi_psi_omegas)
        seq_length_2 = len(aa_seq)

        torsion_seq = ''

        if (dssp_seq_in_pdb_seq_flag == 1): #only covers cases where the dssp sequence is shorter than the pdb sequence
                starting_index = aa_seq_2.find(aa_seq)
        else:
                starting_index = 0

        if (exclamation_present_flag == 1):

                exclamation_indices = []
                for i in range(seq_length):
                        if (aa_seq_2[i] == '!'):
                                exclamation_indices += [i]

                n_exclamation_indices = len(exclamation_indices)

                for i in range(seq_length):

                        if (aa_seq_2[i] == '!'): #chain breaks or discontuity in sequence
                                continue

                        phi = phi_psi_omegas[i][0]
                        psi = phi_psi_omegas[i][1]
                        omega = phi_psi_omegas[i][2]

                        n_exclamation_indices_less_than_i = 0
                        for j in range(n_exclamation_indices):
                                if (exclamation_indices[j] < i):
                                        n_exclamation_indices_less_than_i += 1

                        k = n_exclamation_indices_less_than_i

                        if ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i-k] == 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i-k] != 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i-k] == 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i-k] != 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) > 100)):
                                torsion_seq += 'E'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) <= 100)):
                                torsion_seq += 'G'
                        elif (abs(omega) < 90):
                                torsion_seq += 'O'

        elif (exclamation_present_flag == 0):
                for i in range(seq_length_2):

                        phi = phi_psi_omegas[i+starting_index][0]
                        psi = phi_psi_omegas[i+starting_index][1]
                        omega = phi_psi_omegas[i+starting_index][2]

                        if ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i] == 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i] != 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i] == 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i] != 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) > 100)):
                                torsion_seq += 'E'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) <= 100)):
                                torsion_seq += 'G'
                        elif (abs(omega) < 90):
                                torsion_seq += 'O'

        return torsion_seq

def convert_torsion_angles_to_7_state(phi_psi_omegas, aa_seq, ss_seq_3_state, aa_seq_2, exclamation_present_flag, dssp_seq_in_pdb_seq_flag):

        seq_length = len(phi_psi_omegas)
        seq_length_2 = len(aa_seq)

        torsion_seq = ''

        if (dssp_seq_in_pdb_seq_flag == 1): #only covers cases where the dssp sequence is shorter than the pdb sequence
                starting_index = aa_seq_2.find(aa_seq)
        else:
                starting_index = 0

        if (exclamation_present_flag == 1):

                exclamation_indices = []
                for i in range(seq_length):
                        if (aa_seq_2[i] == '!'):
                                exclamation_indices += [i]

                n_exclamation_indices = len(exclamation_indices)

                for i in range(seq_length):

                        if (aa_seq_2[i] == '!'): #chain breaks or discontuity in sequence
                                continue

                        phi = phi_psi_omegas[i][0]
                        psi = phi_psi_omegas[i][1]
                        omega = phi_psi_omegas[i][2]

                        n_exclamation_indices_less_than_i = 0
                        for j in range(n_exclamation_indices):
                                if (exclamation_indices[j] < i):
                                        n_exclamation_indices_less_than_i += 1

                        k = n_exclamation_indices_less_than_i

                        if ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i-k] == 'L')):
                                torsion_seq += 'L'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i-k] != 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i-k] == 'L')):
                                torsion_seq += 'M'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i-k] != 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) > 100)):
                                torsion_seq += 'E'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) <= 100)):
                                torsion_seq += 'G'
                        elif (abs(omega) < 90):
                                torsion_seq += 'O'

        elif (exclamation_present_flag == 0):
                for i in range(seq_length_2):

                        phi = phi_psi_omegas[i+starting_index][0]
                        psi = phi_psi_omegas[i+starting_index][1]
                        omega = phi_psi_omegas[i+starting_index][2]

                        if ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i] == 'L')):
                                torsion_seq += 'L'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi > -125) and (psi <= 50)) and (ss_seq_3_state[i] != 'L')):
                                torsion_seq += 'A'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i] == 'L')):
                                torsion_seq += 'M'
                        elif ((abs(omega) >= 90) and (phi < 0) and ((psi <= -125) or (psi > 50)) and (ss_seq_3_state[i] != 'L')):
                                torsion_seq += 'B'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) > 100)):
                                torsion_seq += 'E'
                        elif ((abs(omega) >= 90) and (phi >= 0) and (abs(psi) <= 100)):
                                torsion_seq += 'G'
                        elif (abs(omega) < 90):
                                torsion_seq += 'O'

        return torsion_seq

def convert_pdb_id_to_lowercase_format( pdb_id ):

        chain_id = pdb_id[-1]
        pdb_id_wo_chain = pdb_id[0:4]
        new_pdb_id = pdb_id_wo_chain.lower() + chain_id

        return new_pdb_id

def give_me_ss_number_3_state(ss_label):

        if ((ss_label != 'H') and (ss_label != 'E') and (ss_label != 'L')):
                print 'Error: ss label is not one of HEL labels'

        return ((ss_label == 'H')*0 + (ss_label == 'E')*1 + (ss_label == 'L')*2)

def give_me_ss_number_8_state(ss_label):

        if ((ss_label != 'H') and (ss_label != 'G') and (ss_label != 'I') and (ss_label != 'E') and (ss_label != 'B') and (ss_label != 'S') and (ss_label != 'T') and (ss_label != ' ') and (ss_label != 'L') and (ss_label != 'C')):
                print 'Error: ss label is not one of 8 state DSSP labels'

        return ((ss_label == 'H')*0 + (ss_label == 'G')*1 + (ss_label == 'I')*2 + ((ss_label == 'E') or (ss_label == 'e'))*3 + (ss_label == 'B')*4 + (ss_label == 'S')*5 + (ss_label == 'T')*6 + ((ss_label == ' ') or (ss_label == 'L') or (ss_label == 'C'))*7 )

def give_me_torsion_number_ABEGO(torsion_label):

        if ((torsion_label != 'A') and (torsion_label != 'B') and (torsion_label != 'E') and (torsion_label != 'G') and (torsion_label != 'O')):
                print 'Error: torsion label is not one of ABEGO labels'

        return ((torsion_label == 'A')*0 + (torsion_label == 'B')*1 + (torsion_label == 'E')*2 + (torsion_label == 'G')*3 + (torsion_label == 'O')*4)

def give_me_torsion_number_LAMBEGO(torsion_label):

        if ((torsion_label != 'L') and (torsion_label != 'A') and (torsion_label != 'M') and (torsion_label != 'B') and (torsion_label != 'E') and (torsion_label != 'G') and (torsion_label != 'O')):
                print 'Error: torsion label is not one of LAMBEGO labels'

        return ((torsion_label == 'L')*0 + (torsion_label == 'A')*1 + (torsion_label == 'M')*2 + (torsion_label == 'B')*3 + (torsion_label == 'E')*4 + (torsion_label == 'G')*5 + (torsion_label == 'O')*6)

def give_me_sa_number_2_state(sa_label):

        if ((sa_label != 'e') and (sa_label != 'b')):
                print 'Error: sa label is not one of eb labels'

        return ((sa_label == 'e')*0 + (sa_label == 'b')*1)

def read_blosum_matrix(blosum_matrix_filename):

        blosum_matrix_file = open(blosum_matrix_filename, 'r')

        done = 0

        blosum_matrix = [0]*20
        for i in range(20):
                blosum_matrix[i] = [0]*20

        i = 0
        while not done:

                line = blosum_matrix_file.readline()

                if (len(line) == 0):
                        done = 1
                        continue

                if (line[0] == '#'): #skip the comment lines
                        continue

                token_list = line.split()
                if (token_list[1] == 'R'):
                        continue

                for j in range(20):
                        blosum_matrix[i][j] = int(token_list[j+1])

                i += 1

                if (i >= 20):
                        done = 1
                        continue

        blosum_matrix_file.close()

        return blosum_matrix

if __name__ == "__main__":

	if (len(sys.argv) != 20):
		sys.stderr.write(usage)
		sys.exit(1)

	prediction_task = int(sys.argv[1])
	n_classes = int(sys.argv[2])
	fasta_filename = sys.argv[3]
	hmm_alignments_dir = sys.argv[4]
	pdb_files_dir = sys.argv[5]
	dssp_files_dir = sys.argv[6]
	torsion_files_dir = sys.argv[7]
	hmm_alignments_list_filename = sys.argv[8]
	get_pdb_script_filename = sys.argv[9]
        get_dssp_script_filename = sys.argv[10]
	extract_aa_ss_sa_labels_script_filename = sys.argv[11]
	compute_torsion_angles_script_filename = sys.argv[12]
	#rosetta_database_dir = sys.argv[10]
	hhblits_score_percentage_threshold = float(sys.argv[13])
	percentage_identity_threshold = float(sys.argv[14])
	power_parameter = int(sys.argv[15])
	blosum_matrix_filename = sys.argv[16]
	output_files_dir = sys.argv[17]
	max_sequence_identities_filename = sys.argv[18]
	sp_1_filename = sys.argv[19]

        hmm_alignments_list_file = open(hmm_alignments_list_filename, 'r')
        hmm_alignments_list_file_inside = hmm_alignments_list_file.readlines()
        hmm_alignments_list_file.close()
        n_alignments = len(hmm_alignments_list_file_inside)

        #blosum_matrix = read_blosum_matrix(blosum_matrix_filename)

	max_sequence_identities_file = open(max_sequence_identities_filename, 'w')

        protein_counter = 0
        for i in range(n_alignments):

		##assuming that the list file does not include the full path
                #hmm_alignment_filename = hmm_alignments_dir + hmm_alignments_list_file_inside[i].rstrip()

                #assuming that the list file includes the full path
                hmm_alignment_filename = hmm_alignments_list_file_inside[i].rstrip()
		
                hmm_alignment_file = open(hmm_alignment_filename, 'r')
                hmm_alignment_file_inside = hmm_alignment_file.readlines()
                hmm_alignment_file.close()
                n_lines = len(hmm_alignment_file_inside)

                if (n_lines == 0):
                        continue

                query_line = hmm_alignment_file_inside[0].rstrip()
                token_list = query_line.split()
                query_id = token_list[1]

                query_id_orig = query_id

		query_id_core = query_id[0:4]
                query_id_core_small_case = query_id_core.lower()
		
                #if ('3ecsA' not in query_id):
                #        continue

		print 'query = %s' % query_id
	
		fasta_file = open(fasta_filename, 'r')
		fasta_id = fasta_file.readline().rstrip()
		aa_sequence_query = fasta_file.readline().rstrip()
		fasta_file.close()
		n_aas_query = len(aa_sequence_query)

                query_id = convert_pdb_id_to_lowercase_format(query_id)

                protein_counter += 1

                hit_line_start_index = 9
                done = 0
                j = 0
                hit_ids = []
                probability_scores = []
                hit_indices = []
                max_probability_score = 0.0
                while not done:

                        hit_line_index = hit_line_start_index + j
                        hit_line = hmm_alignment_file_inside[hit_line_index].rstrip()

                        if (len(hit_line) == 0):
                                done = 1
                                continue

                        token_list2 = hit_line.split()
                        hit_id = token_list2[1]
                        chain_id = hit_id[-1]
                        hit_id_raw = hit_id[0:4]
                        #probability_score = float(token_list2[2])
			probability_score = float(hit_line[34:40])
                        probability_scores += [ probability_score ]
                        hit_index = int(token_list2[0])
                        hit_indices += [ hit_index ]

                        if (probability_score > max_probability_score):
                                max_probability_score = probability_score

                        hit_ids.append([hit_id])
                        j += 1

                probability_score_threshold = max_probability_score - hhblits_score_percentage_threshold

                #scan through the alignments one by one
                n_hits = j
                alignments_start_index = hit_line_start_index + n_hits + 1
                n_hits_processed = 0

                pssm = [1.0]*n_aas_query
                average_confidence_scores = [0.0]*n_aas_query
                position_specific_hit_counts = [0]*n_aas_query
		pssm_update_flag = 0

                for j in range(n_aas_query):
                        pssm[j] = [1.0]*n_classes
			average_confidence_scores[j] = [0.0]*n_classes
	                position_specific_hit_counts[j] = [0]*n_classes

                indices_covered = [0]*n_aas_query
		
		max_seq_identity = 0

                for j in range(n_hits):

                        hit_id = hit_ids[j][0]

                        hit_id_orig = hit_id

                        string_to_check_1 = 'No ' + str(hit_indices[j])
                        string_to_check_2 = '>' + hit_id_orig
                        string_to_check_3 = 'Probab='

                        string_to_check_4 = 'No '
                        string_to_check_5 = '>'
                        string_to_check_6 = 'Done!'

	                hit_id_core = hit_id[0:4]
			hit_id_core_small_case = hit_id_core.lower()
			chain_id = hit_id[-1]
	                #hit_id = hit_id_core.lower() + chain_id

                        hit_id = convert_pdb_id_to_lowercase_format(hit_id)

                        #if (hit_id in query_id) or (hit_id_core_small_case in query_id_core_small_case):
                        #        continue
		
			#print hit_id

                        pdb_filename = hit_id + '.pdb'
                        pdb_filename_hit = pdb_files_dir + hit_id + '.pdb'
                        torsion_filename = hit_id + '.torsion'
                        torsion_filename_hit = torsion_files_dir + hit_id + '.torsion'
			dssp_filename = hit_id + '.dssp'
			dssp_filename_hit = dssp_files_dir + hit_id + '.dssp'
			temp_dssp_dir = dssp_filename
			temp_aa_ss_sa_filename = hit_id + '.aa_ss_sa'
			exclamation_filename = hit_id + '.exc'

			#print pdb_filename_hit

                	#${get_dssp_files_dir}mkdssp -i $pdb_id_lower_core${chain_id}.pdb -o $pdb_id_lower_core${chain_id}.dssp
			#get_dssp_script_filename -i pdb_filename_hit -o dssp_filename_hit

                        chain_id_hit = hit_id[-1]
                        hit_id_raw = hit_id[0:4]
                        fasta_filename_hit = hit_id + '.fasta'

			if (not os.path.exists(pdb_filename_hit)):
				os.system('%s %s %s > %s' % (get_pdb_script_filename, hit_id_raw, chain_id_hit, fasta_filename_hit))
				os.system('\mv %s %s' % (pdb_filename, pdb_files_dir))
				os.system('rm -rf %s' % (fasta_filename_hit))

                        if (not os.path.exists(pdb_filename_hit)):
				print 'Could not obtain the pdb file for %s. Skipping this hit' % hit_id
				continue

                        if (not os.path.exists(dssp_filename_hit)):
                                os.system('%s -i %s -o %s' % (get_dssp_script_filename, pdb_filename_hit, dssp_filename_hit))

                        if (not os.path.exists(dssp_filename_hit)):
                                print 'Could not obtain the dssp file for %s. Skipping this hit' % hit_id
                                continue
		
			if (os.path.exists(exclamation_filename)):
	                        #os.system('\mv *.exc %s' % (dssp_files_dir))
                                os.system('rm -rf %s' % (exclamation_filename))

			if (os.path.exists(pdb_filename)):
				os.system('\mv %s %s' % (pdb_filename, pdb_files_dir))

			#compute torsion angle labels
			if (prediction_task == 3) or (prediction_task == 4):

				#if (not os.path.exists(torsion_filename_hit)):
				os.system("%s -p %s > %s" % (compute_torsion_angles_script_filename, pdb_filename_hit, torsion_filename_hit) )
				(phi_psi_omegas, aa_sequence_hit) = retrieve_phi_psi_omega(torsion_filename_hit)

				seq_length_pdb = len(phi_psi_omegas)
				seq_length_dssp = len(aa_sequence_hit)

				exclamation_present_flag = 0
				dssp_seq_in_pdb_seq_flag = 0

				if (seq_length_pdb != seq_length_dssp):
					exc_filename = dssp_files_dir + hit_id + '.exc'
					if (os.path.exists(exc_filename)):
						exc_file = open(exc_filename, 'r')
						aa_seq_2 = exc_file.readline()
						exc_file.close()
						exclamation_present_flag = 1
					else:
						#I will include the cases where the dssp based aa seq is a shifted version of the pdb based aa sequence
						if (aa_sequence_hit in aa_seq_2):
							dssp_seq_in_pdb_seq_flag = 1
						else:
							#output_file_2.write("%s\n" % protein_id_lowercase)
							continue
				else:
					aa_seq_2 = aa_sequence_hit

                        os.system('mkdir -p %s' % temp_dssp_dir)
                        os.system('cp %s %s' % (dssp_filename_hit, temp_dssp_dir))
                        os.system('python %s %s %s' % (extract_aa_ss_sa_labels_script_filename, temp_dssp_dir, temp_aa_ss_sa_filename))
                        os.system('rm -rf %s' % (temp_dssp_dir))
			os.system('rm -rf %s' % (exclamation_filename))
			if (os.path.exists(fasta_filename_hit)):
				os.system('rm -rf %s' % (fasta_filename_hit))

			temp_aa_ss_sa_file = open(temp_aa_ss_sa_filename, 'r')
			temp_aa_ss_sa_file_inside = temp_aa_ss_sa_file.readlines()
			temp_aa_ss_sa_file.close()

                        os.system('rm -rf %s' % (temp_aa_ss_sa_filename))

			#extract secondary structure and solvent accessibility labels
			fasta_id = temp_aa_ss_sa_file_inside[0].rstrip()
			protein_id = fasta_id[1:]
			chain_id = protein_id[4:5]
			protein_id_core = protein_id[:4]
			protein_id_lowercase = protein_id_core.lower() + chain_id
			aa_seq_hit_from_dssp = temp_aa_ss_sa_file_inside[1].rstrip()

                        ss_seq = temp_aa_ss_sa_file_inside[2].rstrip()
                        sa_seq = temp_aa_ss_sa_file_inside[3].rstrip()
                        if (prediction_task == 5):
                                sa_seq_2 = convert_sa_to_2_state(sa_seq, aa_seq_hit_from_dssp)

			if (prediction_task == 1) or (prediction_task == 3) or (prediction_task == 4):
	                        ss_seq_3_state = reduce_from_8_to_3_state_EHL_mapping(ss_seq)

			#print aa_seq_hit_from_dssp
			aa_sequence_hit = aa_seq_hit_from_dssp
	                        #if (len(ss_seq_3_state) != len(aa_sequence_hit)):
				#	continue

                        #if (prediction_task == 3) or (prediction_task == 4):
	                #        if (len(ss_seq_3_state) != len(aa_sequence_hit)) or (len(ss_seq_3_state) != len(phi_psi_omegas)):
        	        #                continue

			#if ('12e8H' in query_id):
			#	print torsion_filename_hit
			#	print hit_id
			#	print len(ss_seq_3_state)
			#	print len(phi_psi_omegas)
			#	print len(aa_sequence_hit)
			#	print len(aa_seq_2)

			if (prediction_task == 3):
				torsion_label_sequence_hit = convert_torsion_angles_to_5_state(phi_psi_omegas, aa_sequence_hit, ss_seq_3_state, aa_seq_2, exclamation_present_flag, dssp_seq_in_pdb_seq_flag)
				torsion_seq_length_hit = len(torsion_label_sequence_hit)

				if (torsion_seq_length_hit == 0):
					continue

			elif (prediction_task == 4):
                                torsion_label_sequence_hit = convert_torsion_angles_to_7_state(phi_psi_omegas, aa_sequence_hit, ss_seq_3_state, aa_seq_2, exclamation_present_flag, dssp_seq_in_pdb_seq_flag)
				torsion_seq_length_hit = len(torsion_label_sequence_hit)

				if (torsion_seq_length_hit == 0):
					continue

                        #aa_sequence_hit = extract_aa_sequence(torsion_filename_hit)
                        #torsion_label_sequence_hit = extract_torsion_label_sequence(torsion_filename_hit)

                        #starting row index of the j^th alignment
                        for k in range(alignments_start_index, n_lines-2, 1):
                                temp_line = hmm_alignment_file_inside[k].rstrip()
                                temp_line_2 = hmm_alignment_file_inside[k+1].rstrip()
                                temp_line_3 = hmm_alignment_file_inside[k+2].rstrip()
                                if ((string_to_check_1 in temp_line) and (string_to_check_2 in temp_line_2) and (string_to_check_3 in temp_line_3)):
                                        alignment_start_index = k
                                        break

                        file_corrupt_flag = 0
                        next_alignment_successfully_found_flag = 0
                        #ending row index of the j^th alignment
                        for k in range(alignment_start_index+2, n_lines-2, 1):
                                temp_line = hmm_alignment_file_inside[k].rstrip()
                                temp_line_2 = hmm_alignment_file_inside[k+1].rstrip()
                                temp_line_3 = hmm_alignment_file_inside[k+2].rstrip()

                                #if ('2h88D' in query_id):
                                #       print temp_line
                                #       print temp_line_2
                                #       print temp_line_3

                                if (string_to_check_6 in temp_line):
                                        alignment_end_index = k
                                        next_alignment_successfully_found_flag = 1
                                        break

                                if (string_to_check_6 in temp_line_2):
                                        alignment_end_index = k+1
                                        next_alignment_successfully_found_flag = 1
                                        break

                                if (string_to_check_6 in temp_line_3):
                                        alignment_end_index = k+2
                                        next_alignment_successfully_found_flag = 1
                                        break

                                if (string_to_check_4 in temp_line):
                                        token_list_temp = temp_line.split('No ')
                                        next_hit_id_number = int(token_list_temp[-1])
                                        #if ((j+1 > n_hits-1) or (next_hit_id_number != hit_indices[j+1])):
                                        #       file_corrupt_flag = 1
                                        #       break

                                if (string_to_check_5 in temp_line_2):
                                        token_list_temp = temp_line_2.split('>')
                                        next_hit_id = token_list_temp[1]
                                        #if ((j+1 > n_hits-1) or (next_hit_id not in hit_ids[j+1][0])):
                                        #       file_corrupt_flag = 1
                                        #       break
                                        #else:
                                        #       alignment_end_index = k
                                        #       next_alignment_successfully_found_flag = 1
                                        #       break
                                        alignment_end_index = k
                                        next_alignment_successfully_found_flag = 1
                                        break

                                        #if (string_to_check_5 in temp_line_2):
                                        #       token_list_temp = temp_line_2.split('>')
                                        #       next_hit_id = token_list_temp[1]
                                        #       if ((j+1 > n_hits-1) or (next_hit_id not in hit_ids[j+1][0])):
                                        #               file_corrupt_flag = 1
                                        #               break
                                        #       else:
                                        #               alignment_end_index = k
                                        #               next_alignment_successfully_found_flag = 1
                                        #               break

                        #print alignment_start_index
                        #print alignment_end_index

			#print file_corrupt_flag
			#print next_alignment_successfully_found_flag

                        if (file_corrupt_flag == 1):
                                break

                        if (next_alignment_successfully_found_flag == 0):
                                continue

                        #process the alignment
                        probability_line = hmm_alignment_file_inside[alignment_start_index+2].rstrip()
                        token_list_prob = probability_line.split()
                        probability_str_raw = token_list_prob[0]
                        token_list_prob2 = probability_str_raw.split('=')
                        probability_score = float(token_list_prob2[1])

                        if (probability_score < probability_score_threshold):
                                break

                        identity_str_raw = token_list_prob[4]
                        token_list_identity_str = identity_str_raw.split('=')
                        identity_str_raw_2 = token_list_identity_str[1]
                        perc_identity = float(identity_str_raw_2[0:-1]) / 100.0 #discards the % symbol at the end of the string

			#skip the hit if it is too similar to the query
			if (perc_identity > (percentage_identity_threshold / 100.0)):
				continue 

			if (perc_identity > max_seq_identity):
				max_seq_identity = perc_identity
	
			#print perc_identity
			#print probability_score
	
                        evalue_str_raw = token_list_prob[1]
                        token_list_evalue_str = evalue_str_raw.split('=')
                        e_value = float(token_list_evalue_str[1])

                        sim_score_str_raw = token_list_prob[2]
                        token_list_sim_score_str = sim_score_str_raw.split('=')
                        sim_score = float(token_list_sim_score_str[1])

                        n_aligned_cols_str_raw = token_list_prob[3]
                        token_list_n_aligned_cols_str = n_aligned_cols_str_raw.split('=')
                        n_aligned_cols = float(token_list_n_aligned_cols_str[1])

                        sim_score2_str_raw = token_list_prob[5]
                        token_list_sim_score2_str = sim_score2_str_raw.split('=')
                        sim_score2 = float(token_list_sim_score2_str[1])

                        sum_probs_str_raw = token_list_prob[6]
                        token_list_sum_probs_str = sum_probs_str_raw.split('=')
                        sum_probs = float(token_list_sum_probs_str[1])

                        #FROM ZHANG's MUSTER PAPER
                        e_value_weight = 1.0
                        if (e_value < 1e-10):
                                e_value_weight = 1.0
                        else:
                                e_value_weight = -0.05*math.log10(e_value) + 0.5

                        if (e_value_weight < 0):
                                e_value_weight = 0.0

                        read_query_flag = 0
                        read_hit_flag = 0
                        done = 0
                        max_n_lines_to_consider = alignment_end_index - alignment_start_index + 1
                        k = alignment_start_index

			#print "evalue weight: %f" % e_value_weight

                        average_confidence_score = 0.0
                        total_aligned_region_length = 0
                        torsion_sim_score = 0.0

                        while not done:

                                if (k >= alignment_end_index):
                                        done = 1
                                        continue

                                line = hmm_alignment_file_inside[k].rstrip()
                                #print line
                                #print query_id

                                if ((query_id_orig in line) and ('>' not in line)):
                                        token_list3 = line.split()
                                        starting_position_query = int(token_list3[2])
                                        ending_position_query = int(token_list3[4])
                                        aligned_sequence_query = token_list3[3]
                                        read_query_flag = 1

					#print hmm_alignment_file_inside[k+4]
					#print hmm_alignment_file_inside[k+5]
					#print hmm_alignment_file_inside[k+6]

                                        line2 = hmm_alignment_file_inside[k+4].rstrip() #this line should contain the hit
                                        line3 = hmm_alignment_file_inside[k+6].rstrip('\r\n') #this line should contain confidence

					if ('ss_pred' in line3):
						line3 = hmm_alignment_file_inside[k+7].rstrip('\r\n')

                                        if ('Confidence' not in line3):
                                                k += 1
                                                continue

                                        if (hit_id_orig not in line2):
                                                k += 1
                                                continue

                                        token_list4 = line2.split()
                                        token_list5 = line3.split()

                                        if (len(token_list4) != 6):
                                                k += 1
                                                continue

                                        starting_position_hit = int(token_list4[2])
                                        ending_position_hit = int(token_list4[4])
                                        aligned_sequence_hit = token_list4[3]

                                        confidence_sequence = line3[22:]

					#print confidence_sequence

                                        if (len(confidence_sequence) != len(aligned_sequence_hit)):
                                                print "*********************************************************"
                                                print confidence_sequence
                                                print aligned_sequence_hit
                                                print len(confidence_sequence)
                                                print len(aligned_sequence_hit)
                                                print hmm_alignment_filename

                                        k += 7

                                        aligned_region_length = len(aligned_sequence_query)

                                        position_indices_query = [0]*aligned_region_length
                                        position_indices_hit = [0]*aligned_region_length

                                        last_position_index_query = starting_position_query
                                        last_position_index_hit = starting_position_hit

                                        for m in range(aligned_region_length):
                                                if (aligned_sequence_query[m] != '-'):
                                                        position_indices_query[m] = last_position_index_query
                                                        last_position_index_query += 1
                                                if (aligned_sequence_hit[m] != '-'):
                                                        position_indices_hit[m] = last_position_index_hit
                                                        last_position_index_hit += 1

                                        total_aligned_region_length += aligned_region_length

					#must_be_corrected_flag = 0
					correction_factor = 0
					besli = ''
					dssp_seq_length = len(aa_sequence_hit)
					besli_matched_flag = 0

                                        for m in range(aligned_region_length-5):

						if (aligned_sequence_hit[m] != '-') and (aligned_sequence_hit[m+1] != '-') and (aligned_sequence_hit[m+2] != '-'):
							 position_index_hit = position_indices_hit[m]-1
                                                         position_index_hit_2 = position_indices_hit[m+1]-1
                                                         position_index_hit_3 = position_indices_hit[m+2]-1
                                                         position_index_hit_4 = position_indices_hit[m+3]-1
                                                         position_index_hit_5 = position_indices_hit[m+4]-1
							 
							 if (position_index_hit >= dssp_seq_length) or (position_index_hit_2 >= dssp_seq_length) or (position_index_hit_3 >= dssp_seq_length) or (position_index_hit_4 >= dssp_seq_length) or (position_index_hit_5 >= dssp_seq_length):
								continue

                                                         aa_residue_hit = aa_sequence_hit[position_index_hit]
                                                         aa_residue_hit_2 = aa_sequence_hit[position_index_hit_2]
                                                         aa_residue_hit_3 = aa_sequence_hit[position_index_hit_3]
                                                         aa_residue_hit_4 = aa_sequence_hit[position_index_hit_4]
                                                         aa_residue_hit_5 = aa_sequence_hit[position_index_hit_5]

							 besli += aligned_sequence_hit[m]
                                                         besli += aligned_sequence_hit[m+1]
                                                         besli += aligned_sequence_hit[m+2]
                                                         besli += aligned_sequence_hit[m+3]
                                                         besli += aligned_sequence_hit[m+4]

							 if (aa_residue_hit != aligned_sequence_hit[m]) or (aa_residue_hit_2 != aligned_sequence_hit[m+1]) or (aa_residue_hit_3 != aligned_sequence_hit[m+2]) or (aa_residue_hit_4 != aligned_sequence_hit[m+3]) or (aa_residue_hit_4 != aligned_sequence_hit[m+4]):
								correct_position_index = aa_sequence_hit.find(besli)
								correction_factor = position_index_hit - correct_position_index
								besli_matched_flag = 1
								break

					#print 'correction_factor = %d' % correction_factor

					if (besli_matched_flag == 0):
						continue

                                        for m in range(aligned_region_length):

                                                if ((aligned_sequence_query[m] != '-') and (aligned_sequence_hit[m] != '-')):

					
                                                        combined_posteriors = [0.0]*n_classes

                                                        position_index_query = position_indices_query[m]-1
                                                        position_index_hit = position_indices_hit[m]-1-correction_factor

                                                        if (position_index_query > len(aa_sequence_query)-1):
                                                                break

                                                        if (position_index_hit > len(aa_sequence_hit)-1):
                                                                break

                                                        aa_residue_query = aa_sequence_query[position_index_query]
                                                        aa_residue_hit = aa_sequence_hit[position_index_hit]

                                                        aa_number_query = give_me_aa_number(aa_residue_query)
                                                        aa_number_hit = give_me_aa_number(aa_residue_hit)

                                                        #residue_sim_score = blosum_matrix[aa_number_query][aa_number_hit]

							if (confidence_sequence[m] == ' '):
								confidence_score = 0
							else:
                                                        	confidence_score = int(confidence_sequence[m])

                                                        average_confidence_score += confidence_score

							if (prediction_task == 1):
								label_hit = ss_seq_3_state[position_index_hit]
								label_number_hit = give_me_ss_number_3_state(label_hit)
							elif (prediction_task == 2):
                                                                label_hit = ss_seq[position_index_hit]
                                                                label_number_hit = give_me_ss_number_8_state(label_hit)
                                                        elif (prediction_task == 3):
                                                        	label_hit = torsion_label_sequence_hit[position_index_hit]
	                                                        label_number_hit = give_me_torsion_number_ABEGO(label_hit)
                                                        elif (prediction_task == 4):
	                                                        label_hit = torsion_label_sequence_hit[position_index_hit]
                                                        	label_number_hit = give_me_torsion_number_LAMBEGO(label_hit)
                                                        elif (prediction_task == 5):
                                                                label_hit = sa_seq_2[position_index_hit]
                                                                label_number_hit = give_me_sa_number_2_state(label_hit)

                                                        pssm[position_index_query][label_number_hit] += e_value_weight*sim_score*pow(perc_identity, power_parameter)*confidence_score
							average_confidence_scores[position_index_query][label_number_hit] += confidence_score
							position_specific_hit_counts[position_index_query][label_number_hit] += 1
							#pssm[position_index_query][label_number_hit] += sim_score*pow(perc_identity, power_parameter)*aligned_region_length*sum_probs*confidence_score
                                                        #pssm[position_index_query][torsion_label_number_hit] += e_value_weight*sim_score*sim_score2*sum_probs*pow(perc_identity, power_parameter)*confidence_score
                                                        indices_covered[position_index_query] = 1 #those are the positions with match
							pssm_update_flag = 1
                                k += 1
                        n_hits_processed += 1
                        #print 'hit_id = %s' % hit_id

                        if (average_confidence_score != 0.0):
                                average_confidence_score /= total_aligned_region_length



		max_sequence_identities_file.write('%f\n' % max_seq_identity)

		print "PSSM update flag for %s: %d" % (query_id_orig, pssm_update_flag)

                #normalize pssm
                for j in range(n_aas_query):
                        row_sum = 0.0
                        for k in range(n_classes):
                                row_sum += pssm[j][k]
                        for k in range(n_classes):
                                pssm[j][k] /= row_sum
				if (position_specific_hit_counts[j][k] != 0):
					average_confidence_scores[j][k] /= float(position_specific_hit_counts[j][k])

		token_list = fasta_filename.split('/')
		fasta_filename_wo_path = token_list[-1]
		protein_id = fasta_filename_wo_path.replace('.fasta', '')

		#write the pssm to the output file
		#if (prediction_task == 1):
	        #        output_filename = output_files_dir + protein_id + '_ss3_sp_1.struct'

		#elif (prediction_task == 2):
	        #        output_filename = output_files_dir + protein_id + '_ss8_sp_1.struct'

		#elif (prediction_task == 3):
                #        output_filename = output_files_dir + protein_id + '_ta5_sp_1.struct'

		#elif (prediction_task == 4):
                #        output_filename = output_files_dir + protein_id + '_ta7_sp_1.struct'

		#elif (prediction_task == 5):
                #        output_filename = output_files_dir + protein_id + '_sa2_sp_1.struct'

		#if ('3ecsA' not in query_id_orig):
		#	continue
		#print output_filename

		output_filename = sp_1_filename

		output_file = open(output_filename, 'w')
                average_confidence_scores_filename = output_files_dir + query_id_orig + '.conf'
                average_confidence_scores_file = open(average_confidence_scores_filename, 'w')

		for j in range(n_aas_query):
			for k in range(n_classes):
				output_file.write('%.11f ' % pssm[j][k])
	                        average_confidence_scores_file.write('%.3f ' % average_confidence_scores[j][k])

			output_file.write('\n')
                        average_confidence_scores_file.write('\n')

		output_file.close()
		average_confidence_scores_file.close()

	max_sequence_identities_file.close()

