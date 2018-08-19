__author__ = "Zafer Aydin"
__project__ = "Protein Structure Prediction"
__date__ = "23 April 2018"

usage = """
	Usage: python extract_identities.py <identities_for_protein_filename> <max_identity_filename>
"""

import sys

if __name__ == "__main__":

	if (len(sys.argv) != 5):
		sys.stderr.write(usage)
		sys.exit(1)

	identities_for_protein_filename = sys.argv[1]
	pc_thr = float(sys.argv[2])
	max_identity_filename = sys.argv[3]
	skip_top_hit_flag = int(sys.argv[4])

	identities_for_protein_file = open(identities_for_protein_filename, 'r')
	identities_for_protein_file_inside = identities_for_protein_file.readlines()
	identities_for_protein_file.close()
	n_lines = len(identities_for_protein_file_inside)

	max_identity_file = open(max_identity_filename, 'w')

	identities = []

	for i in range(n_lines):

		line = identities_for_protein_file_inside[i].rstrip()

                token_list = line.split()

		aligned_cols_str = token_list[3]
		token_list_2 = aligned_cols_str.split('=')
		aligned_cols_value = int(token_list_2[1])

		if (aligned_cols_value < 10):
			continue

		identity_str = token_list[4]
		token_list_2 = identity_str.split('=')
		identity_str_2 = token_list_2[1]
		identity_value = float(identity_str_2[:-1])

		if (identity_value == 0):
			continue

		if (identity_value > pc_thr):
			continue

                if ((i == 0) and (skip_top_hit_flag == 1) and (identity_value == 100)):
                        continue

		identities += [ identity_value ]

	if identities:
		max_identity = max(identities)
	else:
		max_identity = 0.0

	max_identity_file.write("%.3f\n" % max_identity)
	max_identity_file.close()

