__author__ = "Zafer Aydin"
__date__ = "12/05/2008"

#FILE: generate_postp_files.py
#PROJECT: Protein Secondary Structure Prediction

usage = """

USAGE: python generate_postp_files.py <number_of_states> <test_postp_files_filename> <jt_output_filename>

"""
import sys
import os
from  struct import *
from  math   import *

if __name__ == "__main__":


	if ( len(sys.argv) != 4):
		
		sys.stderr.write(usage)
		sys.exit(1)
	
	number_of_states = int(sys.argv[1])
        jt_output_filename = sys.argv[2]
	test_postp_files_filename = sys.argv[3]

	jt_output_file = open( jt_output_filename, 'r')
	test_postp_files_file = open(test_postp_files_filename, 'r')                                
	done = 0
	number_of_seqs_read = 0

	while not done:

		a_line = jt_output_file.readline()
		
		if (len(a_line) == 0):

			done = 1
			continue
		
		if ( not a_line.startswith("Segment") ):

			continue

		test_postp_filename = test_postp_files_file.readline()
		test_postp_filename = test_postp_filename.rstrip()

		test_postp_file = open( test_postp_filename, 'w')

		done2 = 0

		while not done2:

			a_line = jt_output_file.readline()

			if (len(a_line) == 0):
				
				done = 1 
	                        done2 = 1
        	                continue 

			if ( (a_line.startswith("Number") or a_line.startswith("Segment") or a_line.startswith("###")) ):
		
				done2 = 1
				continue

			if ( not a_line.startswith("Partition") ):                        
                        	continue 
					
			done3 = 0
			
			prob_values = [0] * number_of_states

			while not done3:

				a_line = jt_output_file.readline()

				if (len(a_line) == 0):

                                	done = 1
                                	done2 = 1
					done3 = 1
                                	continue

				if ( (a_line.startswith("Number") or a_line.startswith("Segment") or a_line.startswith("###"))):

					done2 = 1
					done3 = 1
					continue

				if ( a_line.startswith("-----") ):
				
					done3 = 1
					continue

				token_list = a_line.split()
        			
				if ( len(token_list) != 3 ):
            				
					print ' ERROR in getCliqueProbs ??? more than 3 tokens ??? ', len(token_list)
            				print a_line, token_list
            				sys.exit(-1)

				postp_value = float ( token_list[1] )
				t_list = token_list[2].split('=')

				state_class = int( t_list[-1] )

				prob_values[state_class] = postp_value

			for i in range(0, number_of_states, 1):

				test_postp_file.write('%.5f\t' % prob_values[i])
			
			test_postp_file.write('\n')
		
		test_postp_file.close()
		number_of_seqs_read += 1
	
	if ( number_of_seqs_read == 0):

		print 'Error: No Posterior Probability is Read !'

	test_postp_files_file.close()
	jt_output_file.close()

