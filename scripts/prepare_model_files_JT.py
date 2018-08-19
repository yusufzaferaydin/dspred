__author__ = "Zafer Aydin"
__date__ = "11/24/2009"
__project__ = "Protein Structure Prediction"

usage="""

        Usage: prepare_model_files_JT.py <n_states> <lss> <Dmax> <new_indices_filename> <model filename JT>
"""

import sys
import random
import math

if __name__ == "__main__":

        if ( len(sys.argv) != 12 ):

                sys.stderr.write(usage)
                sys.exit(1)

	n_states = int(sys.argv[1])
	lss = int(sys.argv[2])
	Dmax = int(sys.argv[3])
        new_indices_filename = sys.argv[4]
	models_path = sys.argv[5]
	fixed_params_filename = sys.argv[6]
        learned_means_filename = sys.argv[7]
        learned_diag_covs_filename = sys.argv[8]
        learned_dlinks_filename = sys.argv[9]
        learned_dense_cpts_filename = sys.argv[10]
	model_filename_JT = sys.argv[11]

        new_indices_file = open(new_indices_filename, 'r')
        new_indices_file_inside = new_indices_file.readlines()
        new_indices_file.close()

	ngc = pow(n_states, (lss+1))
 
	model_file_JT = open(model_filename_JT, 'w')

	model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
	model_file_JT.write('%\n')
	model_file_JT.write('% This is the master file for the model_train graph\n')
	model_file_JT.write('%\n')
	model_file_JT.write('\n')
	model_file_JT.write('#include \"common_params.txt\"\n')
	model_file_JT.write('\n')
	model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
	model_file_JT.write('%                              DIRICHLET TABLES                               %\n')
	model_file_JT.write('% none needed at this time ...\n')
	model_file_JT.write('\n')
	model_file_JT.write('DIRICHLET_TAB_IN_FILE   inline\n')
	model_file_JT.write('0\n')	
	model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                               DECISION TREEs                                %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% a decision tree is needed whenever we use a DeterministicCPT or a mapping;  %\n')
        model_file_JT.write('% here we list first those DTs that are paired with a DeterministicCPT of the %\n')
        model_file_JT.write('% same name (w/o the _DT ending), then those that are used as part of a       %\n')
        model_file_JT.write('% mapping of some sort ...\n')
        model_file_JT.write('\n')
        model_file_JT.write('DT_IN_FILE      inline\n')
        model_file_JT.write('4\n')
	model_file_JT.write('\n')
        model_file_JT.write('0                               % index\n')
        model_file_JT.write('start_history_DT                % name\n')
        model_file_JT.write('1                               % one parent (state_class)\n')

	
        model_file_JT.write('    0    %d    ' % n_states)
	for i in range(n_states-1):
		model_file_JT.write('%d ' % i)
	model_file_JT.write('default\n')

	for i in range(n_states):

		k = ngc-n_states+i
		model_file_JT.write('        -1   %d\n' % k)

	model_file_JT.write('1                               % index\n')
	model_file_JT.write('update_history_DT                % name\n')
        model_file_JT.write('2                               % two parents (current state_class, state class history (this also includes current state)\n')
	model_file_JT.write('    0    %d  ' % ngc)
	for i in range(ngc-1):
                model_file_JT.write('%d ' % i)
        model_file_JT.write('default\n')
	
	temp_value = ngc / n_states

	for i in range(n_states):

		for j in range(temp_value):

			k = j*n_states

			if (k != 0):
				model_file_JT.write('        -1   { p1 + %d }\n' % k)

			else:
				model_file_JT.write('        -1   { p1 }\n')

        model_file_JT.write('\n')
        model_file_JT.write('2                               % index\n')
        model_file_JT.write('1D_map                          % name\n')
        model_file_JT.write('1                               % one parent\n')
        model_file_JT.write('        -1 {p0}                 % returns the parent\n')
        model_file_JT.write('\n')
        model_file_JT.write('3                               % index\n')
        model_file_JT.write('tuple_map                          % name\n')
        model_file_JT.write('1                               % one parent\n')
        model_file_JT.write('    0    %d  ' % ngc)
        for i in range(ngc-1):
                model_file_JT.write('%d ' % i)
        model_file_JT.write('default\n')
        for i in range(ngc):
                new_tuple_index = int(new_indices_file_inside[i].rstrip())
                model_file_JT.write('        -1 %d        \n' % new_tuple_index)
        model_file_JT.write('\n')

        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                             DETERMINISTIC CPTs                              %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('\n')
        model_file_JT.write('DETERMINISTIC_CPT_IN_FILE       inline\n')
        model_file_JT.write('\n')
        model_file_JT.write('2                               % total number of Det CPTs to follow\n')
        model_file_JT.write('\n')
        model_file_JT.write('0                               % index\n')
        model_file_JT.write('start_history_CPT               % name\n')
        model_file_JT.write('1                               % number of parents\n')
        model_file_JT.write('STATE_CARD STATE_CLASS_HISTORY_CARD  % parent and self cardinalities\n')
        model_file_JT.write('start_history_DT                % DT name\n')
        model_file_JT.write('\n')
        model_file_JT.write('1                               % index\n')
        model_file_JT.write('update_history_CPT              % name\n')
        model_file_JT.write('2                               % number of parents\n')
        model_file_JT.write('STATE_CLASS_HISTORY_CARD STATE_CARD STATE_CLASS_HISTORY_CARD  % parent and self cardinalitie\n')
        model_file_JT.write('update_history_DT                % DT name\n')
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                                 DENSE CPTs                                  %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('\n')
        model_file_JT.write('DENSE_CPT_IN_FILE       %s                ascii\n' % (models_path+learned_dense_cpts_filename))
	model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                               DLINK STRUCTURES                              %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% Dlink structures specify the directed dependencies that may exist between   %\n')
        model_file_JT.write('% the individual elements of observation vectors.  Given a particular element %\n')
        model_file_JT.write('% of an observation vector, a dependency might exist with a parent in the     %\n')
        model_file_JT.write('% same frame as the current element or with a parent in some previous or      %\n')
        model_file_JT.write('% future frame.                                                               %\n')
        model_file_JT.write('% ( see section 4.5 in GMTK reference documentation )                         %\n')
        model_file_JT.write('\n')
        model_file_JT.write('DLINK_IN_FILE                  %s           ascii\n' % (models_path+fixed_params_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                  MEANs and COVARs for GAUSSIAN COMPONENTS                   %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('\n')
        model_file_JT.write('MEAN_IN_FILE                   %s           ascii\n' % (models_path+learned_means_filename))
        model_file_JT.write('COVAR_IN_FILE                  %s           ascii\n' % (models_path+learned_diag_covs_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                               DLINK MATRICES                                %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% A dlink matrix has a coefficient for each dependency specified in a dlink   %\n')
        model_file_JT.write('% structure.\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% We will have as many dlink matrices as we have emission classes ...\n')
        model_file_JT.write('\n')
        model_file_JT.write('DLINK_MAT_IN_FILE               %s          ascii\n' % (models_path+learned_dlinks_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                               GAUSSIAN COMPONENTS                           %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% each mixture component is specified on a single line with the following     %\n')
        model_file_JT.write('% pieces of information :                                                     %\n')
        model_file_JT.write('%       a) mixture component index\n')
        model_file_JT.write('%       b) mixture component dimensionality\n')
        model_file_JT.write('%       c) type of component (1 means uses a dlink, 0 means just mean, covar)\n')
        model_file_JT.write('%       d) name of mixture component\n')
        model_file_JT.write('%       e) name of mean vector\n')
        model_file_JT.write('%       f) name of covariance vector\n')
        model_file_JT.write('%       g) name of dlink structure (if type=1)\n')
        model_file_JT.write('\n')
        model_file_JT.write('MC_IN_FILE              %s      		ascii\n' % (models_path+fixed_params_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                    DENSE PROBABILITY MASS FUNCTIONS                         %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% these are needed as weight vectors to combine one or more components into a \n')
        model_file_JT.write('% mixture of components -- we are only using single-component mixtures, so    %\n')
        model_file_JT.write('% we will only need a single weight vector with one weight=1                  %\n')
        model_file_JT.write('\n')
        model_file_JT.write('DPMF_IN_FILE            %s      ascii\n' % (models_path+fixed_params_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                          MIXTURES OF GAUSSIAN COMPONENTS                    %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% Here we define the \'mixtures\' of componenets (which in this case are not    %\n')
        model_file_JT.write('% really mixtures at all since they contain only one componenet) ...          %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% Each mixture is defined on a single line with the following pieces of       %\n')
        model_file_JT.write('% information :\n')
        model_file_JT.write('%       a) mixture index\n')
        model_file_JT.write('%       b) mixture dimensionality\n')
        model_file_JT.write('%       c) name of mixture\n')
        model_file_JT.write('%       d) number of components\n')
        model_file_JT.write('%       e) name of vector of mixing weights\n')
        model_file_JT.write('%       f) name of 1st component\n')
        model_file_JT.write('%       g) name of 2nd component, etc ...\n')
        model_file_JT.write('\n')
        model_file_JT.write('MX_IN_FILE              %s   		ascii\n' % (models_path+fixed_params_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                              NAME COLLECTIONs                               %\n')
        model_file_JT.write('%                                                                             %\n')
        model_file_JT.write('% and finally we define a collection of gaussian mixtures ...                 %\n')
        model_file_JT.write('NAME_COLLECTION_IN_FILE  %s   ascii\n' % (models_path+fixed_params_filename))
        model_file_JT.write('\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('%                                    END                                      %\n')
        model_file_JT.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        model_file_JT.write('\n')
	

	model_file_JT.close()

	
