__author__ = "Zafer Aydin"
__date__ = "04/29/2010"
__project__ = "Protein Structure Prediction"

import sys
import os
#import Numeric

#usage = """USAGE: dssp2fasta_detailed.py [options] <DSSP directory> <sequence fasta> <label fasta>
usage = """USAGE: dssp2fasta_detailed.py [options] <DSSP directory> <output dataset filename>

  Reads each file in <DSSP directory>.  Adds the amino acid sequence, secondary structure sequence, solvent accessibility sequence and base pairing sequences in fasta format to the dataset file

  The secondary structure labels are as follows:

     H = alpha helix
     e = extended strand, participates in one beta bridge
     E = extended strand, participates in two beta bridges
     b = bulge (i.e., strand with zero beta bridges)
     B = residue in isolated beta-bridge
     G = 3-helix (3/10 helix)
     I = 5 helix (pi helix)
     T = hydrogen bonded turn
     S = bend
     L = loop

  The script outputs solvent accessibility labels.
  These are computed on an absolute DSSP scale (0-200), or relative to
  the per-protein maximum (0-1).  The alphabet goes alphabetically
  from A, using thresholds defined by the user.

  Options:
    --accessibility absolute|relative
    --thresholds <int>[,<int>]+

"""

#############################################################################
def print_fasta (id, chain_ids, sequence, ss_seq, acc_seq, bp1_seq, bp2_seq, position_seq, outfile):

  # Refuse to print the same chain ID twice.
  printed_ids = {}
  
  if ("!" in sequence):
    chain_id_list = chain_ids.split("!")

    ss_seq_list = ss_seq.split("!")
    #acc_seq_list = acc_seq.split("!")
    #bp1_seq_list = bp1_seq.split("!")
    #bp2_seq_list = bp2_seq.split("!")
    
    seq_list = sequence.split("!")

    n_seqs = len(seq_list)
    starts_ends = [0]*n_seqs
    for i in range(n_seqs):
       starts_ends[i] = [0]*2

    seq_start_index = 0
    seq_end_index = -2

    for i in range(n_seqs):
       temp_seq = seq_list[i]
       seq_end_index += len(temp_seq) + 1
       starts_ends[i][0] = seq_start_index
       starts_ends[i][1] = seq_end_index
       seq_start_index = seq_end_index + 2

    seq_counter = 0    
    for k in range(n_seqs):

      seq = seq_list[k]
      seq_start_index = starts_ends[k][0]
      seq_end_index = starts_ends[k][1]
 
      ss_seq_selected = ss_seq_list[k]
      acc_seq_selected = acc_seq[seq_start_index:(seq_end_index+1)]
      bp1_seq_selected = bp1_seq[seq_start_index:(seq_end_index+1)]
      bp2_seq_selected = bp2_seq[seq_start_index:(seq_end_index+1)]
      position_seq_selected = position_seq[seq_start_index:(seq_end_index+1)]

      chain_id = chain_id_list[0][0]
      if (printed_ids.has_key(chain_id)):
        sys.stderr.write("Skipping repeated ID (%s_%s).\n" % (id, chain_id))
      else:
        outfile.write(">%s_%s\n%s\n%s\n" % (id, chain_id, seq, ss_seq_selected))

        n_res = len(seq)
        n_accs = len(acc_seq_selected)
        n_bp1s = len(bp1_seq_selected)
        n_bp2s = len(bp2_seq_selected)
	n_position_seqs = len(position_seq_selected)

        if ((n_accs != n_res) or (n_bp1s != n_accs) or (n_position_seqs != n_res)):
            print n_res
            print n_accs
            print n_bp1s
            print seq
            print acc_seq_selected
            print bp1_seq_selected
	    print n_position_seqs
            sys.stderr.write("The lengths do not match for protein (%s_%s).\n" % (id, chain_id))
            
        for i in range(n_accs):
           outfile.write("%d " % int(acc_seq_selected[i]))
        outfile.write("\n")
        #for i in range(n_bp1s):
        #   outfile.write("%d " % int(bp1_seq_selected[i]))
        #outfile.write("\n")
        #for i in range(n_bp2s):
        #   outfile.write("%d " % int(bp2_seq_selected[i]))
        #outfile.write("\n")
	#for i in range(n_position_seqs):
	#   outfile.write("%d " % int(position_seq_selected[i]))
	#outfile.write("\n")

        printed_ids[chain_id] = 1

      chain_id_list = chain_id_list[1:]
      seq_counter += 1

  else:
    outfile.write(">%s\n%s\n%s\n" % (id, sequence, ss_seq))
    n_res = len(sequence)
    n_accs = len(acc_seq)
    n_bp1s = len(bp1_seq)
    n_bp2s = len(bp2_seq)
    n_position_seqs = len(position_seq)

    if ((n_accs != n_res) or (n_bp1s != n_accs)):
       sys.stderr.write("The lengths do not match for protein (%s_%s).\n" % (id, chain_id))
    for i in range(n_accs):
       outfile.write("%d " % int(acc_seq[i]))
    outfile.write("\n")
    #for i in range(n_bp1s):
    #   outfile.write("%d " % int(bp1_seq[i]))
    #outfile.write("\n")
    #for i in range(n_bp2s):
    #   outfile.write("%d " % int(bp2_seq[i]))
    #outfile.write("\n")
    #for i in range(n_position_seqs):
    #   outfile.write("%d " % int(position_seq[i]))
    #outfile.write("\n")

#############################################################################
def retrieve_dssp_sequences(accessibility_type, dssp_directory, dssp_filename):

  # Initialize the return variables.
  chain_ids = ""
  sequence = ""
  sequence2 = ""  
  labels = ""
  acc_seq = []
  bp1_seq = []
  bp2_seq = []
  position_seq = []

  # Open the file for reading.
  dssp_file = open(os.path.join(dssp_directory, dssp_filename), "r")

  exclamation_flag = 0

  # Read the header.
  for line in dssp_file:
    if (line[0:3] == "  #"):
      break

  # Read the rest of the file.
  prev_label = ""
  for line in dssp_file:
    chain_id = line[11]
    amino = line[13]
    label = line[16]
    accessibility = line[35:38]
    bp1 = line[25:29]
    bp2 = line[29:33]
    position_index = line[0:5]
 
    if (bp1 == "****"):
        bp1 = "10000" #Since we don't know the actual value of the residue index that is greater than 10,000
    if (bp2 == "****"):
        bp2 = "10000" #Since we don't know the actual value of the residue index that is greater than 10,000

    #print line
    #print amino
    #print label
    #print accessibility
    #print bp1
    #print bp2

    # Chain breaks are indicated by "!*".  These are different from
    # "!" alone, which only indicates a discontinuity of backbone
    # coordinates.
    if (amino == "!"):
      exclamation_flag = 1
      #print 'Exclamation symbol is present for %s' % dssp_filename

      if (line[14] == "*"):
        chain_ids = chain_ids + "!"
        sequence = sequence + "!"
        labels = labels + "!"
        prev_label = "!"
	acc_seq += [ "!" ] 
        bp1_seq += [ "!" ]
        bp2_seq += [ "!" ]
	position_seq += [ "!" ]
      sequence2 = sequence2 + "!"
      continue

    # Disulfide bridges are indicated by lowercase.
    if (amino.islower()):
      amino = "C"

    # Loop or irregular is indicated by a space.
    if (label == " "):
      label = "L"

    # Make sure the label is valid.
    if (label not in "HBEGITSL"):
      sys.stderr.write("ERROR: Invalid label (%s) in %s\n%s" %
                       (label, dssp_filename, line))
      sys.exit(1)

    # Make sure the amino acid is valid.
    if (amino not in "ACDEFGHIKLMNOPQRSTVWYXBZ"):
      sys.stderr.write("ERROR: Invalid amino acid (%s) in %s\n%s" %
                       (amino, dssp_filename, line))
      sys.exit(1)

    # Discriminate between strands with one or two beta bridges.
    if (label == "E"):
#      sys.stderr.write("<%s> <%s>\n" % (line[25:29], line[29:33]))

      # Count the number of bridges.
      num_bridges = 0
      for bridge in (line[25:29], line[29:33]):
        if ((bridge == "****") or # Handle index >10,000.
            (int(bridge) != 0)):
          num_bridges = num_bridges + 1

      # If we got no bridges, then we're in a bulge.
      if (num_bridges == 0):
        label = "b"

      # Otherwise, indicate a single bridge with lowercase.
      elif (num_bridges == 1):
        label = "e"
      
    # Extend the two sequences.
    chain_ids = chain_ids + chain_id
    sequence = sequence + amino
    sequence2 = sequence2 + amino
    labels = labels + label
    acc_seq += [ str(int(accessibility)) ]
    bp1_seq += [ str(int(bp1)) ]
    bp2_seq += [ str(int(bp2)) ] 
    position_seq += [ str(int(position_index)) ]

    #if (accessibility_type == "none"):
    #  labels = labels + label
    #else:
    #  labels = labels + " " + accessibility + " "

    prev_label = label

  # Close the file and return.
  dssp_file.close()

  if (exclamation_flag == 1):
	  token_list = dssp_filename.split('.dssp')
	  protein_id = token_list[0]
	  exclamations_filename = protein_id + '.exc'
	  exclamations_file = open(exclamations_filename, 'w')
	  exclamations_file.write("%s\n" % sequence2)
	  exclamations_file.close()
  return(chain_ids, sequence, labels, acc_seq, bp1_seq, bp2_seq, position_seq)
                
#############################################################################
ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
def accessibility_alphabet(accessibility, threshold_string, labels):


  # Get the maximum accessibility score for this protein.
  accessibility_scores = labels.split()
  max_score = 0.0
  for score in accessibility_scores:
    if (score != "!") and (float(score) > max_score):
      max_score = int(score)

  # Convert the given thresholds into floats.
  threshold_strings = threshold_string.split(",")
  if (len(threshold_strings) > len(ALPHABET)):
    sys.stderr.write("Too many solvent accessibility thresholds (%d).\n"
                     % len(threshold_strings))
    sys.exit(1)
  thresholds = []
  for string in threshold_strings:
    thresholds.append(float(string))
  #thresholds = Numeric.sort(thresholds)
  thresholds = thresholds.sort()
  alphabet = ALPHABET[0:len(threshold_strings) + 1]

  # Convert each observed value into a letter.
  return_value = ""
  for score in accessibility_scores:
    letter = ""
    if (score == "!"):
      letter = "!"
    else:

      # Use the relative or absolute score.
      converted_score = float(score)
      if (accessibility == "relative"):
        converted_score = converted_score / max_score

      # Find the corresponding score threshold.
      for i_threshold in range(0, len(thresholds)):
        if (converted_score < thresholds[i_threshold]):
          letter = alphabet[i_threshold]
      if (letter == ""):
        letter = alphabet[-1]

    return_value = return_value + letter

  return(return_value)

#############################################################################
# MAIN
#############################################################################

# Set defaults.
#accessibility = "none"
accessibility = "absolute"
thresholds = ""

# Parse the command line.
sys.argv = sys.argv[1:]
while (len(sys.argv) > 3):
  next_arg = sys.argv[0]
  sys.argv = sys.argv[1:]
  if (next_arg == "--accessibility"):
    accessibility = sys.argv[0]
    sys.argv = sys.argv[1:]
  elif (next_arg == "--thresholds"):
    thresholds = sys.argv[0]
    sys.argv = sys.argv[1:]
  else:
    sys.stderr.write("Invalid option (%s).\n" % next_arg)
    sys.exit(1)

if (len(sys.argv) != 2):
  sys.stderr.write(usage)
  sys.exit(1)
dssp_directory = sys.argv[0]
dataset_filename = sys.argv[1]
#label_filename = sys.argv[2]

# If thresholds weren't set, use defaults.
if ((accessibility != "none") and (thresholds == "")):
  if (accessibility == "absolute"):
    #thresholds = "15"
    thresholds = "25"
    #sys.stderr.write("Using default threshold %s.\n" % thresholds)
  elif (accessibility == "relative"):
    #thresholds = "0.15"
    thresholds = "0.25"
    #sys.stderr.write("Using default threshold %s.\n" % thresholds)
  else:
    sys.stderr.write("Invalid accessibility option (%s).\n" % accessibility)
    sys.exit(1)

# Open the two output files.
dataset_file = open(dataset_filename, "w")
#label_file = open(label_filename, "w")

# Traverse the directory.
for dssp_filename in (os.listdir(dssp_directory)):
  #sys.stderr.write("%s\n" % dssp_filename)
  
  # Check that the filename is OK.
  if (dssp_filename[-5:] != ".dssp"):
    sys.stderr.write("Skipping file %s.\n" % dssp_filename)
    continue

  # Parse the DSSP file.
  (chain_ids, sequence, ss_seq, acc_seq, bp1_seq, bp2_seq, position_seq) = retrieve_dssp_sequences(accessibility,
                                                          dssp_directory,
                                                          dssp_filename)

  # If necessary, convert accessibility to an alphabet.
  #if (accessibility != "none"):
  #  labels = accessibility_alphabet(accessibility, thresholds, labels)

  # Write out the sequences.
  protein_id = dssp_filename[:-5]
  print_fasta(protein_id, chain_ids, sequence, ss_seq, acc_seq, bp1_seq, bp2_seq, position_seq, dataset_file)
  #print_fasta(protein_id, chain_ids, labels, label_file)

# Close the output files.
dataset_file.close()
#label_file.close()

# Local variables:
# mode: python
# py-indent-offset: 2
# End:
