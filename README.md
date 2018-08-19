# DSPRED
DSPRED is a two-stage hybrid classifier for predicting one-dimensional structure of proteins such as secondary structure, solvent accessibility and torsion angle classes. It employs a committee of dynamic Bayesian networks at the first stage and a support vector machine at the second stage. The feature sets of DSPRED include position-specific scoring matrices (PSSMs) and structural profile matrices derived by PSI-BLAST and HHblits methods. 

# Operating System Requirements
DSPRED runs in Linux/Unix operating system.

# Software Requirements
The following software must be installed for DSPRED along with the protein databases for PSI-BLAST and HHblits:
1. PSI-BLAST from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
2. HHblits: https://github.com/soedinglab/hh-suite 
3. Rosetta (get_pdb.py and amino_acids.py scripts are required only to download PDB files): https://www.rosettacommons.org/software 
4. Graphical Models Toolkit (GMTK): https://melodi.ee.washington.edu/gmtk/
5. libSVM: https://www.csie.ntu.edu.tw/~cjlin/libsvm/

For PSI-BLAST, the NR database and for HHblits, the Uniprot as well as PDB70 databases should be downloaded. Detailed installation instructions for HHblits can be found [here](install_hhsuite). HMMER3 software is used to compute the second structural profile matrix and its binaries are already under code/scripts/pfamscan/hmmer3/. If that does not work you can install HMMER by visiting http://www.hmmer.org and setting the path variable called hmmer_dir inside code/scripts/dspred file.

# Downloading the DSPRED files

1. cd to the directory where you want to install DSPRED
2. git clone git clone https://github.com/yusufzaferaydin/dspred.git
3. cd dspred/scripts

# Setting the paths and running DSPRED

1. Set the paths in dspred file 
2. ./dspred -h will print the help menu and parameter options of DSPRED
3. ./run_dspred includes example command line for running dspred

# Feedback and comments

You are welcome to send comments to zafer.aydin@agu.edu.tr

