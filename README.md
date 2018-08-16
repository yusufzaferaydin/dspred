# DSPRED
DSPRED is a two-stage hybrid classifier for predicting one-dimensional structure of proteins such as secondary structure, solvent accessibility and torsion angle classes. It employs a committee of dynamic Bayesian networks at the first stage and a support vector machine at the second stage. The feature sets of DSPRED include position-specific scoring matrices (PSSMs) and structural profile matrices derived by PSI-BLAST and HHblits methods. 

# Operating System Requirements
DSPRED runs in Linux/Unix operating system.

# Software Requirements
The following software must be installed for DSPRED along with the protein databases for PSI-BLAST and HHblits:
1. PSI-BLAST from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
2. HHblits: https://github.com/soedinglab/hh-suite 
3. HMMER3 (required for DSPRED v2 only): http://www.hmmer.org
4. Graphical Models Toolkit (GMTK): https://melodi.ee.washington.edu/gmtk/
5. libSVM: https://www.csie.ntu.edu.tw/~cjlin/libsvm/

For PSI-BLAST, the NR database and for HHblits, the Uniprot as well as PDB70 databases should be downloaded. Detailed installation instructions for HHblits can be found [here](install_hhblits).

# Setting Paths

