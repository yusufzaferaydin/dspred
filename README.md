# DSPRED
DSPRED is a two-stage hybrid classifier for predicting one-dimensional structure of proteins such as secondary structure, solvent accessibility and torsion angle classes. It employs a committee of dynamic Bayesian networks at the first stage and a support vector machine at the second stage. The feature sets of DSPRED include position-specific scoring matrices (PSSMs) and structural profile matrices derived by PSI-BLAST and HHblits methods. 

# Operating System Requirements
DSPRED runs in Linux/Unix operating system.

# Software Requirements
The following software must be installed for DSPRED along with the protein databases for PSI-BLAST and HHblits:
1. PSI-BLAST from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
2. HHblits: https://github.com/soedinglab/hh-suite 
3. Rosetta (get_pdb.py and amino_acids.py scripts are required only to download PDB files): https://www.rosettacommons.org/software 
4. Install Moose for Perl using the following command (needs admin previleges): sudo cpan Moose
5. Graphical Models Toolkit (GMTK): https://melodi.ee.washington.edu/gmtk/
6. libSVM: https://www.csie.ntu.edu.tw/~cjlin/libsvm/

For PSI-BLAST, the NR database and for HHblits, the Uniprot as well as PDB70 databases should be downloaded. Detailed installation instructions for HHblits can be found [here](install_hhsuite). 

## Note about PfamScan and SP2

The second structural profile matrix (SP2) is computed using the PfamScan program (https://www.ebi.ac.uk/seqdb/confluence/display/THD/PfamScan), which is available under scripts/pfamscan/ folder.  It uses HMMER3 software to compute sequence-HMM-profile alignments, as well as Perl, Moose for Perl, BioPerl, blastp binary of BLAST, T-Coffee and hhmake binary of HHblits. HMMER's binaries are available by default under scripts/pfamscan/hmmer3/, BioPerl under scripts/pfamscan/bioperl-1.4/ and T-Coffee under scripts/pfamscan/tcoffee/. We already provide all the binaries for PfamScan for convenience. However, if a problem occurs during SP2 computation you can install PfamScan and its dependencies by reading [here](scripts/pfamscan/README_pfamscan). PfamScan program depends on the following software:

1. HMMER3: http://www.hmmer.org
2. Perl and Moose for Perl
3. BioPerl (must be installed under PfamScan dir): https://metacpan.org/pod/BioPerl 
4. blastp binary of NCBI's BLAST
5. hhmake binary of HHblits
6. T-Coffee: http://tcoffee.crg.cat/
7. get_pdb.py and amino_acids.py of Rosetta: https://www.rosettacommons.org/software 

If you had to re-install BioPerl you need to install it under scripts/pfamscan/ folder. If you had to re-install T-Coffee and HMMER3 then you need to set the path variables called hmmer_dir and tcoffee_msa_dir inside scripts/dspred file.

# Downloading DSPRED 

1. cd to the directory where you want to install DSPRED
2. Type the following command<br/>
    git clone https://github.com/yusufzaferaydin/dspred.git
3. Download the model files for SVM under [models](models) directory of DSPRED from the following Google drive: https://drive.google.com/open?id=1hqzf2E4QWoKc9HCpnI9Ytmi3crwAAonw
4. Download the pfamscan directory under [scripts](scripts) directory of DSPRED from the following Google drive: https://drive.google.com/open?id=1CvKlPE2tek12x6FU_2KRzgU9kdV63kyl
5. cd dspred/scripts

# Setting the paths and running DSPRED

1. Set the paths in file named [dspred](scripts/dspred) 
2. ./dspred -h will print the help menu and parameter options of DSPRED
3. ./run_dspred includes example command line for running dspred
4. Outputs will be saved under the output directory (default: outputs/) in protein_id_dspred.out format.

# Feedback and comments

You are welcome to send comments to zafer.aydin@agu.edu.tr

