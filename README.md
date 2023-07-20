PRGpred: Deep neural network based Plant Resistance Gene Prediction
==========

PRGpred classify Plant resistance genes in 8 classes: CNL, KIN, LYK, LECRK, RLP, RLK, TIR, TNL.

Input
-----

The PRGpred server requires amino acid sequence(s) in fasta format, and does not support nucleic acid sequences.

Development Environment and Prerequisite
----------------------------------------
This source code was developed in Linux, and has been tested on Linux and OS X. The only prerequisite is to have Python 3.7 or above installed. 

Installation
------------
    The installation of PRGpred can be done in two ways:

    1. Create a dedicated miniconda3 environment
    
        Download PRGpred from:
            
            https://bioinfo.usu.edu/PRGpred/download/PRGpred.tar.gz

        Download the Miniconda installer: 
            
            (https://docs.conda.io/en/latest/miniconda.html#linux-installers)

        Extract the downloaded file:

            tar -xvzf PRGpred.tar.gz

        cd PRGpred
        
        Create and activate a conda environment

            conda env create -f environment.yml

            conda activate PRGpred

            pip3 install .
    
    2. Intall using system Python3

        Download PRGpred from: 

             https://bioinfo.usu.edu/PRGpred/download/PRGpred.tar.gz

        Extract the downloaded file:

            tar -xvzf PRGpred.tar.gz

        cd PRGpred

            pip3 install .

Input and Execution
-------------------

    PRGpred will be installed under the name 'PRGpred' with three possible arguments:

    -i, --fasta_file. Protein sequence input in FASTA format.
    -od, --output_dir. Output directory name.
    -l, --level. 

Output Files
------------

The output directory will have three files in tab-delimted format:

1. Phase1_dnn_log.txt
    * 1st column: Protein ID.
    * 2nd column: Prediction.
    * 3rd-4th column: Probability for Rgenes and Non-Rgenes.
    
2. Phase2_dnn.log.txt
    * 1st column: Protein ID.
    * 2nd column: Prediction.
    * 3rd-14th column: Probability for each class.

Queries and Contact
----------------------

Written by Naveen Duhan (naveen.duhan@usu.edu),

Kaundal Bioinformatics Lab, Utah State University,

Released under the terms of GNU General Public Licence v3

In case of technical problems (bugs etc.) please contact Naveen Duhan naveen.duhan@usu.edu.

For any Questions on the scientific aspects of the PRGpred method please contact:

Rakesh Kaundal, rkaundal@usu.edu

Naveen Duhan, naveen.duhan@usu.edu



