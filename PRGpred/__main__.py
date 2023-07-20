#!/usr/bin/python
"""
Title: Main Script for plant resistance gene prediction
Author: Naveen Duhan
Lab: KAABiL(Kaundal Artificial Intelligence & Advanced Bioinformatics Lab)
Version: 0.1
"""

import os, time
import argparse
import pandas as pd
from PRGpred import nn_prediction
from PRGpred import utils
from PRGpred import __version__

def DNN(fasta_file,output_dir,level):
    mess=""
    global result
    if level == "Phase1":
       
            enz1 = utils.preprocess(fasta_file, 'DPC', [0, 0], 400)
            phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
            phase1, result = nn_prediction.pred_Phase1(enz1, phase1_out)

    if level == "Phase2":

            enz1 = utils.preprocess(fasta_file, 'DPC', [0, 0], 400)
            phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
            phase1, result = nn_prediction.pred_Phase1(enz1, phase1_out)

            input2 = "%s/phase_2_input.fasta" % (output_dir)
            utils.fasta_process(fasta_file, input2, phase1)
            enz2 = utils.preprocess(input2, 'DPC', [0, 0], 400)
            phase2_out = '%s/Phase_2_dnn_log.txt' % (output_dir)
            result = nn_prediction.pred_Phase2(enz2, phase2_out)

    return result, mess

    
def main():
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

    start=time.time()
    parser =utils.argument_parser(version=__version__)
    options = parser.parse_args()
    fasta_file = options.fasta_file
    output_dir = options.output_dir
    output= options.output_file
    level = options.level

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    print("\n Created output directory\n")

    res, mess = DNN(fasta_file,output_dir,level)
    out ='%s/%s'% (output_dir,output)
    if mess=='' :
            res.to_csv(out, sep="\t",index=False)  
    else:
        print(mess)
        with open(out,'w') as fp:
            fp.write(mess)
            fp.close()


