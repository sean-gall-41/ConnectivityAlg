#!./venv/bin/python

import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# define global path vars
ROOT_DIR = '/home/seang/Dev/Git/UT_Austin/ConnectivityAlg/'
DATA_DIR = ROOT_DIR + 'data/'
DATA_INPUT_DIR = DATA_DIR + 'inputs/'
CELL_PARAM_FILE = DATA_INPUT_DIR + 'cell_params.json'
CON_PARAM_FILE = DATA_INPUT_DIR + 'con_params_test_1.json'


# open input param files for cell/connectivity info for analysis
# NOTE: only doing golgi so far (05/01/22)
with open(CELL_PARAM_FILE, 'r') as cell_param_file:
    global cell_param_dict
    cell_param_dict = json.load(cell_param_file)['GO'];

with open(CON_PARAM_FILE, 'r') as con_param_file:
    global con_param_dict
    con_param_dict = json.load(con_param_file)['GOGO']

num_tot_src_cell = int(cell_param_dict['numTotCells'])
num_tot_dest_cell = num_tot_src_cell
num_con_to_make = int(con_param_dict['numConToMake'])

# obtain outputs from connectivity algorithm, run analysis on them
with open(sys.argv[1], 'rb') as f:
    data_array = np.fromfile(f, dtype=np.intc)

    # get slices from total data input
    dest_num_con_arr = data_array[0:num_tot_dest_cell]
    src_ind_offset = num_tot_dest_cell * (1 + num_con_to_make)
    src_num_con_arr = data_array[src_ind_offset:(src_ind_offset + num_tot_src_cell)]
   
    # compute histogram of convergence array
    _ = plt.hist(
            dest_num_con_arr,
            bins=(num_con_to_make),
            weights=np.full(num_tot_dest_cell,
                1 / num_tot_dest_cell)
            ) 

    plt.title("golgi-golgi convergence")
    plt.show()
