#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

with open(sys.argv[1], 'rb') as f:
    data_array = np.fromfile(f, dtype=np.intc)
    print(data_array)
