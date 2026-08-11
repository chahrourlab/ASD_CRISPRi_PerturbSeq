#!/usr/bin/env python

import os
import re
import sys
import collections
import argparse
import itertools
import matplotlib
import glob
import math
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import scipy.sparse as sp_sparse
import scanpy as sc 
import scanpy.external as sce
from cnmf import cNMF
from collections import defaultdict
from scipy import sparse, io


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


output_directory = '/output_directory'
run_name = 'Neuron_Prod1_cNMF_run'
seed = 14 
cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
cnmf_obj.combine()

cnmf_obj.k_selection_plot(close_fig=False)
plt.savefig('K_Selection_plot.png')

density_threshold = 0.1
cnmf_obj.consensus(k=50, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
cnmf_obj.consensus(k=100, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
cnmf_obj.consensus(k=200, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
cnmf_obj.consensus(k=250, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
cnmf_obj.consensus(k=500, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)

