#!/bin/bash

conda activate spyceatac

spyce-create \
      --setup_file "/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/sPYce/spyce.embed" \
      --k 6 \
      --norm "none" \
      --n_policy "remove" \
      --verbosity 2 \
      --data_path "/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/results/sPYce/" \
      --fig_path "/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC/plots/sPYce/" \
      --save_prefix "k6" \
      --correct_gc \
      --n_jobs 36
