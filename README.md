# Readme file with instructions to reproduce the numerical experiments of the papers about Mean Field Games for the SIR model of Josu Doncel, Nicolas Gast and Bruno Gaujal

## Numerical experiments of the Netgcoop 2020 paper

To create Figure 2: execute the code of the notebook Figure2_doplot.ipynb

To create Figure 3: execute the code of the notebook Figure3_doplot.ipynb


## Notes:

- Some of the notebooks load the data of all the policies (mean-field equilibrium, 
N-players equilibrium, mean-field optimum and global optimum with N finite), 
which is contained in the folder data. To generate this data , execute the
code of the notebook compute_policy_t0.ipynb and to create new data, use the same
notebook

- The code of the file computation_mfg.py is used in most of the notebooks.

- Some npy files in the data folder must be uncompressed. 
