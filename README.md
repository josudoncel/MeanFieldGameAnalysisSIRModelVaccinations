# Readme file with instructions to reproduce the numerical experiments of the papers about Mean Field Games for the SIR model of Josu Doncel, Nicolas Gast and Bruno Gaujal

## Numerical experiments of [1]

To create Figure 2: execute the code of the notebook Figure2_doplot.ipynb

To create Figure 3: execute the code of the notebook Figure3_doplot.ipynb

## Numerical experiments of [2]

To create Figure 2: execute the code of the notebook Figure2_doplot.ipynb

To create Figure 3: execute the code of the notebook Figure3_doplot.ipynb

To create Figure 4: execute the code of the notebook Figure4_doplot.ipynb

To create Figure 5: execute the code of the notebook Figure5_doplot.ipynb

To create Figure 6: execute the code of the notebook Figure6_doplot.ipynb

To create Figure 7: execute the code of the notebook Figure7_doplot.ipynb

To create Figure 8 and Figure 9: execute the code of the notebook Figure8_and_Figure9_doplot.ipynb

## Notes:

- Some of the notebooks load the data of all the policies (mean-field equilibrium, 
N-players equilibrium, mean-field optimum and global optimum with N finite), 
which is contained in the folder data. To generate this data , execute the
code of the notebook compute_policy_t0.ipynb and to create new data, use the same
notebook

- The code of the file computation_mfg.py is used in most of the notebooks.

- Some npy files in the data folder must be uncompressed. 

## References

[1] Josu Doncel and Nicolas Gast and Bruno Gaujal, "Virus Propagation in a Large Population: Mean
Field Equilibrium versus Social Optimum". 10th International Conference on NETwork Games, COntrol and OPtimization (NETGCOOP). Cargèse, Corsica, France. Sept 2021.

[2] Josu Doncel and Nicolas Gast and Bruno Gaujal, "A Mean Field Game Analysis of SIR Dynamics with Vaccination". Submitted to Probability in the Engineering and Informational Sciences.
