#####

This repository contains all the codes to reproduce 
the numerical experiments of the paper A Mean Field Game 
Analysis for SIR Models with Vaccinations.

## Experiment 1: Dynamics of the population under social optimum and mean-field equilibrium

In this experiment, we fix the parameters and we observe evolution of the dynamics of the system for the social optimum and the mean-field equilibrium. 

### Required Software

* Matlab 2014 or higher
* Octave 3.8.1 or higher

### Running the examples

Execute region_plot.m file

### Output 

A figure with three areas that determine when the mean-field equilibrium and the social optimum vaccinate with maximum rate and with minimum rate and with a dashed line (resp. solid line) the evolution of the mean-field equilibrium (resp. social optimum).

## Experiment 2: Comparison of the thresholds of the social optimum and of the mean-field equilibrium when c_V varies

### Required Software

* Python 2.7.6

### Running the examples

Execute scr.sh in a terminal

### Output

A file called data.txt with 5 columns and in each row the output obtained for each c_V. In the first column, there is the value of c_V, in the second and the third one, the value of the threshold and the cost of the mean-field equilibrium and in the forth and the fifth one, the value of the threshold and the cost of the social optimum. From this data, one can create figures comparing the thresholds and the costs of both policies with the desired software.
