#!/bin/bash          

##################################################################################
#
# Script to compare the threshold of the social optimum and the mean-field 
#   equilibrium when c_V varies from 0.01 to 1.21 with step 0.05
#
# Required file: comp_t_eq_and_t_opt.py
#
# Output: data.txt file,  where 
#       - the first column is the value of c_V
#	- the second column, the threshold value of the mean-field equilibrium
#	- the third column, the cost under the mean-field equilibrium vaccination
#	- the forth column, the threshold value of the social optimum
#	- the fifth column, the cost under the social optimum vaccination
#
##################################################################################

for a in 0.01 0.06 0.11 0.16 0.21 0.26 0.31 0.36 0.41 0.46 0.51 0.56 0.61 0.66 0.71 0.76 0.81 0.86 0.91 0.96 1.01 1.06 1.11 1.16 1.21
do
   python comp_t_eq_and_t_opt.py $a >> data.txt
done            
