function cost_Nfinite_mfe()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab code to generate the figure where we compare the cost
%   for the N-player equilibrium with the cost of the mean-field
%    equilibrium. The values are obtained from 
%    SIR-MFG-vs-Nplayergame.ipynb
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the cost of the mean-field equilibrium and the N-players equilibrium
mfe_cost=0.6824;
eq_cost=[0.6345,0.6566,0.6694,0.6737,0.6772];

vN=[5,10,20,30,50];
plot(vN,ones(1,length(vN))*mfe_cost)
hold on;
plot(vN,eq_cost,'*-');
plot(vN,mfe_cost-0.26./vN)
end
