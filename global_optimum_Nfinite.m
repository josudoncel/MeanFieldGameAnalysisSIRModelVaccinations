function global_optimum_Nfinite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab code to calculate the cost of the global optimum with 
%   N objects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model parameters
rho=36.5;
gamma=73;
theta=10;
T=1;
cI=36.5;
cV=0.5;

cost_opt=compute_cost_eq_Nfinite(10,rho,gamma,theta,T,cI,cV);
save('data/Nfinite_opt_cost_N10')

cost_opt=compute_cost_eq_Nfinite(20,rho,gamma,theta,T,cI,cV);
save('data/Nfinite_opt_cost_N20')

end



function cost_opt=compute_cost_eq_Nfinite(N,rho,gamma,theta,T_,cI,cV)

C=100*max([rho,gamma,theta]);
T=C*T_*N;


J=zeros(N+1,N+1,T);
bestpol=zeros(N+1,N+1,T);

for t=T:-1:2
    % CASE: no susceptibles, no infected
    mS=1; mI=1;
    J(mS,mI,t-1)=J(mS,mI,t);
    bestpol(mS,mI,t-1)=0;
    % CASE: no susceptibles, but infected
    for mI=2:N+1
      J(mS,mI,t-1)=cI*(mI-1)/N+rho/C*(mI-1)/N*J(mS,mI-1,t)+(1-rho/C*(mI-1)/N)*J(mS,mI,t);
      bestpol(mS,mI,t-1)=0;
    end
    for mS=2:N+1
          % CASE: no infected, but susceptibles
          mI=1;
	  bestpol(mS,mI,t-1)=0;
          J(mS,mI,t-1)=J(mS,mI,t);
	  % CASE: infected and susceptibles
          for mI=2:N-mS+2           
	          if cV+(1/C)*J(mS-1,mI,t)-(1/C)*J(mS,mI,t)<0		
		          J(mS,mI,t-1)=cV*theta*(mS-1)/N+cI*(mI-1)/N+(theta/C)*(mS-1)/N*J(mS-1,mI,t)...
					              +(rho/C)*(mI-1)/N*J(mS,mI-1,t)...
					              +(gamma/C*(mI-1)/N*(mS-1)/N)*J(mS-1,mI+1,t)...
					              +(1-theta/C*(mS-1)/N-rho/C*(mI-1)/N-gamma/C*(mI-1)/N*(mS-1)/N)*J(mS,mI,t);
 		          bestpol(mS,mI,t-1)=theta;
	          else
		          J(mS,mI,t-1)=cI*(mI-1)/N+rho/C*(mI-1)/N*J(mS,mI-1,t)+gamma/C*(mI-1)/N*(mS-1)/N*J(mS-1,mI+1,t)...
					                               +(1-rho/C*(mI-1)/N-gamma/C*(mI-1)/N*(mS-1)/N)*J(mS,mI,t);
 		          bestpol(mS,mI,t-1)=0;
	          end
          end
    end
end


ind=1;
for i=2:N+1
  for j=2:N-i+2
     cost_opt(i-1,j-1)=J(i,j,1)/(N*C);
  end
end


end


