function cost_Nfinite_mfopt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab code to generate the figure where we compare the cost
%   for the optimum with N objects and the cost of the mean-field
%    optimum. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model parameters
rho=36.5;
gamma=73;
theta=10;
T=1;
cI=36.5;
cV=0.5;

figure; 

x_N=[10 20 30 50 70 100];
ind=1;
for N=x_N
   res1(ind)=compute_Nfinite_mS_mI(N,rho,gamma,theta,T,cI,cV,0.4,0.4);
   ind=ind+1;
end
res1_mf=compute_meanfield_mS_mI(rho,gamma,theta,T,cI,cV,0.4,0.4);
plot(x_N,res1,x_N,ones(1,length(x_N))*res1_mf,x_N,0.6818-0.39./x_N)
title('Optimal cost')
xlabel('N')
legend('N finite','Mean field','0.6818-0.39/N')
end

function res=compute_Nfinite_mS_mI(N,rho,gamma,theta,T_,cI,cV,m_S,m_I)

C=max([rho,gamma,theta]);
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
res=J(m_S*N+1,m_I*N+1,1)/(C*N);
end



function cost=compute_meanfield_mS_mI(rho,gamma,theta,T_,cI,cV,m_S,m_I)

% the mean-field solution
h=0.0001;
H=T_/h;

%%%%%%%%%%
INF(1)=m_S;
SUS(1)=m_I;
for t=1:H
    SUS(t+1)=SUS(t)+h*(-INF(t)*SUS(t)*gamma);
    INF(t+1)=INF(t)+h*(INF(t)*SUS(t)*gamma-INF(t)*rho);
end
cost_novac=sum(cI*INF.*h); cost=cost_novac;
bestpol_mf_=0;

for tc=2:H
  for t=1:H
    if t<tc
       vac(t)=theta;
    else
       vac(t)=0;
    end
    SUS(t+1)=SUS(t)+h*(-INF(t)*SUS(t)*gamma-vac(t)*SUS(t));
    INF(t+1)=INF(t)+h*(INF(t)*SUS(t)*gamma-INF(t)*rho);
  end
  vac(t+1)=0;
  cost_=sum(SUS.*vac.*cV.*h+cI*INF.*h);
  if cost_<cost
     bestpol_mf_=theta;
     cost=cost_;
  end
end


end
