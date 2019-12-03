function switchingcurve_Nfinite_mfopt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab code to generate the figure where we compare the switching curve
%   for the optimum with N objects and the switching curve of the mean-field
%    optimum. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('ternaryplotfiles')

% model parameters
rho=36.5;
gamma=73;
theta=10;
T=1;
cI=36.5;
cV=0.5;



compute_Nfinite(10,rho,gamma,theta,T,cI,cV);
%pause
compute_Nfinite(20,rho,gamma,theta,T,cI,cV);
compute_Nfinite(50,rho,gamma,theta,T,cI,cV);
compute_Nfinite(100,rho,gamma,theta,T,cI,cV);
compute_Nfinite(150,rho,gamma,theta,T,cI,cV);
compute_Nfinite(200,rho,gamma,theta,T,cI,cV);
%pause
ternlabel('Proportion of Infected Population',...
    'Proportion of Susceptible Population',...
    'Proportion of Recovered Population')
compute_meanfield(500,rho,gamma,theta,T,cI,cV);
%legend('N=10','N=50','N=100','N=150','Mean Field')
%xlabel('m_S')
%ylabel('m_I')
%title('Switching Curve')
%savefig('switching_curve3.fig')
end

function compute_Nfinite(N,rho,gamma,theta,T_,cI,cV)

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
%t/T
%J(:,:,t-1)
%pause
end

%for i=1:N+1
%  for j=1:N-i+2
%    if bestpol(i,j,1)==theta
%      plot((i-1)/N,(j-1)/N,'g*')
%    else
%      plot((i-1)/N,(j-1)/N,'r*')
%    end
%  end
%end

ind=1;
for j=2:N+1
  for i=2:N-j+2
    if bestpol(i,j,1)==theta && bestpol(i-1,j,1)==0
       x_(ind)=(i-1)/N;
       y_(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end

ind2=2;
x(1)=x_(1); y(1)=y_(1);
for c=2:ind-1
   if x_(c-1)~=x_(c) && y_(c-1)~=y_(c)
     x(ind2)=x_(c);
     y(ind2)=y_(c);
     ind2=ind2+1;
   end
end

ternplot(x,y)
vx=[x(end) x_(end)];
vy=[y(end) y_(end)];
hold on;
ternplot(vx,vy);
vx=[x(1) 1-y(1)];
vy=[y(1) y(1)];
ternplot(vx,vy);
end



function compute_meanfield(M,rho,gamma,theta,T_,cI,cV)

% the mean-field solution
val=zeros(M+1,M+1);

h=0.001;
H=T_/h;

%%%%%%%%%%
%INF(1)=0.4;
%SUS(1)=0.4;
%for t=1:H
%    SUS(t+1)=SUS(t)+h*(-INF(t)*SUS(t)*gamma);
%    INF(t+1)=INF(t)+h*(INF(t)*SUS(t)*gamma-INF(t)*rho);
%end
%cost_novac=sum(cI*INF.*h); cost=cost_novac;
%bestpol_mf_=0;
%
%for tc=2:H
%  for t=1:H
%    if t<tc
%       vac(t)=theta;
%    else
%       vac(t)=0;
%    end
%    SUS(t+1)=SUS(t)+h*(-INF(t)*SUS(t)*gamma-vac(t)*SUS(t));
%    INF(t+1)=INF(t)+h*(INF(t)*SUS(t)*gamma-INF(t)*rho);
%  end
%  vac(t+1)=0;
%  cost_=sum(SUS.*vac.*cV.*h+cI*INF.*h);
%  if cost_<cost%cost_novac
%     bestpol_mf_=theta;
%     cost=cost_;
%     %break;
%  end
%end
%cost
%pause

for mS=1:M+1
for mI=1:M-mS+2
INF(1)=(mI-1)/M;
SUS(1)=(mS-1)/M;

for t=1:H
    SUS(t+1)=SUS(t)+h*(-INF(t)*SUS(t)*gamma);
    INF(t+1)=INF(t)+h*(INF(t)*SUS(t)*gamma-INF(t)*rho);
end
cost_novac=sum(cI*INF.*h);
bestpol_mf(mS,mI)=0;

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
  if cost_<cost_novac
     bestpol_mf(mS,mI)=theta;
     %break;
  end
end



%[val(mS,mI) pos(mS,mI)]=min(cost);
%if pos(mS,mI)>1
%   bestpol_mf(mS,mI)=theta; 
%else
%   bestpol_mf(mS,mI)=0;
%end
end
end
 
%figure; hold on;
%
%for i=1:M+1
%  for j=1:M-i+2
%    if bestpol_mf(i,j)==theta
%      plot((i-1)/M,(j-1)/M,'g*')
%    else
%      plot((i-1)/M,(j-1)/M,'r*')
%    end
%  end
%end
%text(0.7,0.7,'VAC MAX','Color','green','FontSize',14)
%text(0.7,0.6,'NO VAC','Color','red','FontSize',14)
%title('Mean field')
%xlabel('m_S')
%ylabel('m_I')
%rectangle('Position',[0.67,0.55,0.25,0.2])

ind=1;
for j=2:M+1
  for i=2:M-j+2
    if bestpol_mf(i,j,1)==theta && bestpol_mf(i-1,j,1)==0
       x_(ind)=(i-1)/M;
       y_(ind)=(j-1)/M;
       ind=ind+1;
       
    end
  end
end

ind2=2;
x(1)=x_(1); y(1)=y_(1);
for c=2:ind-1
   if x_(c-1)~=x_(c) && y_(c-1)~=y_(c)
     x(ind2)=x_(c);
     y(ind2)=y_(c);
     ind2=ind2+1;
   end
end


ternplot(x,y)
vx=[x(1) 1-y(1)];
vy=[y(1) y(1)];
ternplot(vx,vy);


end
