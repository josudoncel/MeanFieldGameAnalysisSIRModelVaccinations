function switching_curve

figure; hold on;
mat=dlmread('eqpol_t0_n5.txt')
ind=1;
N=5;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end

ternplot(x,y)


mat=dlmread('eqpol_t0_n10.txt')
ind=1;
N=10;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end

plot(x,y)
mat=dlmread('eqpol_t0_n20.txt')
ind=1;
N=20;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end

plot(x,y)

mat=dlmread('eqpol_t0_n30.txt')
ind=1;
N=30;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end
plot(x,y)

mat=dlmread('eqpol_t0_n40.txt')
ind=1;
N=40;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end
plot(x,y)

mat=dlmread('eqpol_t0_n50.txt')
ind=1;
N=50;
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y(ind)=(i-1)/N;
       x(ind)=(j-1)/N;
       ind=ind+1;
    end
  end
end
plot(x,y)

mat=dlmread('eqpol_t0_mf.txt')
ind=1;
N=size(mat,1);
for j=1:N+1
  for i=2:N-j+2
    if mat(j,i)==10 && mat(j,i-1)==0
       y_(ind)=(i-1)/N;
       x_(ind)=(j-1)/N;
       ind=ind+1;
	break;
    end
  end
end
x_
y_

plot(x_,y_)

end
