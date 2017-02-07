function t=expectedtimefullness(p,C)

%Chosen Parameter Values
p=p;
r=1-p;

%Transition Probability Matrix 
P=[r,p/2,0,0,p/2,0,0,0,0,0,0,0,0,0,0,0;
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,r,p/2,0,0,0,0,0,p/2,0,0,0,0,0,0;
    0,0,p,r,0,0,0,0,0,0,0,0,0,0,0,0;
    p/2,0,0,0,r,0,0,0,p/2,0,0,0,0,0,0,0;
    0,0,0,0,0,r,p/2,0,0,p/2,0,0,0,0,0,0;
    0,0,p/4,0,0,p/4,r,p/4,0,0,p/4,0,0,0,0,0;
    0,0,0,0,0,0,p/2,r,0,0,0,p/2,0,0,0,0;
    0,0,0,0,p/3,0,0,0,r,p/3,0,0,p/3,0,0,0;
    0,0,0,0,0,p/2,0,0,p/2,r,0,0,0,0,0,0;
    0,0,0,0,0,0,p/2,0,0,0,r,0,0,0,p/2,0;
    0,0,0,0,0,0,0,p/2,0,0,0,r,0,0,0,p/2;
    0,0,0,0,0,0,0,0,p/2,0,0,0,r,p/2,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,p,r,0,0;
    0,0,0,0,0,0,0,0,0,0,1/2,0,0,0,0,1/2;
    0,0,0,0,0,0,0,0,0,0,0,p/2,0,0,p/2,r];
    %Eigenvalue decomposition to find the stationary vector
[V,D,W]=eig(P);
x=1;
t=1;
while abs(x)>=.0001 
    x=1-D(t,t);
    t=t+1;
end
pi=W(:,t-1)/sum(W(:,t-1));

%Expected Time to C pieces of cheese eaten (fullness)
t=absorptiontime(p)+(C-1)/pi(12);