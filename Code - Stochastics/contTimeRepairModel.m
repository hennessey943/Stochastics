function [lambda,mu,theta,pi]=contTimeRepairModel(N,M,R,Tr,muOn,muOff)
%This function calculates the stationary distribution for a continuous time
%machine repair model. 
%N is the total number of machines in the scenario
%M is the maximum number of machines that can work at any time t
%R is the maximum number of machines that can be repaired at any time t
%Tr is the average amount of time it takes to repair a machine
%muOn is the rate a machine breaks when it is in use
%muOff is the rate a machine breks when it is not in use

%Calculate the lambda_n
lambda=zeros(1,N);
for n=1:N-R+1
    lambda(n)=R*Tr;
end
for n=N-R+2:N
    lambda(n)=Tr*(N-n+1);
end

%Calculate the mu_n
mu=zeros(1,N);
for n=1:M
    mu(n)=muOn*n;
end
for n=M+1:N
    mu(n)=muOn*M+muOff*(n-M);
end

%Calculate the theta_n
theta=zeros(1,N+1);
theta(1)=1;
prod=1;
for k=2:N+1
    for j=2:k
        prod=prod*lambda(k-1)/mu(k-1);
    end
theta(k)=prod;
end

%Calculate the pi vector
sum=0;
pi=zeros(1,N+1);
for k=1:N+1
    sum=sum+theta(k);
end
for j=1:N+1
    pi(j)=theta(j)/sum;
end

%Plot the stationary distribution
t=0:N;
scatter(t,pi)