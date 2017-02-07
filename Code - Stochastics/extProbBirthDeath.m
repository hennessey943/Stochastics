function [u,w,lambda,mu]=extProbBirthDeath(K,lmb,m)
%This function calculates the probability of and the expected time until
%extinction of an asexual population with carrying capacity K 

%Form transition-rate matrix A
A=zeros(K+1,K+1);
A(1,1)=1;
for i=2:K+1
    A(i,i-1)=m*(i-1);
    A(i,i+1)=lmb*(i-1)*(1-(i-1)/K);
    A(i,i)=-m*(i-1)-lmb*(i-1)*(1-(i-1)/K);
end
A=A(1:K+1,1:K+1);
%Calculate the absorption probability vector u whose entries correspond to
%absorption from beginning population corresponding to the entry#-1

u=-A(2:K+1,2:K+1)\A(2:K+1,1);
u=u(2);

%Calculate the lambda_n
lambda=zeros(1,K+1);
for n=1:K+1
    lambda(n)=lmb*(n-1)*(1-(n-1)/K);
end

%Calculate the mu_n
mu=zeros(1,K);
for n=1:K
    mu(n)=m*n;
end


%Calculate the mean time to absorption
prod=1;
sum=0;
lambda(K+1)=1;
for i=1:K
    for j=1:i
        prod=prod*lambda(j+1)/mu(j);
    end
    sum=sum+prod;
end

w=sum;

%Iterate over a number of values of K with suitable values of lambda and
%mu, plot u and w
U=zeros(1,K);
U(1,1)=1;
W=zeros(1,K);
W(1,1)=1;
for k=2:K
    %Form transition-rate matrix A
    A=zeros(k+1,k+1);
    A(1,1)=1;
    for i=2:k+1
        A(i,i-1)=m*(i-1);
        A(i,i+1)=lmb*(i-1)*(1-(i-1)/k);
        A(i,i)=-m*(i-1)-lmb*(i-1)*(1-(i-1)/k);
    end
    A=A(1:k+1,1:k+1);
    %Calculate the absorption probability vector u whose entries correspond to
    %absorption from beginning population corresponding to the entry#-1

    absprob=-A(2:k+1,2:k+1)\A(2:k+1,1);
    U(k)=1-abs(1-abs(absprob(2)));
    

    %Calculate expected time until extinction
    %Calculate the lambda_n
    lambda=zeros(1,k+1);
    for n=1:k+1
        lambda(n)=lmb*(n-1)*(1-(n-1)/k);
    end

    %Calculate the mu_n
    mu=zeros(1,k);
    for n=1:k
        mu(n)=m*n;
    end


    %Calculate the mean time to absorption
    prod=1;
    sum=0;
    lambda(k+1)=1;
    for i=1:k
        for j=1:i
            prod=prod*lambda(j+1)/mu(j);
        end
        sum=sum+prod;
    end
    W(k)=sum;
end


K=1:K;
plot(K,U)
figure
plot(K,W)


        