function [qA,qB]=detExtinctProb(PA,PB,N)
% %We use the iterative formula to determine the extinction probability of a
% %two type branching process markov chain with joint pmfs PA and PB
% %N is the  maximum epoch of interest
%qA is the extinction probability of the scenario beginning with A0=1
%qB is the extinction probability of the scenario beginning with B0=1

%Iterate the iteration and find extinction probability vectors
for i=1:N
    [phiA,phiB]=iterativePGF(PA,PB,i);
    qA(i)=phiA;
    qB(i)=phiB;
end

%Need to iterate as letting s=t=0 in pgf function only gives first
%entry of matrix
function [phiA,phiB]=iterativePGF(PA,PB,m)
phiA=0;
phiB=0;
for j=1:m
    [phiA,phiB]=pgf(PA,PB,phiA,phiB);
end

%probability generating function for any s and t
function [phiA,phiB]=pgf(PA,PB,s,t)
[~,M]=size(PA);
M=M-1;
phiA=0;
phiB=0;
for k=1:M+1
    for l=1:M+1
        phiA=phiA+s^(k-1)*t^(l-1)*PA(k,l);
        phiB=phiB+s^(k-1)*t^(l-1)*PB(k,l);
    end
end




        