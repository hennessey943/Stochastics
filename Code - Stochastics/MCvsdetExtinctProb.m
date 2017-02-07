function MCvsdetExtinctProb(PA,PB,N)
%we compare the deterministic calculation with the average of the 
%monte carlo simulations which go extinct
%for ease of calculation, we let 0<N<=20
%normA and normB are the normed differences between the estimated and
%deterministic solutions

%Calculate deterministic values
[qA,qB]=detExtinctProb(PA,PB,N);


[~,M]=size(PA);
M=M-1;
estqA=zeros(1,N);
estqB=zeros(1,N);
%Iterate over the simulations
for j=1:N
    %Initialize the counters
    countA=0;
    countB=0;
    %Run 500 Monte Carlo simulations and find % that go extinct to estimate qA
    for i=1:500
        [A,B]=twoTypeBranchingProcess(1,0,PA,0,PB,0,j,M);
        if A(j)==0&&B(j)==0
            countA=countA+1;
        end
    end
    estqA(j)=countA/500;

    for i=1:500
        [A,B]=twoTypeBranchingProcess(0,1,PA,0,PB,0,j,M);
        if A(j)==0&&B(j)==0
            countB=countB+1;
        end
    end
    estqB(j)=countB/500;
end
x=1:N;
t=0:N-1;
plot(x,qA,t,estqA);
hold off
figure

plot(x,qB,t,estqB);


