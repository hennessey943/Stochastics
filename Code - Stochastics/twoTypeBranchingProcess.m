function [A,B]=twoTypeBranchingProcess(a,b,PA,indA,PB,indB,N,M)
%Monte Carlo Simulation of General Two Generational Branching Process
%a and b stands for the intial number of Type A and Type B members
%respectively
%A and B stands for the number of Type A and Type B members on the n-th day
%PA is the joint probability mass function for the offspring of Type A
%PB is the joint probability mass function for the offspring of Type B
%PA and PB are M+1xM+1 matrices whose row sums are the Type A offspring
%marginal distribution and whose column sums are the Type B offspring
%marginal distribution.
%Note that PA and PB do not need to be normalized for input
%indA and indB are Boolean inputs to simulate if type of offspring produced
%is independent for PA and PB respectively.
%The column index j=Type B offspring+1 
%and the row index i=Type A offspring+1
%N is the end time
%M is the maximum number of offspring from one individual
%dist** is the pmf for the Y^** outlined in the setup of the question
A=zeros(1,N);
B=zeros(1,N);
A(1)=a;
B(1)=b;
pAA=zeros(1,N);
pAB=zeros(1,N);
pBA=zeros(1,N);
pBB=zeros(1,N);
%Normalize PA
PA=PA/(sum(sum(PA)));
%Normalize PB
PB=PB/(sum(sum(PB)));
%Formation of Probability Vector for Type A offspring
PAvec=reshape(PA',1,(M+1)^2);
%Formation of Probability Vector for Type A offspring
PBvec=reshape(PB',1,(M+1)^2);
%Formation of marginal pmfs
margpAA=sum(PA);
margpAB=sum(PA');
margpBA=sum(PB);
margpBB=sum(PB');       
            
for n=1:N-1
    %Initialize offspring generators at epoch n
    pAA(n)=0;
    pAB(n)=0;
    pBA(n)=0;
    pBB(n)=0;
    yAB=zeros(1,A(n));
    yAA=zeros(1,A(n));
    yBB=zeros(1,B(n));
    yBA=zeros(1,B(n));
    %Calculate number of offspring produced by Type A's
    for k=1:A(n)
       
        if indA==0; %If type of offspring of A are dependent.
          
            %Simulation of offspring produced by Type A's
            x=rand();
            i=1;
            while x>0
                x=x-PAvec(i);
                i=i+1;
            end
            if i>M+1
                yAA(k)=mod(i,M+1);
            else
                yAA(k)=i-2;
            end
            
            j=1;
            t=1;
            while j>0
                j=i-1-t*(M+1);
                t=t+1;
            end
            yAB(k)=t-2;
        
        else %If type of offspring of A are independent
        
            %Simulation of Type A offspring produced by Type A's
            x=rand();
            i=1;
            while x>0
                x=x-margpAA(i);
                i=i+1;
            end
            yAA(k)=i-2;
        
            %Simulation of Type B offspring produced by Type A's
            z=rand();
            j=1;
            while z>0
                z=z-margpAB(j);
                j=j+1;
            end
            yAB(k)=j-2;
        end
        
        pAA(n)=yAA(k)+pAA(n);
        pAB(n)=yAB(k)+pAB(n);
    end        
    
    for k=1:B(n)
        
        
        if indB==0 %If type of offspring of B are not independent
 
            %Simulation of offspring produced by Type A's
            x=rand();
            i=1;
            while x>0
                x=x-PBvec(i);
                i=i+1;
            end
            if i>M+1
                yBA(k)=mod(i,M+1);
            else
                yBA(k)=i-2;
            end
            j=1;
            t=1;
            while j>0
                j=i-1-t*(M+1);
                t=t+1;
            end
            yBB(k)=t-2;
        
        else %If type of offspring of B are independent
                   
            %Simulation of Type A offspring produced by Type B's
            x=rand();
            i=1;
            while x>0
                x=x-margpBA(i);
                i=i+1;
            end
            yBA(k)=i-2;
        
            %Simulation of Type B offspring produced by Type A's
            z=rand();
            j=1;
            while z>0
                z=z-margpBB(j);
                j=j+1;
            end
            yBB(k)=j-2;
        end
        pBA(n)=yBA(k)+pBA(n);
        pBB(n)=yBB(k)+pBB(n);
    end
       
    %Stochastic Update Rule
    A(n+1)=pAA(n)+pBA(n);
    B(n+1)=pAB(n)+pBB(n);
end
% n=1:N;
% stairs(n,A);
% hold on
% stairs(n,B);
% legend('# of Type A individuals', '# of Type B individuals');
% title('Two Type Branching Process');

%Calculate mean number of offspring of each type produced

% %Type A/B produces a mean number of Type A&B offspring, as calculated from the
% %marginal pmfs;
% meanAA=0;
% meanAB=0;
% meanBA=0;
% meanBB=0;
% for i=1:M+1
%     meanAA=(i-1)*margpAA(i)+meanAA;
%     meanAB=(i-1)*margpAB(i)+meanAB;
%     meanBA=(i-1)*margpBA(i)+meanBA;
%     meanBB=(i-1)*margpBB(i)+meanBB;
% end
% 
% meanA=meanAA+meanBA;
% meanB=meanAB+meanBB;

