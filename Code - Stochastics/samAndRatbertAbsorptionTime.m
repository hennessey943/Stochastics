function tau=samAndRatbertAbsorptionTime(pr,ps)
%tau is the mean time it takes for Sam and Ratbert to meet in the maze
%pr is the probability for Ratbert to move
%ps is the probability for Sam to move
rr=1-pr;
rs=1-ps;
%Probability transition matrix R for Ratbert adjusted for new indexing scheme
R=zeros(16,16);
R(1,1)=rr;
R(1,2)=pr/2;
R(1,5)=pr/2;
R(2,1)=pr/3;
R(2,2)=rr;
R(2,3)=pr/3;
R(2,6)=pr/3;
R(3,2)=pr/2;
R(3,3)=rr;
R(3,4)=pr/2;
R(4,3)=pr/2;
R(4,4)=rr;
R(4,8)=pr/2;
R(5,1)=pr;
R(5,5)=rr;
R(6,2)=pr/2;
R(6,6)=rr;
R(6,7)=pr/2;
R(7,6)=pr/2;
R(7,7)=rr;
R(7,11)=pr/2;
R(8,4)=1;
R(9,10)=.5;
R(9,13)=.5;
R(10,9)=pr/2;
R(10,10)=rr;
R(10,11)=pr/2;
R(11,7)=pr/4;
R(11,10)=pr/4;
R(11,11)=rr;
R(11,12)=pr/4;
R(11,15)=pr/4;
R(12,11)=pr/2;
R(12,12)=rr;
R(12,16)=pr/2;
R(13,9)=pr/2;
R(13,13)=rr;
R(13,14)=pr/2;
R(14,13)=pr/2;
R(14,14)=rr;
R(14,15)=pr/2;
R(15,11)=pr/2;
R(15,14)=pr/2;
R(15,15)=rr;
R(16,12)=pr;
R(16,16)=rr;

%Transition probabilty matrix S for Sam
S=zeros(16,16);
S(1,1)=rs;
S(1,2)=ps/2;
S(1,5)=ps/2;
S(2,1)=ps/3;
S(2,2)=rs;
S(2,3)=ps/3;
S(2,6)=ps/3;
S(3,2)=ps/2;
S(3,3)=rs;
S(3,4)=ps/2;
S(4,3)=ps/2;
S(4,4)=rs;
S(4,8)=ps/2;
S(5,1)=ps;
S(5,5)=rs;
S(6,2)=ps/2;
S(6,6)=rs;
S(6,7)=ps/2;
S(7,6)=ps/2;
S(7,7)=rs;
S(7,11)=ps/2;
S(8,4)=1;
S(9,10)=.5;
S(9,13)=.5;
S(10,9)=ps/2;
S(10,10)=rs;
S(10,11)=ps/2;
S(11,7)=ps/4;
S(11,10)=ps/4;
S(11,11)=rs;
S(11,12)=ps/4;
S(11,15)=ps/4;
S(12,11)=ps/2;
S(12,12)=rs;
S(12,16)=ps/2;
S(13,9)=ps/2;
S(13,13)=rs;
S(13,14)=ps/2;
S(14,13)=ps/2;
S(14,14)=rs;
S(14,15)=ps/2;
S(15,11)=ps/2;
S(15,14)=ps/2;
S(15,15)=rs;
S(16,12)=ps;
S(16,16)=rs;

%Formation of transition matrix P for the movement of both Sam and Ratbert
P=zeros(256,256);
for i=1:256;
    for j=1:256;
        if mod(i,16)==0&&mod(j,16)~=0
            P(i,j)=R(ceil(i/16),ceil(j/16))*S(16,mod(j,16));
        elseif mod(i,16)==0&&mod(j,16)==0
            P(i,j)=R(ceil(i/16),ceil(j/16))*S(16,16);
        elseif mod(i,16)~=0&&mod(j,16)==0
            P(i,j)=R(ceil(i/16),ceil(j/16))*S(mod(i,16),16);
        else
            P(i,j)=R(ceil(i/16),ceil(j/16))*S(mod(i,16),mod(j,16));
        end
    end
end

 %Make the states (x,x) absorbing
 for i=1:256;
         if ceil(i/16)==mod(i,16)
             P(i,:)=zeros(1,256);
             P(i,i)=1;
         elseif mod(i,16)==0&&ceil(i/16)==16
             P(i,:)=zeros(1,256);
             P(i,i)=1;
         end
 end
 
 %Remove Absorbing states of P to form Q
for i=256:-1:1
    if P(i,i)==1
        P(i,:)=[];
        P(:,i)=[];
    end
end
Q=P;
%Form M
M=(eye(240,240)-Q)^(-1);
%f is the function which counts the number of steps
f=ones(1,240);
w=M*f';
%As Ratbert starts at the cheese (state 14) and Sam starts at state 16
tau=w(14*16-14);