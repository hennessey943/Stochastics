function stoch_disease_model_scenario1(p,q,r,n)

%Probability Transition Matrix for Scenario 1
P1=zeros(16,16);
P1(1,1)=1;
P1(2,1)=q;
P1(2,2)=1-q;
P1(3,1)=q^2;
P1(3,2)=2*q*(1-q);
P1(3,3)=(1-q)^2;
P1(4,1)=q^3;
P1(4,2)=3*q^2*(1-q);
P1(4,3)=3*q*(1-q)^2;
P1(4,4)=(1-q)^3;
P1(5,5)=1;
P1(6,5)=q;
P1(6,6)=1/3*(1-q)*(3-2*p);
P1(6,7)=2/3*p*(1-q);
P1(7,5)=q^2;
P1(7,6)=2/3*q*(1-q)+4/3*q*(1-q)*(1-p);
P1(7,7)=1/3*(1-q)^2+4/3*q*(1-q)*p+2/3*(1-q)^2*(1-p);
P1(7,8)=2/3*(1-q)^2*p;
P1(8,5)=q^3;
P1(8,6)=2*q^2*(1-q)*(1-p)+q^2*(1-q);
P1(8,7)=(1-q)*q*(3-3*q+p*(4*q-2));
P1(8,8)=(1-q)^3+2*p*q*(1-q)^2;
P1(9,9)=1;
P1(10,9)=q;
P1(10,10)=1/3*(1-q)*(1-p)*(3-p);
P1(10,11)=2/3*(1-q)*p*(2-p);
P1(10,12)=1/3*(1-q)*p^2;
P1(11,9)=q^2;
P1(11,10)=4/3*q*(1-q)*(1-p)+2/3*q*(1-q)*(1-p)^2;
P1(11,11)=4/3*q*(1-q)*p+2/3*(1-q)^2*(1-p)+4/3*q*(1-q)*(1-p)*p+1/3*(1-q)^2*(1-p)^2;
P1(11,12)=2/3*(1-q)^2*p+2/3*q*(1-q)*p^2+1/3*(1-q)^2*(1-(1-p)^2);
P1(12,9)=q^3;
P1(12,10)=q^2*(1-q)*(1-p)*(3-p);
P1(12,11)=2*q*(1-q)^2*(1-p)+2*q^2*(1-q)*(1-p)*p+2*q^2*(1-q)*p+q*(1-q)^2*(1-p)^2;
P1(12,12)=(1-q)^3+2*q*(1-q)^2*p+q*(1-q)^2*(1-(1-p)^2)+q^2*(1-q)*p^2;
P1(13,13)=1;
P1(14,13)=q;
P1(14,14)=(1-p)^2*(1-q);
P1(14,15)=2*(1-p)*(p)*(1-q);
P1(14,16)=(1-q)*(p)^2;
P1(15,13)=q^2;
P1(15,14)=2*q*(1-q)*(1-p)^2;
P1(15,15)=(1-p)^2*(1-q)^2+4*q*p*(1-q)*(1-p);
P1(15,16)=(1-q)^2*(1-(1-p)^2)+2*q*(1-q)*(p)^2;
P1(16,13)=q^3;
P1(16,14)=3*q^2*(1-q)*(1-p)^2;
P1(16,15)=3*(1-q)^2*q*(1-p)^2+6*q^2*p*(1-q)*(1-p);
P1(16,16)=(1-q)^3+3*q*(1-q)^2*(1-(1-p)^2)+3*q^2*(1-q)*(p)^2;

rvec=zeros(1,16);
rvec(2)=(1-r)^3;
rvec(6)=3*r*(1-r)^2;
rvec(10)=3*r^2*(1-r);
rvec(14)=r^3;

%Expected Value of P1 using weighted sums%
prob1=eye(16);
expect1=zeros(1,n);
for i=1:n
    prob1=prob1*P1;
    vec=rvec*prob1;
    expect1(i)=0*vec(1)+vec(2)+2*vec(3)+3*vec(4)+0*vec(5)+vec(6)+2*vec(7)+3*vec(8)+0*vec(9)+vec(10)+2*vec(11)+3*vec(12)+0*vec(13)+vec(14)+2*vec(15)+3*vec(16);
end
%plot(1:n,expect1);

%Monte Cristo simulation of scenario 1

x=rand;
