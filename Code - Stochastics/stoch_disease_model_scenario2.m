function stoch_disease_model_scenario2(p,q,r,n)

%Probability Transition Matrix for Scenario 2
P2=[1,0,0,0;
    q,(1-p*r)^2*(1-q),2*(1-p*r)*(p*r)*(1-q),(1-q)*(p*r)^2;
    q^2,2*q*(1-q)*(1-p*r)^2,(1-p*r)^2*(1-q)^2+4*q*p*r*(1-q)*(1-p*r),(1-q)^2*(1-(1-p*r)^2)+2*q*(1-q)*(p*r)^2;
    q^3,3*q^2*(1-q)*(1-p*r)^2,3*(1-q)^2*q*(1-p*r)^2+6*q^2*p*r*(1-q)*(1-p*r),(1-q)^3+3*q*(1-q)^2*(1-(1-p*r)^2)+3*q^2*(1-q)*(p*r)^2];


%Expected Value of P2 using weighted sums
prob2=eye(4);
expect=zeros(1,n);
for i=1:n
    prob2=prob2*P2;
    expect(i)=0*prob2(2,1)+prob2(2,2)+2*prob2(2,3)+3*prob2(2,4);
end
plot(1:n,expect);
title('Expected Number of Infections at Day n: Scenario 2, p=.3,q=.15,r=.6');
xlabel('Day');
ylabel('Number of Infections Expected');

Monte Carlo Simulation of P2 scenario

% t=2;
% inf=zeros(1,n);
% inf(1)=1;
% i=0;
% while t<n
% x=rand;
%     while x>0
%        x=x-P2(inf(t-1)+1,i+1);
%        i=i+1;
%     end
% inf(t)=i-1;    
% t=t+1;
% i=0;
% end
% 
% stairs(1:1:n,inf)
% title('Disease Model Scenario 2');
% 
% xlabel('Day');
% ylabel('Number of Infections');



