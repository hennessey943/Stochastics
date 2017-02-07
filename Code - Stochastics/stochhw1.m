p=.3;
q=.15;
r=.5;

P=[1,0,0,0;q,(1-p*r)^2*(1-q),2*(1-p*r)*(p*r)*(1-q),(1-q)*(p*r)^2;
    q^2,2*q*(1-q)*(1-p*r)^2,(1-p*r)^2*(1-q)^2+4*q*p*r*(1-q)*(1-p*r),(1-q)^2*(1-(1-p*r)^2)+2*q*(1-q)*(p*r)^2;
    q^3,3*q^2*(1-q)*(1-p*r)^2,3*(1-q)^2*q*(1-p*r)^2+6*q^2*p*r*(1-q)*(1-p*r),(1-q)^3+3*q*(1-q)^2*(1-(1-p*r)^2)+3*q^2*(1-q)*(p*r)^2]

%Expected Value of P using weighted sums
prob=P;
expect=zeroes(1,n);
for i=i:n
    prob=prob*P;
    expect(i)=0*prob(2,1)+prob(2,2)+2*prob(2,3)+3*prob(2,4);
end
plot(1:n,expect)