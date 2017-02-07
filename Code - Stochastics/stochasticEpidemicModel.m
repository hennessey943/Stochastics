function [timeMat,susceptMat,infectMat,recovMat]=stochasticEpidemicModel(S0,I0,lambda,meanIll,sDevIll,maxTime)
%This function simulates a continuous time stochastic epidemic model
%using the next-event update method.
%S0 and I0 are the initial populations of susceptible and infected
%individuals and N is the size of the total population
%lambda is the rate at which the infection is communicated to susceptibles
%we will choose a normal distribution with mean meanill and std. deviation
%sDevIll as the recovery time distribution.

I=I0;
S=S0;
N=S+I;
%Assign infection threshold for each susceptible individual which are
%exponentially distributed random variables with mean 1
phi=exprnd(1,[1,S0]);
%order phi
phi=sort(phi);

%Assign simulated Recovery Times to each susceptible and infected individual
T = meanIll + sDevIll.*randn(1,I0+S0);
  
%Simulation of the model
count=1;
time=0;
%Save initial state
timeMat(count,1)=time;
susceptMat(count,1)=S;
infectMat(count,1)=I;
recovMat(count,1)=N-S-I;
while I>0&&time<maxTime
    if I==0
        break
    end
    count=count+1;
    newTime=exprnd(1/lambda);
    b=min(I,numel(T));
    if newTime<min(T(1,1:b))
        %Calculate Infection Pressure
        sumPress=0;
        b=min(I,numel(T));
        %remove the elapsed time from the recovery times
        T(1,1:b)=T(1,1:b)-newTime*ones(1,b);
        for i=1:b
            sumPress=T(i)+sumPress;
        end
        infPress=(lambda*sumPress)/N;
        %Decide how many susceptibles get infected during the maximum infection
        %pressure
        sumPhi=0;
        i=1;
        while sumPhi<infPress&&phi(i)>0
            sumPhi=phi(i)+sumPhi;
            i=i+1;
            if i>numel(phi)
                break
            end
        end
    
        %Remove the newly infected person(s) from the infection threshold
        if i>S
            phi=0;
        else
            phi=phi(1,i:S);
        end
        %Increment I, S and time
        I=I+i-1;
        S=S-i+1;
        time=time+newTime;
        %Save the times that events happen, along with the population size at
        %those times
        timeMat(count,1)=time;
        susceptMat(count,1)=S;
        infectMat(count,1)=I;
        recovMat(count,1)=N-S-I;
    else
        %Let the soonest recovery happen, then restart the loop
        b=min(I,numel(T));
        time=time+min(T(1,b));
        m=numel(T);
        T(1,1:b)=T(1,1:b)-min(T(1,1:b))*ones(1,b);
      
        %Decrement I from the soonest recovery
        I=I-1;
    
        %Remove the infection from the infection times
        [~,j]=min(T(1,1:b));
        T=T([1:j-1,j+1:m]);
        %Save the times that events happen, along with the population size at
        %those times
        timeMat(count,1)=time;
        susceptMat(count,1)=S;
        infectMat(count,1)=I;
        recovMat(count,1)=N-S-I;
    end
    
end

%plot SIR classes over time
plot(timeMat,susceptMat,timeMat,infectMat,timeMat,recovMat)
legend('Susceptible','Infected', 'Recovered')
xlabel('Time')
ylabel('Population')
title('Stochastic Epidemic Model')

    
