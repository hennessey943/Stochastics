function iterateStochasticEpidemicModel(S0,I0,lambda,meanIll,sDevIll,maxTime)
for i=1:1000
    [timeMat,susceptMat,infectMat,recovMat]=stochasticEpidemicModel(S0,I0,lambda,meanIll,sDevIll,maxTime);
    
    n=numel(infectMat);
    m=numel(recovMat);
    everInfected(1,i)=recovMat(m)+infectMat(n);
end
histogram(everInfected)