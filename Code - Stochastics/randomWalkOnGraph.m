function randomWalkOnGraph(n,p)

%Monte-Carlo Simulation of Random Walk on Graph

P=[0,p,0,0,1-p;
    1-p,0,p,0,0;
    0,1-p,0,p,0;
    0,0,1-p,0,p;
    p,0,0,1-p,0];

pi=[1/5,1/5,1/5,1/5,1/5];

%Forward Simulation Expected Value

% for m=1:1000
%     
%     pos=zeros(n);
%     t=1;
%     while t<=n+1
%     
%     x=rand;
%     i=1;
% 
%     if t==1
%         while x>0
%             x=x-pi(i);
%             i=i+1;
%         end
%     else
%         while x>0
%             x=x-P(pos(t-1),i);
%             i=i+1;
%         end
%     end
%     
%     pos(t)=i-1;
%     t=t+1;
% 
%     end
%     for j=1:5
%         expf(j,m)=sum(pos(:)==j);
%     end
% end
% 
% avgf=mean(expf')
% 
% %Forward Simulation
% pos=zeros(n+1);
% t=1;
% while t<=n+1
% x=rand;
% i=1;
% 
%     if t==1
%     while x>0
%        x=x-pi(i);
%        i=i+1;
%     end
%     else
%         while x>0
%             x=x-P(pos(t-1),i);
%             i=i+1;
%         end
%     end
%     
% pos(t)=i-1;
% t=t+1;
% 
% end
% figure
% stairs(0:1:n,pos)
% title('Forward Random Walk on Graph');
% 
% xlabel('Epoch');
% ylabel('Position on Graph');
% 
% Reverse Prob-Trans Matrix
 Pr=P';
% 
% Reverse Simulation
% 
% t=n+1;
% while t>=1
%     x=rand;
%     j=1;
%     if t==n+1
%         while x>0
%             x=x-pi(j);
%             j=j+1;
%         end
%     else
%         while x>0
%             x=x-Pr(negpos(n-t+1),j);
%             j=j+1;
%         end
%     end
%     negpos(n+2-t)=j-1;
%     t=t-1;
% end
% figure
% stairs(n:-1:0,negpos)
% title('Reverse Random Walk on Graph');
% 
% xlabel('Epoch');
% ylabel('Position on Graph');
% 
% Reverse Simulation Expected Values
% 
% for m=1:5000
%     
%     pos=zeros(n);
%     t=n+1;
%     while t>=1
%         x=rand;
%         j=1;
%         if t==n+1
%             while x>0
%                 x=x-pi(j);
%                 j=j+1;
%             end
%         else
%             while x>0
%                 x=x-Pr(pos(t+1),j);
%                 j=j+1;
%             end
%         end
%         pos(t)=j-1;
%         t=t-1;
%     end
%     
%     for j=1:5
%         expr(j,m)=sum(pos(:)==j);
%     end
% end
% 
% avgr=mean(expr')
% 
% diff=norm(avgf-avgr)

%Backward to -n Simulation

t=n;
while t>=-100
    x=rand;
    j=1;
    if t==n
        while x>0
            x=x-pi(j);
            j=j+1;
        end
    else
        while x>0
            x=x-Pr(negpos(n+2+t),j);
            j=j+1;
        end
    end
    negpos(n+1+t)=j-1;
    t=t-1;
end
% t=0;
% while(t>=-n && t<=0)
%     x=rand;
%     j=1;
%     if t==0
%     while x>0
%         x=x-P(pos(t+1),j);
%         j=j+1;
%     end
%     else
%         while x>0
%             x=x-P(negpos(-t),j);
%             j=j+1;
%         end
%     end
%     negpos(-t+1)=j-1;
%     t=t-1;
% end
% 
% position=zeros(2*n+1);
% position(1:n+1)=negpos;
% position(n+2:2*n+1)=pos;
% 
figure
 stairs(-n:n,negpos)
title('Reverse Random Walk from 100 to -100')
xlabel('Epoch')
ylabel('Position')

%Expected Values



