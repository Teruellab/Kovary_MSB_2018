clear all
close all
%Same model function as in Figure 1A
% estimation of signal R needed to create 95% accuracy in signal
% transmission; accuracy of 95% accuracy=.05; var 7 or 15 %, .07 or .15;
% NumberSignalingComponents is 10 in the Pathway in Figure1A
% Factor of 2 is because the assumption is for same variation in basal
% and stimulated state

%CDF=1/2*(1+erf(x/CV/sqrt(2))  cummulative distribution function and erf
%erfinv(2*CDF-1)=x/CV/sqrt(2)  to integrate a one-sided 5% tail one has to
%take erfinv(.9), the CV that gives this tail is then relative
%1/erfinv(.9)/sqrt(2), about 0.6
%if there are ten random variables in the system and the tail should be at
%an R0 2-fold stimulus,
%the corresponding R0=exp(sqrt(10)*erfinv(.9)*sqrt(2)*CV)
% x=3:.5:45;
NumberSignalingComponents=10;
y0=exp(norminv(.95)*.05*sqrt(10)*sqrt(2)*2);
yy(1)=y0;

y0=exp(norminv(.95)*.1*sqrt(10)*sqrt(2)*2);
yy(2)=y0;

y0=exp(norminv(.95)*.15*sqrt(10)*sqrt(2)*2);
yy(3)=y0;

y0=exp(norminv(.95)*.25*sqrt(10)*sqrt(2)*2);
yy(4)=y0;

figure,bar(log2(yy))
title('Rel. Detection limits for 5, 10, 15, and 25% variation')

R0=yy(2);

int=.5; %interval used for histogram plots
tspan = [0:15]; % time is in hours
y0 = [1 1 1 1 1]; % initial values
errv=0.15;
for j=1:15000
    for i=1:2
        if i==1
            R=1;
        elseif i==2
            R=R0;
        end
        e(1:10)=exp(randn(10,1)*errv); %Uncorreated equally strong variation
        [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
        VV(j,i)=y(end,5);
    end
end

for i=1:2
    [p x]=hist(log2(VV(:,i)),-2.5:.1:5);
    P(1:length(x),i)=p;
end
figure,plot(x,P(:,1),'b-')
hold on, plot(x,P(:,2),'r-')
axis([-2.5 5 0 max(P(:))*1.1])
line([log2(R)/2 log2(R)/2],[0 2000],'LineStyle','--')

title(['Histogram for R=1 vs R=2.83; Detection limit at 10%'])
sum(log2(VV(:,1))>log2(yy(2))/2)/ length(VV(:,1))
sum(log2(VV(:,2))<log2(yy(2))/2)/ length(VV(:,2))

function dydt = vd(t,y,R,e)
dydt = zeros(5,1);    % a column vector
dydt(1)=R*e(1)-e(2)*y(1);
dydt(2)=e(3)*y(1)-e(4)*y(2);
dydt(3)=e(5)*y(2)-e(6)*y(3);
dydt(4)=e(7)*y(3)-e(8)*y(4);
dydt(5)=e(9)*y(4)-e(10)*y(5);
end
