clear all
close all
%Same model Function as in Figure 1A
%Exploration of analog and binary signaling accuracy

%Plot 1 calculates Population-level Binary accuracy using an apparent Hill
%coefficient

%fIDL and aHC
tspan = [0:15]; % time is in hours
y0 = [1 1 1 1 1]; % initial values
Thr=10; % Threshold used for binary signaling, cells with y5 above Thr are active and cells with y5 below
% Thr are inactive
R1=.5:.01:12; %fold-increase in stimulation

CV=100*log(R1)/norminv(Thr/(1+Thr))/sqrt(10); %calculates how large the CV must be for the 
%tail region to be approximately 9% of cells (1-1/(Thr +1)) for a given R, The value (1-1/(Thr +1)) 
%corresponds to the fraction of cells that a Hill coefficient of 1 reaches when R=Thr

n=find(CV>5,1);
figure, plot(CV,log(Thr)./log(R1),'k-')
hold on
line([0;5],[log(Thr)/log(R1(n));log(Thr)/log(R1(n))],'Color', 'r','LineStyle','-')
plot(5,log(Thr)/log(R1(n)),'ro')
line([5 5],[0 log(Thr)/log(R1(n))],'Color', 'r','LineStyle','--')

n=find(CV>10,1);
hold on
line([0;10],[log(Thr)/log(R1(n));log(Thr)/log(R1(n))],'Color', 'b','LineStyle','-')
plot(10,log(Thr)/log(R1(n)),'bo')
line([10 10],[0 log(Thr)/log(R1(n))],'Color', 'b','LineStyle','--')

n=find(CV>20,1);
hold on
line([0;20],[log(Thr)/log(R1(n));log(Thr)/log(R1(n))],'Color', 'g','LineStyle','-')
plot(20,log(Thr)/log(R1(n)),'go')
line([20 20],[0 log(Thr)/log(R1(n))],'Color', 'g','LineStyle','--')

n=find(CV>40,1);
hold on
line([0;40],[log(Thr)/log(R1(n));log(Thr)/log(R1(n))],'Color', 'm','LineStyle','-')
plot(40,log(Thr)/log(R1(n)),'mo')

line([40 40],[0 log(Thr)/log(R1(n))],'Color', 'm','LineStyle','--')
title('Accuracy of binary population responses (aHC), High cooperativity reduces controllability by input')
xlabel('% Expression Variation')
ylabel('apparent Hill coefficient')
axis([0 43 1 11.5])

%% Plot calculates analog signaling accuracy (fold-Input Detection Limit) as a function of % expression variation
figure, hold on

x=0:.5:45; %Expression Variation
for i=1:length(x)
    y(i)=exp(norminv(.95)*x(i)/100*sqrt(10)*2); %determines fIDL based on 
    %two lognormal distributions, unstimulated and stimulated, 95%
    %accuracy & 10 variable pathway commponents
end
plot(x,log2(y),'k-')

y0=exp(norminv(.95)*5/100*sqrt(10)*2);
line([0 5],[log2(y0) log2(y0)],'LineStyle','-','Color','k')
plot(5,log2(y0),'ko')
line([5 5],[0 log2(y0)],'LineStyle','--','Color','r')

y0=exp(norminv(.95)*10/100*sqrt(10)*2);
line([0 10],[log2(y0) log2(y0)],'LineStyle','-','Color','b')
plot(10,log2(y0),'bo')
line([10 10],[0 log2(y0)],'LineStyle','--','Color','b')

y0=exp(norminv(.95)*20/100*sqrt(10)*2);
line([0 20],[log2(y0) log2(y0)],'LineStyle','-','Color','g')
plot(20,log2(y0),'go')
line([20 20],[0 log2(y0)],'LineStyle','--','Color','g')

y0=exp(norminv(.95)*40/100*sqrt(10)*2);
line([0 40],[log2(y0) log2(y0)],'LineStyle','-','Color','m')
plot(40,log2(y0),'mo')
line([40 40],[0 log2(y0)],'LineStyle','--','Color','m')

title('Analog single-cell signaling accuracy: fold-Input Detection Limits')
xlabel('% Expression Variation')
ylabel('fold-Input Detection Limit (fIDL) (log2)')
axis([0 43 0 6.3])
%% Covariance plot, codependence between fIDL and aHC

figure,hold on
CV=100*log(R1)/norminv(Thr/(1+Thr))/sqrt(10);
n1=find(CV>10,1);
YY=exp(norminv(.95)*CV/100*sqrt(10)*2);
yy1=exp(norminv(.95)*CV(n1)/100*sqrt(10)*2);
plot(log(Thr)./log(R1), log2(YY),'k-')
plot(log(Thr)/log(R1(n1)),log2(yy1),'ko')

CV=100*log(R1)/norminv(Thr/(1+Thr))/sqrt(50);
%sqrt(5^2 + 5^2) for correlated kinases and correlated phosphatases
n2=find(CV>10,1);
YY=exp(norminv(.95)*CV/100*sqrt(50)*2);
yy2=exp(norminv(.95)*CV(n2)/100*sqrt(50)*2);
plot(log(Thr)/log(R1(n2)),log2(yy2),'mo')

title('Covariance analysis (10% CV), Control (black), positive and negative regulators are co-varied (mangenta)')
ylabel('fold-INput Detection Limit (fIDL)(log2)')
xlabel('apparent Hill coefficient (aHC)')
axis([1 11.5 0 6.3])
%% Position of different expression variation in a fIDL versus aHC plot (N=10 components)
figure, hold on
kk=2;
CV=100*log(R1)/norminv(Thr/(1+Thr))/sqrt(10);
n1=find(CV>5,1);
n2=find(CV>10,1);
n3=find(CV>20,1);
n4=find(CV>40,1);
YY=exp(norminv(.95)*CV/100*sqrt(10)*2);
yy1=exp(norminv(.95)*CV(n1)/100*sqrt(10)*2);
yy2=exp(norminv(.95)*CV(n2)/100*sqrt(10)*2);
yy3=exp(norminv(.95)*CV(n3)/100*sqrt(10)*2);
yy4=exp(norminv(.95)*CV(n4)/100*sqrt(10)*2);

plot(log(Thr)./log(R1), log2(YY),'k-')
plot(log(Thr)/log(R1(n1)),log2(yy1),'ro')
plot(log(Thr)/log(R1(n2)),log2(yy2),'bo')
plot(log(Thr)/log(R1(n3)),log2(yy3),'go')
plot(log(Thr)/log(R1(n4)),log2(yy4),'mo')

title('Analog vs Binary signaling accurcay at 5(red), 10(blue) 20(green) and 40% (mangenta) expression variation')
ylabel('fold-INput Detection Limit (fIDL)(log2)')
xlabel('apparent Hill coefficient (aHC)')
axis([1 11.5 0 6.3])
%% Effect of less or more pathway components on fIDL versus aHC plot 

figure,hold on
Num=[4 10 24];
for kk=1:3
    CV=100*log(R1)/norminv(Thr/(1+Thr))/sqrt(Num(kk));
    n2=find(CV>10,1);    
    YY=exp(norminv(.95)*CV/100*sqrt(Num(kk))*2);
    yy2=exp(norminv(.95)*CV(n2)/100*sqrt(Num(kk))*2);
    if kk==1
        plot(log(Thr)/log(R1(n2)),log2(yy2),'co')  
    elseif kk==2
        plot(log(Thr)./log(R1), log2(YY),'k-')
        plot(log(Thr)/log(R1(n2)),log2(yy2),'ko')       
    elseif kk==3
        plot(log(Thr)/log(R1(n2)),log2(yy2),'mo')
    end
end

title('Accurcay for different Numbers of Components (N) (Var 10%): Cyan, N=4, Black 10 and Mangenta 24')
ylabel('fold-INput Detection Limit (fIDL)(log2)')
xlabel('apparent Hill coefficient (aHC)')
axis([1 11.5 0 6.3])
%% 


function dydt = vd(t,y,R,e)
dydt = zeros(5,1);    % a column vector
dydt(1)=R*e(1)-e(2)*y(1);
dydt(2)=e(3)*y(1)-e(4)*y(2);
dydt(3)=e(5)*y(2)-e(6)*y(3);
dydt(4)=e(7)*y(3)-e(8)*y(4);
dydt(5)=e(9)*y(4)-e(10)*y(5);
end
