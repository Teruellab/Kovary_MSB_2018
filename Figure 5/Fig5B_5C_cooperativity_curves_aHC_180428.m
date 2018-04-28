clear all
close all
%Same model circuit as in Figure 1A, but now a threshold has been added
%at the last step so that cells end up either active (above the threshold)
%or inactive (below the threshold)
%Plot showing simulated fraction of cells activated at for different input
%stimuli
%Hill coefficients are fit and plotted to the data 
% 5%, 1.3; 10% 2.6; 20% 5.5; 40% 10.5
tspan = [0:15]; % time is in hours
y0 = [1 1 1 1 1]; % initial values

figure,hold on
for k=1:5
    if k==1
        errv=.01 %minimal noise case as reference added
    elseif k==2
        errv=.05; %5% variation of pathway components
    elseif k==3
        errv=.1; %10% variation of pathway components
    elseif k==4
        errv=.2; %20% variation of pathway components
    elseif k==5
        errv=.4; %40% variation of pathway components
    end
    
    for jj=1:100
        R =.5*2^((jj-1)*.1); %
        for j=1:1000
            e(1:10)=exp(randn(10,1)*errv); %Uncorrelated variation
            [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
            VV(j)=y(end,5);
        end
        
        VV0(jj,k)=mean(VV(:)>10);
        if k==1
            plot(log2(R),mean(VV(:)>10),'ko')
        elseif k==2
            plot(log2(R),mean(VV(:)>10),'ro')
        elseif k==3
            plot(log2(R),mean(VV(:)>10),'bo')     
        elseif k==4
            plot(log2(R),mean(VV(:)>10),'go')
        elseif k==5
            plot(log2(R),mean(VV(:)>10),'mo')
        end
    end
end
for i=1:24
    plot(log2(10),i*.04,'ko')
end

xlabel('Log2(R)')
ylabel('Percentage of cells that switch to active state at increasing input strength')
axis([-.3 4.3 -0.05 1.05])

for jj=1:100
        RR(jj) =.5*2^((jj-1)*.1); %      
end

ft=fittype('x^n/(x^n+10^n)');
for i=2:5
    cfit=fit(RR',squeeze(VV0(:,i)),ft);
    cfit
    n0(i)=cfit.n;
    if i==2
        hold on,plot(log2(RR),RR.^n0(i)./(RR.^n0(i)+10^n0(i)),'r')
    elseif i==3
        hold on,plot(log2(RR),RR.^n0(i)./(RR.^n0(i)+10^n0(i)),'b')
    elseif i==4
        hold on,plot(log2(RR),RR.^n0(i)./(RR.^n0(i)+10^n0(i)),'g')
    elseif i==5
        hold on,plot(log2(RR),RR.^n0(i)./(RR.^n0(i)+10^n0(i)),'m') 
    end
end
    %
    %%
    figure; bar(n0(2:5))  %barplot of apparent Hill coefficients
    %% 
function dydt = vd(t,y,R,e)
dydt = zeros(5,1);    % a column vector
dydt(1)=R*e(1)-e(2)*y(1);
dydt(2)=e(3)*y(1)-e(4)*y(2);
dydt(3)=e(5)*y(2)-e(6)*y(3);
dydt(4)=e(7)*y(3)-e(8)*y(4);
dydt(5)=e(9)*y(4)-e(10)*y(5);
end
