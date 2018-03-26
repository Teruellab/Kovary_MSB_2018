clear all
close all
%Same model circuit as in Figure 1A.
%Plot showing the importance of co-variation in broadening the stimulus window
%over which receptor inputs can regulate the fraction of cells that switch
%to the on state
%The simulation assumes that cells switch to the on-state when the output
%y(5) becomes greater than 8
%paramter cc=2 for no co-variation, cc=1 for co-variation in the same
%signaling circuit
tspan = [0:15]; % time is in hours
y0 = [1 1 1 1 1]; % initial values
R1=0:.01:log2(12); %over an 8-fold increase in R, the fraction of cells that switch should increase 8-fold
%for optimal linear population-level sensitivy; 50% at
%R==8, 6.25% at R==1
% Rmax=log(8);
% % 
% % figure,hold on
% CV=100*R1/sqrt(10)/erfinv(.875)/sqrt(2);
% n1=find(CV>10,1);
% YY=exp(erfinv(.9)*sqrt(10)*sqrt(2)*CV/100*2);
% yy1=exp(erfinv(.9)*sqrt(10)*sqrt(2)*CV(n1)/100*2);
% x=-log2(R1/Rmax);
% y=log2(YY);
% fill([5 x(2:end) 5],[6.5 y(2:end) 6.5],[.9 .9 .9])
% plot(-log2(R1/Rmax), log2(YY),'k-')
% plot(-log2(R1(n1)/Rmax),log2(yy1),'ko')
% 
% CV=100*R1/5/sqrt(2)/erfinv(.875)/sqrt(2);
% n2=find(CV>10,1);
% YY=exp(erfinv(.9)*5*sqrt(2)*sqrt(2)*CV/100*2);
% yy2=exp(erfinv(.9)*5*sqrt(2)*sqrt(2)*CV(n2)/100*2);
% plot(-log2(R1(n2)/Rmax),log2(yy2),'mo')
% title('Competing accuracy limits: 5 step signaling at 10% var with (red) and without (black) covariance ')
% ylabel('Analog accuracy')
% xlabel('log2 Binary Population accuracy: cooperativity reduces input sensitivity')
% 
% % axis([0 5.5 0 6.5])

figure,hold on
for k=1:3
    if k==1
        errv=.01 %minimal noise case as reference added
        cc=2;
    elseif k==2
        errv=.2; %random variation of pathway components
        cc=2;
    elseif k==3
        errv=.2; %co-variant noise in pathway
        cc=1;
    end
 
    for jj=1:50
        R =.5*2^((jj-1)*.1); %
        RR(jj)=R;
        for j=1:1000
            if cc==1
                e([1 3 5 7 9])=exp(randn(1,1)*errv); % correlated variation in positive and negative regulators
                e([2 4 6 8 10])=exp(randn(5,1)*errv); %uncorrelated
            else
                e(1:10)=exp(randn(10,1)*errv); %Uncorreated equally strong variation
            end
            [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
            VV(j)=y(end,5);
        end
        
        VV0(jj,k)=mean(VV(:)>10);
        if k==1
            plot(log2(R),mean(VV(:)>10),'ko')
        elseif k==2
            plot(log2(R),mean(VV(:)>10),'bo')
        elseif k==3
            plot(log2(R),mean(VV(:)>10),'mo')
        end
    end
end
for i=1:24
            plot(log2(10),i*.04,'ko')
end

xlabel('Log2(R)')
ylabel('fraction of cells that switch cell fate')
axis([-.5 4.3 -0.05 1.05])
% sel=RR<10;
% RR=RR(sel);
%VV0=VV0(sel,:);
ft=fittype('x^n/(x^n+10^n)');
for i=2:3
    cfit=fit(RR',squeeze(VV0(:,i)),ft);
    cfit
    n0=cfit.n;
    if i==2
    hold on,plot(log2(RR),RR.^n0./(RR.^n0+10^n0),'b')
    elseif i==3
            hold on,plot(log2(RR),RR.^n0./(RR.^n0+10^n0),'m')
    end
        
end

% % 
% clear all
% 
% %signaling circuit
% tspan = [0:15]; % time is in hours
% y0 = [1 1 1 1 1]; % initial values
% R1=0:.01:log(10); %over an 8-fold increase in R, the fraction of cells that switch should increase 8-fold
% cc=2;
% 
% figure,hold on
% for k=1:5
%     if k==1      %Defines a range of variations from very small to 60%
%         errv=.001 %lognormal variation parameter
%     elseif k==2
%         errv=.05
%     elseif k==3
%         errv=.10
%     elseif k==4
%         errv=.2
%     elseif k==5
%         errv=.4
%     end
%     for jj=1:60
%         R =.5*2^((jj-1)*.1); %Defines a range of simulated increasing receptor stimuli from 0.5 on higher
%         RR(jj)=R;
%         for j=1:1000
%             
%             e(1:10)=exp(randn(10,1)*errv);
%             [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
%             VV(j)=y(end,5);
%         end
%         if k==1
%             plot(log2(R),mean(VV(:)>8),'ko')
%         elseif k==2
%             plot(log2(R),mean(VV(:)>8),'ro')
%         elseif k==3
%             plot(log2(R),mean(VV(:)>8),'bo')
%         elseif k==4
%             plot(log2(R),mean(VV(:)>8),'Marker','o','Color','g')
%         elseif k==5
%             plot(log2(R),mean(VV(:)>8),'Marker','o','Color','m')
%         end
%         
%         for j=1:1000
%             if cc==1
%                 e([1 3 5 7 9])=exp(randn(1,1)*errv); % correlated variation in positive and negative regulators
%                 e([2 4 6 8 10])=exp(randn(1,1)*errv);
%             else
%                 e(1:10)=exp(randn(10,1)*errv); %Uncorreated equally strong variation
%             end
%             [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
%             VV(j)=y(end,5);
%         end
%         
%         VV0(jj,k)=mean(VV(:)>8);
%     end
% end
% for i=1:24
%     plot(3,i*.04,'ko')
% end
% title('Fracional differentiation as a function of stimulus, 0,5,10,20,40% variation')
% xlabel('Log2(R)')
% ylabel('fraction of cells that switch cell fate')
% ft=fittype('x^n/(x^n+8^n)');
% for i=2:5
%     cfit=fit(RR',squeeze(VV0(:,i)),ft);
%     cfit
%     n0=cfit.n;
%     if i==2
%     hold on,plot(log2(RR),RR.^n0./(RR.^n0+8^n0),'r')
%     elseif i==3
%             hold on,plot(log2(RR),RR.^n0./(RR.^n0+8^n0),'b')
%             elseif i==4
%             hold on,plot(log2(RR),RR.^n0./(RR.^n0+8^n0),'g')
%             elseif i==5
%             hold on,plot(log2(RR),RR.^n0./(RR.^n0+8^n0),'m')
%     end
%         
% end
% 
% axis([-.5 3.5 -.05 1.05])
% 


function dydt = vd(t,y,R,e)
dydt = zeros(5,1);    % a column vector
dydt(1)=R*e(1)-e(2)*y(1);
dydt(2)=e(3)*y(1)-e(4)*y(2);
dydt(3)=e(5)*y(2)-e(6)*y(3);
dydt(4)=e(7)*y(3)-e(8)*y(4);
dydt(5)=e(9)*y(4)-e(10)*y(5);
end
