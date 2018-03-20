clear all
close all
cc=2; %case cc=2 for no co-variation, cc=1 for co-variation
%interval used for histogram plots
for k=1:3
    if k==1
        errv=.05; %error assumption of the model
        int=.1;
    elseif k==2
        errv=.1;
        int=.2;
    elseif k==3
        errv=.25;
        int=.4;
    end
    tspan = [0:15]; % time is in hours
    y0 = [1 1 1 1 1]; % initial values
    figure,hold on
    for j=1:5000
        for i=[1 2 3]
            R =3^(i-1); %
            if cc==1
                e([1 3 5 7 9])=exp(randn(1,1)*errv); % correlated variation in positive and negative regulators
                e([2 4 6 8 10])=exp(randn(1,1)*errv);
            else
                e(1:10)=exp(randn(10,1)*errv); %Uncorrelated equally strong variation
            end
            [t,y] = ode45(@(t,y) vd(t,y,R,e), tspan,y0);
            VV(j,i)=y(end,5);
            %     plot(t,log(y(:,1))/log(2),'r-')
            if j<10
                if i==1
                    plot(t,log(y(:,5))/log(2),'k-')
                elseif i==2
                    plot(t,log(y(:,5))/log(2),'r-')
                elseif i==3
                    plot(t,log(y(:,5))/log(2),'b-')
                elseif i==4
                    plot(t,log(y(:,5))/log(2),'c-')
                end
            end
        end
    end
    xlabel('Time')
    ylabel('log2, Relative Output intensity')
    title(['Relative stimulus R of 1, 3 and 9-fold, Noise: ' num2str(errv)])
    axis([0 15 -2.5 5.5])
    clear P
    for i=[1 2 3]
        [p x]=hist(log2(VV(:,i)),-3:int:7);
        P(1:length(x),i)=p;
    end
    figure,plot(x,P(:,1),'k-')
    hold on, plot(x,P(:,2),'r-')
    plot(x,P(:,3),'b-')
    %plot(x,P(:,4),'c-')
    xlabel('log2, Relative Output intensity')
    ylabel('Number of cells')
    title(['Relative stimulus R of 1, 3 and 9-fold, Noise: ' num2str(errv)])
    axis([-2.5 5.5 0 max(P(:))*1.1])
end
function dydt = vd(t,y,R,e)
dydt = zeros(5,1);    % a column vector
dydt(1)=R*e(1)-e(2)*y(1);
dydt(2)=e(3)*y(1)-e(4)*y(2);
dydt(3)=e(5)*y(2)-e(6)*y(3);
dydt(4)=e(7)*y(3)-e(8)*y(4);
dydt(5)=e(9)*y(4)-e(10)*y(5);
end
