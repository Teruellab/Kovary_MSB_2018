clear all
close all
[num geneNames] = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/headers1.xlsx');
A = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/XenopusNoCalib300_MT_reordered.xlsx');
A=A(1:294,:);

N=[0 28 58 86 116 146 175 204 234 264 294]; %Start of egg data sets, 0,0,20,20,40,40,60,60,80 and 80 min sets
%Select first 27 eggs per set to have same number in each set
N1=[1:1+26 29:29+26 59:59+26 87:87+26 117:117+26 147:147+26]; %0, 20 and 40 min sets
N2=[176:176+26 205:205+26 235:235+26 265:265+26]; %60 and 80 min sets 
N12=[N1 N2]; % all sets

%Correction of residual small calibration error in total mass per egg
AA(1:294,1:27)=0; %Parameter for normalized correction for each peptide
for k=1:27
    AA(:,k)=(A(:,k))/mean(A(:,k));
end
%23 peptides are used for normalization of each egg in the abundance matrix A
%(except for Cyclin A1, B2, Emi1 and Cdc6 that vary during cell cycle)

x=randn(300,27);
y=exp(.2*x);
c=mean(y,2);
figure,plot(sort(c))
figure,bar(std(y,1)./mean(y,1));

for i=1:300
    y(:,k)=y(:,k)/c(i);
end
figure,bar(std(y,1)./mean(y,1));

N0=204;
for j=1:6
    for k=1:27
        cx(j,k)=std(A(N0+10*(j-1)+1:N0+10*(j),k))/mean(A(N0+10*(j-1)+1:N0+10*(j),k));
    end
end

figure,hold on
for k=1:27
    for i=1:6
        plot(k,cx(i,k),'ok')
    end
    plot(k,mean(cx(:,k)),'or')
end

for j=1:294
        A(j,1:27)=A(j,1:27)/mean(AA(j,[5:27]));
        corrV(j)=mean(AA(j,[5:27]));
end

figure,hist(corrV,-.5:.01:1.5)
%% 

% analysis of batches of 9 consecutive eggs in SRM run for variation analysis,
% generates 30 separate variation results for each protein from the 300
% initial eggs 

figure,hold on
for i=1:10
    for j=1:3
        for k=1:27
            cr(3*(i-1)+j,k)=std(A(N(i)+9*(j-1)+1:N(i)+9*(j),k))/mean(A(N(i)+9*(j-1)+1:N(i)+9*(j),k));
            AV(3*(i-1)+j,k)=mean(A(N(i)+9*(j-1)+1:N(i)+9*(j),k));
        end
    end
end

%Bootstrapping of corrected data at 60 minutes, 60 eggs analyzed

B=A(204:263,1:27);

for i=1:27
    CV(i)=std(B(:,i))/mean(B(:,i));
end
figure,bar(CV)

N0=204;
for j=1:6
    for k=1:27
        cx(j,k)=std(A(N0+10*(j-1)+1:N0+10*(j),k))/mean(A(N0+10*(j-1)+1:N0+10*(j),k));
    end
end

figure,hold on
for k=1:27
    for i=1:6
        plot(k,cx(i,k),'ok')
    end
    plot(k,mean(cx(:,k)),'or')
end
figure,hold on
for k=1:27
    B0=B(:,k);
    a=bootstrp(30,@(x)[std(x) mean(x)],squeeze(B0));
    plot(k*ones(30),a(:,1)./a(:,2),'.');  
    plot(k,mean(a(:,1)./a(:,2)),'ro')
end

%Plots the average variation of each protein at the 60 minute time set  
%as y-axis and the relative abundance of the same protein on the y-axis, Figure 2D
for k=1:27
    if AV(j,k)>0
        plot(log(mean(AV(19:24,k)))/log(10),100*mean(cr(19:24,k)),'ok')
        text(log(mean(AV(19:24,k)))/log(10)+.1,100*mean(cr(19:24,k)),char(geneNames(k)))
    end
end
ylabel('% variation')
xlabel('relative protein abundance, log10 units') 
title('Increased apparent variation for cell cycle regulated and low abundant proteins')
axis([-2.5 2 0 24])

%To test for reproducibility, this plot compares the measured variation of the 54 eggs at the 60 min point with the
%variation measured for the 25 egg set
cr300=cr(19:24,:);
% load 'cr25' cr25
% figure,plot(100*mean(cr300(:,:),1),100*mean(cr25(:,:),1),'ko')
% text(100*mean(cr300(:,:),1)+.2,100*mean(cr25(:,:),1),geneNames)
% axis([0 30 0 30])
% xlabel('%Variation 60 min, 6 x 9 egg batches, 300 egg set')
% ylabel('%Variation 5 x 5 egg batches, 25 egg set')
% title('Correlated variation in independent data sets')

%Variation analysis for the 60 and 80 minute egg sets
figure
for j=3:4
    int=6*j+1:6*j+6; %variance from 9 eggs in each data point, 6 variations measured at 60 and 6 at 80 minutes    
    hold on
    for i=1:27
        if j==3
            plot(i-.15,100*cr(int,i),'b.')
            %plot(i-.15,100*mean(cr(int,i)),'bo') %marks the median variation with a circle
        elseif j==4
            plot(i+.15,100*cr(int,i),'r.')
            %plot(i+.15,100*mean(cr(int,i)),'ro')
        end
    end
end
title(['Variance analysis: Blue 60 min and red 80 min time points 6 batches of 9 eggs each'])
hold on
for i=1:27 %marks the 25-75% range with a box
    xa=100*std(cr(19:24,i))/sqrt(5);
    xb=100*std(cr(25:30,i))/sqrt(5);
    xma=mean(100*cr(19:24,i));
    xmb=mean(100*cr(25:30,i));
    

    rectangle('Position',[i-.25 xma-xa .2 2*xa])
    rectangle('Position',[i+.05 xmb-xb .2 2*xb])

end
axis([0 28 0 35]) 
ylabel('% variance')

[p x]=hist(100*cr(:),0:1.5:35);
figure,bar(x(1:end-1),p(1:end-1))
xlabel('Protein Variation in %')
ylabel('Number of times measured, 27 proteins 25 eggs analyzed for variation')
median(median(cr))
title('Frequency of masured variation')

