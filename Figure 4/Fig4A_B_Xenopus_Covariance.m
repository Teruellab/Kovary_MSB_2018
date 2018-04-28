clear all
close all
[num geneNames] = xlsread('headers.xlsx');
A = xlsread('XenopusNoCalib300.xlsx');
A=A(1:294,:);

N=[0 28 58 86 116 146 175 204 234 264 294];
N2=[176:176+26 205:205+26 235:235+26 265:265+26];
N1=[1:1+26 29:29+26 59:59+26 87:87+26 117:117+26 147:147+26];
N12=[N1 N2];

%Batch Normalization for total protein loading variation, 3 subsets with 9
%sequentially analyzed eaggs
AA(1:294,1:27)=0;
for k=1:27
    AA(:,k)=(A(:,k))/median(A(:,k));
end

for j=1:294
    A(j,1:27)=A(j,1:27)/median(AA(j,[5:27])); %uses all peptides for normalization except Cyclin A1, Cyclin B2, Cdc6, and Emi1
end

% variation analysis by batching 9 consecutive eggs in each SRM run
% generates 30 separate variation results for each protein from the 300
% initial eggs
for i=1:10
    for j=1:3
        for k=1:27
            cr(3*(i-1)+j,k)=std(A(N(i)+9*(j-1)+1:N(i)+9*(j),k))/mean(A(N(i)+9*(j-1)+1:N(i)+9*(j),k));
            AV(3*(i-1)+j,k)=mean(A(N(i)+9*(j-1)+1:N(i)+9*(j),k));
            for h=1:27
                [cr2(3*(i-1)+j,k,h) er2(3*(i-1)+j,k,h)]=corr(A(N(i)+9*(j-1)+1:N(i)+9*(j),h),A(N(i)+9*(j-1)+1:N(i)+9*(j),k));
                if er2(3*(i-1)+j,k,h)>.1
                    cr2(3*(i-1)+j,k,h)=0;
                end
            end
        end
    end
end
for k=1:27
    for g=1:27
        [crx(k,g) erx(k,g)]=corr(A(N2,g),A(N2,k));
        if erx(k,g)>.05
            crx(k,g)=0;
        end
    end
end


%%
k=15;g=16;   %MEK and ERK
figure,plot(A(N2,k)/mean(A(N2,k)),A(N2,g)/mean(A(N2,g)),'ko')
data=fitlm(A(N2,k)/mean(A(N2,k)),A(N2,g)/mean(A(N2,g)));
xx=data.Coefficients.Estimate;
yy=data.Coefficients.SE;
hold on,line([0.8 1.2],[xx(2)+.8*xx(1) xx(2)+1.2*xx(1)])
axis([0.8 1.2 0.8 1.2])
[a b]=corr(A(N2,k),A(N2,g));
title(['Correlation ERK/MEK: ' num2str(a) '  p-value correlation vs no -correlation:  ' num2str(b)])

k=18;g=19;   %MCM5 and MCM7
figure,plot(A(N2,k)/mean(A(N2,k)),A(N2,g)/mean(A(N2,g)),'ko')
data=fitlm(A(N2,k)/mean(A(N2,k)),A(N2,g)/mean(A(N2,g)));
xx=data.Coefficients.Estimate;
yy=data.Coefficients.SE;
hold on,line([0.6 1.6],[xx(2)+.8*xx(1) xx(2)+1.2*xx(1)])
axis([0.6 1.6 0.6 1.4])
[a b]=corr(A(N2,k)/mean(A(N2,k)),A(N2,g)/mean(A(N2,g)));
title(['Correlation MCM5/MCM7: ' num2str(a) '  p-value, correlation vs no correlation:  ' num2str(b)])
%%
clear B
B(1:27,1:27)=0;
B=crx*2;
B(B>1)=1;
B(B<-1)=-1;
figure
for i=1:27
    hold on, rectangle('Position',[i-0.5 i-0.5 1 1],'FaceColor',[0 0 0])
    for j=i+1:27
        if isnan(B(i,j))==0
            if B(i,j)>0
                hold on, rectangle('Position',[i-0.5 j-0.5 1 1],'FaceColor',[B(i,j) 0 0])
                hold on, rectangle('Position',[j-0.5 i-0.5 1 1],'FaceColor',[B(i,j) 0 0])
            else
                hold on, rectangle('Position',[j-0.5 i-0.5 1 1],'FaceColor',[0 0 -B(i,j)])
                hold on, rectangle('Position',[i-0.5 j-0.5 1 1],'FaceColor',[0 0 -B(i,j)])
            end
        end
    end
end
title('Co-variance analysis, combined 60-80 minute time points, 2*6*9 eggs')
xlabel('27 measured proteins')
ylabel('27 measured proteins')
axis([.5 29.5 .5 27.5])
% hold on, figure
rectangle('Position',[28 0 1 10])
for i=1:500
    line([28 29],[5+i*.010 5+i*.010],'Color',[i/500 0 0])
    line([28 29],[5-i*.010 5-i*.010],'Color',[0 0 i/500])
end


