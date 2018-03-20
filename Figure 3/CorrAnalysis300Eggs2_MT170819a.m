clear all
close all
[num geneNames] = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/headers1.xlsx');
A = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/XenopusNoCalib300_MT_reordered.xlsx');
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
        A(j,1:27)=A(j,1:27)/median(AA(j,[5:27])); %uses all peptides for normalization except Cyclin A1 and B2, Emi1
end
%figure,hist(mean(AA(:,4:27),2))
%analysis of 9 consecutive eggs in SRM run for variation analysis,
%generates 30 separate variation results for each protein from the 300
%initial eggs 
 
figure,hold on
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
    if AV(j,k)>0
        plot(log(mean(AV(19:24,k)))/log(10),100*mean(cr(19:24,k)),'ok')
        text(log(mean(AV(19:24,k)))/log(10)+.1,100*mean(cr(19:24,k)),char(geneNames(k)))
    end
end
ylabel('% variation')
xlabel('relative protein abundance, log10 units') 
title('Increased apparent variation for cell cycle regulated and low abundant proteins')
axis([-2.5 2 0 24])
 
cr300=cr(19:24,:);
load 'cr25' cr25
 
figure,plot(100*mean(cr300(:,:),1),100*mean(cr25(:,:),1),'ko')
text(100*mean(cr300(:,:),1)+.2,100*mean(cr25(:,:),1),geneNames)
axis([0 30 0 30])
xlabel('%Variation 60 min, 6 x 9 egg batches, 300 egg set')
ylabel('%Variation 5 x 5 egg batches, 25 egg set')
title('Correlated variation in independent data sets')
 
for i=1:10
    for j=1:3
        for k=1:27
        A(N(i)+9*(j-1)+1:N(i)+9*(j),k)=A(N(i)+9*(j-1)+1:N(i)+9*(j),k)/mean(A(N(i)+9*(j-1)+1:N(i)+9*(j),k));       
        end
    end
end
 
%Directly testing for ERK/MEK and MCM5/7 co-variance in combined 60-80
%minute datasets
k=13;g=14;
figure,plot(A(N2,k),A(N2,g),'ko')
data=fitlm(A(N2,k),A(N2,g));
xx=data.Coefficients.Estimate;
yy=data.Coefficients.SE;
hold on,line([0.8 1.2],[xx(2)+.8*xx(1) xx(2)+1.2*xx(1)])
axis([0.8 1.2 0.8 1.2])
[a b]=corr(A(N2,k),A(N2,g));
title(['Correlation ERK/MEK: ' num2str(a) '  p-value correlation vs no -correlation:  ' num2str(b)])
 
k=24;g=25;
figure,plot(A(N2,k),A(N2,g),'ko')
data=fitlm(A(N2,k),A(N2,g));
xx=data.Coefficients.Estimate;
yy=data.Coefficients.SE;
hold on,line([0.8 1.2],[xx(2)+.8*xx(1) xx(2)+1.2*xx(1)])
axis([0.8 1.2 0.8 1.2])
[a b]=corr(A(N2,k),A(N2,g));
title(['Correlation MCM5/MCM7: ' num2str(a) '  p-value, correlation vs no correlation:  ' num2str(b)])
 
%Plotting variance analysis for 60 and 80 minute time points
figure
for j=3:4
    int=6*j+1:6*j+6;
    hold on
    for i=1:27
        if j==3
            plot(i-.15,100*cr(int,i),'b.')
            plot(i-.15,100*median(cr(int,i)),'bo')
        elseif j==4
            plot(i+.15,100*cr(int,i),'r.')
            plot(i+.15,100*median(cr(int,i)),'ro')
        end
    end
end
title(['Variance analysis: Blue 60 min and red 80 min time points 6 batches of 9 eggs each'])
hold on
for i=1:27 
    x25=prctile(100*cr(19:30,i),25);
    x75=prctile(100*cr(19:30,i),75);
    rectangle('Position',[i-.25 x25 .5 x75-x25])
end
axis([0 28 0 35]) 
ylabel('% variance')
 
cv=cr(19:24,:);
[p x]=hist(cv(:),.0:.02:.36);
figure,bar(100*x ,p)
title('Variance analysis of 27 proteins at 60 minute time point, 6 repeats, batches of 9')
xlabel('% Variance')
ylabel('frequency of variance observed')
 
clear B
B(1:27,1:27)=0;
for i=1:18
    for h=1:27
        for j=1:27
            if er2(i,h,j)<.05
                if cr2(i,h,j)>0
                    B(h,j)=B(h,j)+1;
                elseif cr2(i,h,j)<0
                    B(h,j)=B(h,j)-1;
                end
            end
        end
    end
end
B=B/4.5;
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
title('Co-variance analysis, combined 0-20-40 minute time points, 3*6*9 eggs')
xlabel('27 measured proteins, 13,14 is ERK-MEK, 24,25 is MCM5/7')
ylabel('27 measured proteins')
hold on, rectangle('Position',[28 0 1 10])
for i=1:500
    line([28 29],[5+i*.010 5+i*.010],'Color',[i/500 0 0])
    line([28 29],[5-i*.010 5-i*.010],'Color',[0 0 i/500])
end
axis([.5 29.5 .5 27.5])
 
BB1=mean(cr2(1:18,:,:),1);
 
clear B
B(1:27,1:27)=0;
for i=19:30
    for h=1:27
        for j=1:27
            if er2(i,h,j)<.05
                if cr2(i,h,j)>0
                    B(h,j)=B(h,j)+1;
                elseif cr2(i,h,j)<0
                    B(h,j)=B(h,j)-1;
                end
            end
        end
    end
end
BB2=mean(cr2(19:30,:,:),1);
B=B/3;
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
axis([.5 29.5 .5 27.5])
hold on, rectangle('Position',[28 0 1 10])
for i=1:500
    line([28 29],[5+i*.010 5+i*.010],'Color',[i/500 0 0])
    line([28 29],[5-i*.010 5-i*.010],'Color',[0 0 i/500])
end
 
figure,scatter(BB1(:),BB2(:),'o')
hold on,plot(BB1(1,13,14),BB2(1,13,14),'ro')
hold on,plot(BB1(1,24,25),BB2(1,24,25),'go')
axis([-.5 .5 -.5 .5])
title(['repeatability of co-variance, combined 60-80 vs 0-20-40 timepoints, ERK/MEK pair red, MCM5/7 pair green'])
