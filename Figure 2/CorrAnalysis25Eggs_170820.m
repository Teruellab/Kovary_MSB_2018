%% 
clear all
close all
cut=3;
[num geneNames] = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/headers1.xlsx');
A = xlsread('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/SRM stuff/Xenopus25_input_170819.xlsx');
% A=A(1:294,:); 
% A=eggdata;
% A=A(1:25,[4 6 10 7 1 2 3 8:9 11:28]);

figure, hold on
c1=[92 155 211];
c2=[255 124 47];
c3=[250 190 22];
c4=[165 165 165];
c5=75 133 183];
for j=1:5
    for k=1:27
        cr(j,k)=std(A(5*(j-1)+1:5*(j),k))/mean(A(5*(j-1)+1:5*(j),k));
        Av=A(5*(j-1)+1:5*(j),k)/mean(A(:,k));
        Av0=mean(Av);
        if j==1
            rectangle('Position',[k-.4+(j-1)*.15 0 .15 Av0],'FaceColor',c1./255)
        elseif j==2
            rectangle('Position',[k-.4+(j-1)*.15 0 .15 Av0],'FaceColor',c2./255)
        elseif j==3
            rectangle('Position',[k-.4+(j-1)*.15 0 .15 Av0],'FaceColor',c3./255)
        elseif j==4
            rectangle('Position',[k-.4+(j-1)*.15 0 .15 Av0],'FaceColor',c4./255)
        else
            rectangle('Position',[k-.4+(j-1)*.15 0 .15 Av0],'FaceColor',c5./255)
        end
        for i=1:5
        plot(k-.325+(j-1)*.15,Av(i),'ko','MarkerSize',4,'MarkerFaceColor','k')
        end
    end
end
axis([.5 27.5 0 2])
for i=1:27
    text(i,-.02,geneNames(i),'Rotation',270)
end
%% 
figure,plot(A(1:5,1),A(1:5,2),'r.')
hold on, plot(A(6:10,1),A(6:10,2),'b.')
plot(A(11:15,1),A(11:15,2),'g.')
plot(A(16:20,1),A(16:20,2),'c.')
plot(A(21:25,1),A(21:25,2),'k.')
plot(mean(A(1:5,1)),mean(A(1:5,2)),'ro')
text(mean(A(1:5,1)),mean(A(1:5,2)),' 0 min')
plot(mean(A(6:10,1)),mean(A(6:10,2)),'bo')
text(mean(A(6:10,1)),mean(A(6:10,2)),' 20 min')
plot(mean(A(11:15,1)),mean(A(11:15,2)),'go')
text(mean(A(11:15,1)),mean(A(11:15,2)),' 40 min')
plot(mean(A(16:20,1)),mean(A(16:20,2)),'co')
text(mean(A(16:20,1)),mean(A(16:20,2)),' 60 min')
plot(mean(A(21:25,1)),mean(A(21:25,2)),'ko')
text(mean(A(21:25,1)),mean(A(21:25,2)),' 80 min')
line([mean(A(1:5,1)) mean(A(6:10,1))],[mean(A(1:5,2)) mean(A(6:10,2))],'Color','r')
line([mean(A(6:10,1)) mean(A(11:15,1))],[mean(A(6:10,2)) mean(A(11:15,2))],'Color','b')
line([mean(A(11:15,1)) mean(A(16:20,1))],[mean(A(11:15,2)) mean(A(16:20,2))],'Color','g')
line([mean(A(16:20,1)) mean(A(21:25,1))],[mean(A(16:20,2)) mean(A(21:25,2))],'Color','c')
xlabel('Cyclin A')
ylabel('Cyclin B')
title('Time course of Cyclin A versus Cyclin B change')
figure
hold on
for i=1:27
    plot(i,100*cr(:,i),'k.')
    plot(i,100*median(cr(:,i)),'ro')
end
xlabel('27 proteins analyzed')
ylabel('% variation')
title('Variation analysis of 5 eggs at 5 time points')

[p x]=hist(median(cr),15);
figure,bar(100*x,p)
xlabel('Protein Variation in %')
ylabel('Number of times measured, 27 proteins 25 eggs analyzed for variation')
median(median(cr))
title('Frequency of measured variation')
cr25=cr;
save 'cr25' cr25

  
   