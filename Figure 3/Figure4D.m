clear all
close all
for kkk=1:2
    if kkk==1
        casev=2;celltype=1;  %1 for Hela, 2 for MCF10a
    elseif kkk==2
        casev=2;celltype=2;  %1 for Hela, 2 for MCF10a
    elseif kkk==3
        casev=1;celltype=1; %1 for Hela, 2 for MCF10a
    elseif kkk==4
        casev=1;celltype=2;  %1 for Hela, 2 for MCF10a
    end
    if casev==2
        if celltype==1
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170413 - MCF10A Hela Correlations MEK ERK GAPDH ENO MCM5 MCM7/Tracked3/Hela/AllData.mat') %4 DAPI 3
            %Protein;  1 and 2 are antibodies
        elseif celltype==2
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170413 - MCF10A Hela Correlations MEK ERK GAPDH ENO MCM5 MCM7/Tracked3/MCF10A/AllData.mat')
        end
    elseif casev==1
        if celltype==1
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170429-MEK-ERK-Cor-MCF10A-Hela/Hela/AllData.mat') %1DAPI, 4Protein stain, 2 and 3 antobodies
            
        elseif celltype==2
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170429-MEK-ERK-Cor-MCF10A-Hela/MCF10A/AllData.mat') %1DAPI, 4Protein stain, 2 and 3 antobodies
        end
    elseif casev==3
        if celltype==1
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170415 - MCF10A Hela - Many Correlations/Hela/AllData.mat') %1DAPI, 4Protein stain, 2 and 3 antobodies
        elseif celltype==2
            load('/Users/mary/Documents/MATLAB/Kyle_Xenopus_Paper/imaging/170415 - MCF10A Hela - Many Correlations/MCF10A/AllData.mat') %1DAPI, 4Protein stain, 2 and 3 antobodies
        end
    end
    
    % Summarize data contents
    condition_names = fieldnames(AllData);
    measurement_names = fieldnames(AllData.(condition_names{1}).Measurements);
    disp('Conditions in this dataset:')
    disp(condition_names)
    disp('Measurements made (for cells in each condition):')
    disp(measurement_names)
    disp('- - - - - - - - ')
    
    %%  REORGANIZE DATA to make it easier to pull out & compare selected conditions  - - - - - - - -
    
    % Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
    if casev==1
        [xdata_by_condition, xdata_by_well] = restructuredata(AllData,'IntegratedNuc1'); % 1st measurment we want to compare
        [ydata1_by_condition, ydata1_by_well] = restructuredata(AllData,'IntegratedCell2'); % A 2nd measurement we want to compare
        [ydata2_by_condition, ydata2_by_well] = restructuredata(AllData,'IntegratedCell3');
        [zdata_by_condition, zdata_by_well] = restructuredata(AllData,'IntegratedCell4'); % A 4th measurement we want to compare
        [xcoor_by_condition, xcoor_by_well] = restructuredata(AllData,'CentroidX'); % A 4th measurement we want to compare
        [ycoor_by_condition, ycoor_by_well] = restructuredata(AllData,'CentroidY'); % A 4th measurement we want to compare
    elseif casev==2
        [xdata_by_condition, xdata_by_well] = restructuredata(AllData,'IntegratedNuc4'); % 1st measurment we want to compare
        [ydata1_by_condition, ydata1_by_well] = restructuredata(AllData,'IntegratedCell1'); % A 2nd measurement we want to compare
        [ydata2_by_condition, ydata2_by_well] = restructuredata(AllData,'IntegratedCell2');
        [zdata_by_condition, zdata_by_well] = restructuredata(AllData,'IntegratedCell3'); % A 4th measurement we want to compare
        [xcoor_by_condition, xcoor_by_well] = restructuredata(AllData,'CentroidX'); % A 4th measurement we want to compare
        [ycoor_by_condition, ycoor_by_well] = restructuredata(AllData,'CentroidY'); % A 4th measurement we want to compare
    elseif casev==3
        [xdata_by_condition, xdata_by_well] = restructuredata(AllData,'IntegratedNuc4'); % 1st measurment we want to compare
        [ydata1_by_condition, ydata1_by_well] = restructuredata(AllData,'IntegratedCell1'); % A 2nd measurement we want to compare
        [ydata2_by_condition, ydata2_by_well] = restructuredata(AllData,'IntegratedCell2');
        [zdata_by_condition, zdata_by_well] = restructuredata(AllData,'IntegratedCell3'); % A 4th measurement we want to compare
        [xcoor_by_condition, xcoor_by_well] = restructuredata(AllData,'CentroidX'); % A 4th measurement we want to compare
        [ycoor_by_condition, ycoor_by_well] = restructuredata(AllData,'CentroidY'); % A 4th measurement we want to compare
    end
    
    if casev==1
        if celltype==1
            bound=[160000 280000]; %Hela
        elseif celltype==2
            %bound=[90000 150000]; %MCF10a
            bound=[160000 250000]; %MCF10a
        end
        nlist=1:9;
    elseif casev==2
        if celltype==1
            bound=[47000 60000]; %Hela
        elseif celltype==2
            bound=[27000 33000]; %MCF10a
        end
        nlist=1:9;
        %     elseif casev==3
        %         if celltype==1
        %             bound=[310000 380000]; %Hela
        %         elseif celltype==2
        %             bound=[170000 210000]; %MCF10a
        %         end
        %         nlist=1:9;
    end
    vtot=[];wtot=[];zvtot=[];zwtot=[];
    for n=nlist   % number of conditions
        for m=1:3  % 3 wells per condition
            DAPI1=xdata_by_well{n,1}{m,1};
            POI1=ydata1_by_well{n,1}{m,1};
            POI2=ydata2_by_well{n,1}{m,1};
            PS=zdata_by_well{n,1}{m,1};  %protein stain
            xcoor=xcoor_by_well{n,1}{m,1};
            ycoor=ycoor_by_well{n,1}{m,1};
            if n==1 & m==1
                figure,hist(DAPI1,200)   % plot DAPI histogram for first condition, first well
            end
            sel=DAPI1>bound(1) & DAPI1<bound(2) & ((xcoor-540).^2 + (ycoor-540).^2)<350^2;
            DAPI1=DAPI1(sel);
            POI1=POI1(sel);
            POI2=POI2(sel);
            PS=PS(sel);  %protein stain
            xcoor=xcoor(sel);
            ycoor=ycoor(sel);
            POI2=POI2./PS;
            POI1=POI1./PS;
            ztot1=POI1;
            ztot2=POI2;
            v2=POI2(:)/mean(POI2(:));
            v1=POI1(:)/mean(POI1(:));
            if n==1  % for first condition
                if m==1   % first well for that condition
                    v1add1=v1;
                    v2add1=v2;
                else
                    v1add1=[v1add1;v1];
                    v2add1=[v2add1;v2];
                end
                if m==3
                    figure,plot(v1add1,v2add1,'.')
                    if length(v1add1)>10
                        data=fitlm(v1add1,v2add1,'Robustopts','on');
                        xx=data.Coefficients.Estimate;
                        yy=data.Coefficients.SE;
                        hold on,line([0 2],[xx(1) xx(1)+2*xx(2)])
                        [a b]=corr(v1add1,v2add1);
                        title(['Correlation MCM5/7: ' num2str(a) '  p-value, correlation vs no correlation:  ' num2str(b)])
                        axis([0 2 0 2])
                    end
                end
            end
            if n==2
                if m==1
                    v1add2=v1;
                    v2add2=v2;
                else
                    v1add2=[v1add2;v1];
                    v2add2=[v2add2;v2];
                end
                if m==3
                    figure,plot(v1add2,v2add2,'.')
                    if length(v1add2)>10
                        data=fitlm(v1add2,v2add2);
                        xx=data.Coefficients.Estimate;
                        yy=data.Coefficients.SE;
                        [a b]=corr(v1add2,v2add2);
                        title(['Correlation MCM5/GAPDH: ' num2str(a) '  p-value, correlation vs no correlation:  ' num2str(b)])
                        axis([0 2 0 2])
                    end
                end
            end
            if n==5
                if m==1
                    v1add3=v1;
                    v2add3=v2;
                else
                    v1add3=[v1add3;v1];
                    v2add3=[v2add3;v2];
                end
                if m==3
                    if length(v1add3)>10
                        figure,plot(v1add3,v2add3,'.')
                        data=fitlm(v1add3,v2add3,'Robustopts','on');
                        xx=data.Coefficients.Estimate;
                        yy=data.Coefficients.SE;
                        hold on,line([0 2],[xx(1) xx(1)+2*xx(2)])
                        [a b]=corr(v1add3,v2add3);
                        title(['Correlation ERK/MEK: ' num2str(a) '  p-value, correlation vs no correlation:  ' num2str(b)])
                        axis([0 2 0 2])
                    end
                end
            end
            if n==4 | n==5
                vtot=[vtot;v2];
                zvtot=[zvtot;ztot2];
            end
            if n==5 | n==6
                wtot=[wtot;v1];
                zwtot=[zwtot;ztot1];
            end
            
            sel=v2<.5 | v2>2 | v1<.5 | v1>2; %cut outliers larger than factor 2
            v1(sel)=[];
            v2(sel)=[];
            if length(v1)>10
                [d f]=corr(v1,v2);
                v2=std(v2);
                v1=std(v1);%
            else
                d=0;f=1;v1=0;v2=0;
                [n m kkk]
            end
            R(n,m,1:2)=[v1 v2];
            RR(n,m)=d;
            RRerr(n,m)=f;
        end
    end
    xlabel('1 MCM7; 4,5 MEK 2,3,6-9 GAPDH')
    ylabel('% Variation')
    %figure,bar(100*R(:,:,2))
    xlabel('4,5,7 ERK; 6,9 MEK; 1,3 MCM5; 2 MCM7; 8 Eno')
    ylabel('% Variation')
    
    figure,bar(RR([1 2 5],:))
    if kkk==1
        save('File1.mat', 'R','RR','vtot' ,'wtot','zvtot','zwtot');
    elseif kkk==2
        save('File2.mat', 'R','RR','vtot' ,'wtot','zvtot','zwtot');
    elseif kkk==3
        save('File3.mat', 'R','RR','vtot' ,'wtot','zvtot','zwtot');
    elseif kkk==4
        save('File4.mat', 'R','RR','vtot' ,'wtot','zvtot','zwtot');
    end
end
function [by_condition, by_well, by_image] = restructuredata(AllData, measure_name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [by_condition, by_well, by_image] = structuredata(AllData, measure_name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% RESTRUCTUREDATA splits a single-cell measurement (e.g. 'MeanNuc1') into 3 different cell matricies, organized:
% 1) condition
% 2) WELLS within a given cfigure,bar(R(:,:,1))ondition
% 3) IMAGES within a given condition
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Collect data together
all_cond = fieldnames(AllData);
by_condition = cell(size(all_cond));
by_well = cell(size(all_cond));
by_image = cell(size(all_cond));

for i = 1:length(all_cond)
    if ~isfield(AllData.(all_cond{i}).Measurements,measure_name)
        error(['Error: Measurement field ''', measure_name,''' does not exist in one or more conditions.'])
    end
    by_condition{i} =  AllData.(all_cond{i}).Measurements.(measure_name)(:);
    
    % Group by well
    wells = unique(AllData.(all_cond{i}).CellData(:,1));
    by_well{i} = cell(size(wells));
    for j = 1:length(wells)
        by_well{i}{j} = AllData.(all_cond{i}).Measurements.(measure_name)(AllData.(all_cond{i}).CellData(:,1)==wells(j));
    end
    
    % Group by image
    images = unique(AllData.(all_cond{i}).CellData(:,2));
    by_image{i} = cell(size(images));
    for j = 1:length(images)
        by_image{i}{j} = AllData.(all_cond{i}).Measurements.(measure_name)(AllData.(all_cond{i}).CellData(:,2)==images(j));
    end
    
end
end
