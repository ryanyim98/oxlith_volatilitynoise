clear; clc;
addpath('~/Desktop/BayesM10/');
data_path = '~/Desktop/BayesM10/data/';

%flag: run Bayesian Filter on all data, only pre-treatment or only
%post-treatment? The manuscript used all data.
all = 1;
pre_treatment = 0;
post_treatment = 0;

%this is the file where the patients data were stored in long format
file_name = 'movement_data.csv'; 
movement_data = readtable([data_path file_name]);

%record demographic data separately
demo = unique(table(movement_data.ID,movement_data.Gender,movement_data.Age,movement_data.Allocation,'VariableNames',{'ID','Gender','Age','Allocation'}));

%% 

%which days to keep
if all
    daysKept = [-14 50]; %only keep past one week
elseif pre_treatment
    daysKept = [-49  0];
else
    daysKept = [1  50];
end

% make new struct to store values
for participant=1:size(demo,1)
        alldata.subjectdata(participant).ID = demo.ID(participant);
        alldata.subjectdata(participant).Gender = demo.Gender(participant);
        alldata.subjectdata(participant).Age = demo.Age(participant);
        alldata.subjectdata(participant).Allocation = demo.Allocation(participant);
        alldata.subjectdata(participant).allData = movement_data(string(cell2mat(movement_data.ID)) == string(alldata.subjectdata(participant).ID),[5:13 24:25]); %these are the data columns
        
%     alldata.subjectdata(participant).preTreatmentData = alldata.subjectdata(participant).allData(alldata.subjectdata(participant).allData.Day <= 0,:);
%     alldata.subjectdata(participant).postTreatmentData = alldata.subjectdata(participant).allData(alldata.subjectdata(participant).allData.Day > 0,:);

end

% make data that aligns at day 0 (different patients started at different
% dates pre-randomisation)
for participant=1:size(demo,1)
    alldata.subjectdata(participant).data_aligned.day = nan(100,1);
    alldata.subjectdata(participant).data_aligned.L5.mean = nan(100,1);
    alldata.subjectdata(participant).data_aligned.M10.mean = nan(100,1);
    alldata.subjectdata(participant).data_aligned.Time5.mean = nan(100,1);
    alldata.subjectdata(participant).data_aligned.Time10.mean = nan(100,1);
    alldata.subjectdata(participant).data_aligned.pos.mean = nan(100,1);
    alldata.subjectdata(participant).data_aligned.neg.mean = nan(100,1);
    
    days_num = length(alldata.subjectdata(participant).allData.Day);
    
    %realign the dates
    day0_loc = find(alldata.subjectdata(participant).allData.Day == 0);
    days_loc = (50-day0_loc+1):(50+days_num-day0_loc);
    
    for d = 1:days_num
        alldata.subjectdata(participant).data_aligned.day(days_loc(d)) = alldata.subjectdata(participant).allData.Day(d);
        alldata.subjectdata(participant).data_aligned.L5.mean(days_loc(d)) = alldata.subjectdata(participant).allData.L5(d);
        alldata.subjectdata(participant).data_aligned.M10.mean(days_loc(d)) = alldata.subjectdata(participant).allData.M10(d);
        alldata.subjectdata(participant).data_aligned.Time5.mean(days_loc(d)) = alldata.subjectdata(participant).allData.Time5(d);
        alldata.subjectdata(participant).data_aligned.Time10.mean(days_loc(d)) = alldata.subjectdata(participant).allData.Time10(d);
        alldata.subjectdata(participant).data_aligned.pos.mean(days_loc(d)) = alldata.subjectdata(participant).allData.pos(d);
        alldata.subjectdata(participant).data_aligned.neg.mean(days_loc(d)) = alldata.subjectdata(participant).allData.neg(d);
    end 
end

%% 'OL0002' has weird data
% participant = 1;
% alldata.subjectdata(participant).allData.M10(alldata.subjectdata(participant).allData.Day <= 0,:) = nan;
% alldata.subjectdata(participant).allData.L5(alldata.subjectdata(participant).allData.Day <= 0,:) = nan;
% alldata.subjectdata(participant).data_aligned.M10.mean(1:50) = nan;
% alldata.subjectdata(participant).data_aligned.L5.mean(1:50) = nan;


%% transform L5 and M10 into 0~1 bounded scale
%concatenate all L5 and M10 vals
L5 = []; M10 = [];
for participant=1:size(demo,1)
    L5 = [L5, alldata.subjectdata(participant).data_aligned.L5.mean];
    M10 = [M10, alldata.subjectdata(participant).data_aligned.M10.mean];
end

%% plot the movement and onset time values
L5 = L5((daysKept(1)+50):(daysKept(2)+50),:);
M10 = M10((daysKept(1)+50):(daysKept(2)+50),:);

max_L5 = max(L5(:));
min_L5 = min(L5(:));
max_M10 = max(M10(:));
min_M10 = min(M10(:));

subplot(1,2,1);
hist(L5(:)); title(["L5, max = "+max_L5+", min = "+min_L5]);
subplot(1,2,2);
hist(M10(:));title(["M10, max = "+max_M10+", min = "+min_M10]);

%exclude outliers
outlierL5 = quantile(L5(:),0.99)
outlierM10 = quantile(M10(:),0.99)


% transform L5 and M10 into 0~1 scale
for participant=1:size(demo,1)
    alldata.subjectdata(participant).allData = alldata.subjectdata(participant).allData((alldata.subjectdata(participant).allData.Day >= daysKept(1) & alldata.subjectdata(participant).allData.Day <= daysKept(2)),:);
    alldata.subjectdata(participant).allData.L5 = (alldata.subjectdata(participant).allData.L5 - min_L5)/(max_L5-min_L5)*0.6+0.2;
    alldata.subjectdata(participant).allData.M10 = (alldata.subjectdata(participant).allData.M10 - min_M10)/(max_M10-min_M10)*0.6+0.2;
end

%% concatenate all L5 and M10 vals
Time5 = []; Time10 = [];
for participant=1:size(demo,1)
    Time5 = [Time5, alldata.subjectdata(participant).data_aligned.Time5.mean];
    Time10 = [Time10, alldata.subjectdata(participant).data_aligned.Time10.mean];
end


Time5 = Time5((daysKept(1)+50):(daysKept(2)+50),:);
Time10 = Time10((daysKept(1)+50):(daysKept(2)+50),:);

max_Time5 = max(Time5(:));
min_Time5 = min(Time5(:));
max_Time10 = max(Time10(:));
min_Time10 = min(Time10(:));

subplot(1,2,1);
hist(Time5(:)); title(["Time5, max = "+max_Time5+", min = "+min_Time5]);
subplot(1,2,2);
hist(Time10(:));title(["Time10, max = "+max_M10+", min = "+min_Time10]);

%exclude outliers
outlierTime5 = quantile(Time5(:),0.99)
outlierTime10 = quantile(Time10(:),0.99)

% transform T5 and T10 into 0~1 scale
for participant=1:size(demo,1)
    alldata.subjectdata(participant).allData = alldata.subjectdata(participant).allData((alldata.subjectdata(participant).allData.Day >= daysKept(1) & alldata.subjectdata(participant).allData.Day <= daysKept(2)),:);
    alldata.subjectdata(participant).allData.Time5 = (alldata.subjectdata(participant).allData.Time5 - min_Time5)/(max_Time5-min_Time5)*0.6+0.2;%*0.8+0.1
    alldata.subjectdata(participant).allData.Time10 = (alldata.subjectdata(participant).allData.Time10 - min_Time10)/(max_Time10-min_Time10)*0.6+0.2;
end

%% transform pos and neg into 0~1 scale
%concatenate all L5 and M10 vals
pos = []; neg = [];
for participant=1:size(demo,1)
    pos = [pos, alldata.subjectdata(participant).data_aligned.pos.mean];
    neg = [neg, alldata.subjectdata(participant).data_aligned.neg.mean];
end

pos = pos((daysKept(1)+50):(daysKept(2)+50),:);
neg = neg((daysKept(1)+50):(daysKept(2)+50),:);

max_pos = max(pos(:));
min_pos = min(pos(:));
max_neg = max(neg(:));
min_neg = min(neg(:));

subplot(1,2,1);
hist(pos(:)); title(["pos, max = "+max_pos+", min = "+min_pos]);
subplot(1,2,2);
hist(neg(:));title(["neg, max = "+max_neg+", min = "+min_neg]);

%exclude outliers
outlierPos = quantile(pos(:),0.99)
outlierNeg = quantile(neg(:),0.99)


% transform PA and NA into 0~1 scale

for participant=1:size(demo,1)
    alldata.subjectdata(participant).allData = alldata.subjectdata(participant).allData((alldata.subjectdata(participant).allData.Day >= daysKept(1) & alldata.subjectdata(participant).allData.Day <= daysKept(2)),:);
    alldata.subjectdata(participant).allData.pos = (alldata.subjectdata(participant).allData.pos - min_pos)/(max_pos-min_pos)*0.6+0.2;%*0.6+0.2;
    alldata.subjectdata(participant).allData.neg = (alldata.subjectdata(participant).allData.neg - min_neg)/(max_neg-min_neg)*0.6+0.2;
end

%% save data
save([data_path 'dataraw.mat'],'alldata');
