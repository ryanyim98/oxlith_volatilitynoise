clc; clear;
addpath('/Users/rh/Desktop/BayesM10/');
data_path = "/Users/rh/Desktop/BayesM10/data/";

pre_data = load([data_path+'data_learned_pre_treatment.mat']);
post_data = load([data_path+'data_learned_pre_treatment.mat']);

%% realign data
for participant = 1:size(alldata.subjectdata,2)
        %create an array of 100 days, of which No. 50 corresponds to the
        %randomisation day (day 0)
        alldata.subjectdata(participant).data_aligned.L5.sEst = nan(100,1);
        alldata.subjectdata(participant).data_aligned.L5.vmuEst = nan(100,1);
    
        alldata.subjectdata(participant).data_aligned.L5.muEst = nan(100,1);
        alldata.subjectdata(participant).data_aligned.M10.sEst = nan(100,1);
        alldata.subjectdata(participant).data_aligned.M10.vmuEst = nan(100,1);
        alldata.subjectdata(participant).data_aligned.M10.muEst = nan(100,1);
    
        days_num = length(alldata.subjectdata(participant).allData.Day);
        day0_loc = find(alldata.subjectdata(participant).allData.Day == 0);
        days_loc = (50-day0_loc+1):(50+days_num-day0_loc);
        for d = 1:days_num
            alldata.subjectdata(participant).data_aligned.day(days_loc(d)) = alldata.subjectdata(participant).allData.Day(d);
            alldata.subjectdata(participant).data_aligned.L5.sEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).sEst(d);
            alldata.subjectdata(participant).data_aligned.L5.muEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).muEst(d);
            alldata.subjectdata(participant).data_aligned.L5.vmuEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).vmuEst(d);
        
            alldata.subjectdata(participant).data_aligned.M10.sEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).sEst(d);
            alldata.subjectdata(participant).data_aligned.M10.vmuEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).vmuEst(d);
            alldata.subjectdata(participant).data_aligned.M10.muEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).muEst(d);

        end
end