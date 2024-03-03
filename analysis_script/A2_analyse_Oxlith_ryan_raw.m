clear; clc;
addpath('/Users/rh/Desktop/BayesM10/');
addpath('/Users/rh/Desktop/BayesM10/Oxlith_code');
data_path = '/Users/rh/Desktop/BayesM10/data/';
load([data_path 'dataraw.mat']);

params.vmurange=[1e-15 1];
% params.vmurange=[1e-7 1];
to_learn = "movement";  %"movement" "time" "mood"

dirname=['raw',num2str(datenum(datetime('now')))];

if to_learn == "movement"
    dataCol = [7 9];
elseif to_learn == "time"
    dataCol = [6 8];
elseif to_learn == "mood"
    dataCol = [10 11];
end

for participant = 1:size(alldata.subjectdata,2)
    participant
    for datatypeid= 1:2 %L5 vs. M10
        datatype = dataCol(datatypeid);
        datause = table2array(alldata.subjectdata(participant).allData(:,datatype));
        part_class(participant,datatypeid).out=maglearn_func_vardiff_flat_miss(datause,params);%grid inputs from Mike
    end
end

for participant = 1:size(alldata.subjectdata,2)
    for datatypeid= 1:2 %L5 vs. M10
        datatype = dataCol(datatypeid);
        alldata.subjectdata(participant).learner(datatypeid)=part_class(participant,datatypeid).out;
    end
end

save([data_path+to_learn+'_data_learned_all_-14toEnd_allparticipants_lowvmu.mat'],'alldata');
