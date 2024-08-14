clear; clc;
addpath('/Users/yanyan/Desktop/BayesM10/');
addpath('/Users/yanyan/Desktop/BayesM10/Oxlith_code');
data_path = '/Users/yanyan/Desktop/BayesM10/data/';

%load the matrix formatted by A1 script
load([data_path 'dataraw.mat']);


%%
%set the range of the parameters, otherwise the function
%maglearn_func_vardiff_flat_miss will set the param range

params.vmurange=[1e-7 1]; % the original range
% params.vmurange=[1e-10 1]; %low vmu range; did not make a substantial
% difference

% the data domain; choose among "movement" "time" "mood" 
to_learn = "time";  %"movement" "time" "mood"


dirname=['raw',num2str(datenum(datetime('now')))];

if to_learn == "movement"
    dataCol = [7 9];
elseif to_learn == "time"
    dataCol = [6 8];
elseif to_learn == "mood"
    dataCol = [10 11];
end

%running the Bayesian Filter
for participant = 1:size(alldata.subjectdata,2)
    for datatypeid= 1:2 %L5 vs. M10; time5 vs. time10; PA vs. NA
        datatype = dataCol(datatypeid);
        datause = table2array(alldata.subjectdata(participant).allData(:,datatype));
        part_class(participant,datatypeid).out=maglearn_func_vardiff_flat_miss(datause,params);%grid inputs from Mike
    end
end

for participant = 1:size(alldata.subjectdata,2)
    for datatypeid= 1:2 %L5 vs. M10; time5 vs. time10; PA vs. NA
        datatype = dataCol(datatypeid);
        alldata.subjectdata(participant).learner(datatypeid)=part_class(participant,datatypeid).out;
    end
end

save([data_path+to_learn+'_data_learned_all_-14toEnd_allparticipants.mat'],'alldata');
