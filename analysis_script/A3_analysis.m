clc; clear;
addpath('/Users/rh/Desktop/BayesM10/');
data_path = "/Users/rh/Desktop/BayesM10/data/";

data = "all"; %"post_treatment", "pre_treatment", "all"
if data ~= "post_treatment"
    realign = 1;
else
    realign = 0;
end

to_learn = "movement";  %"movement" "time"

load([data_path+to_learn+'_data_learned_all_-14toEnd_allparticipants_lowvmu.mat']);

%% realign data (day0 = last day before treatment)
if realign
    for participant = 1:size(alldata.subjectdata,2)
        %create an array of 100 days, of which No. 50 corresponds to the
        %randomisation day (day 0)
        if to_learn ~= "mood"
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
            
        elseif to_learn == "mood"
            alldata.subjectdata(participant).data_aligned.pos.sEst = nan(100,1);
            alldata.subjectdata(participant).data_aligned.pos.vmuEst = nan(100,1);
            alldata.subjectdata(participant).data_aligned.pos.muEst = nan(100,1);
        
            alldata.subjectdata(participant).data_aligned.neg.sEst = nan(100,1);
            alldata.subjectdata(participant).data_aligned.neg.vmuEst = nan(100,1);
            alldata.subjectdata(participant).data_aligned.neg.muEst = nan(100,1);
    
            days_num = length(alldata.subjectdata(participant).allData.Day);
            day0_loc = find(alldata.subjectdata(participant).allData.Day == 0);
            days_loc = (50-day0_loc+1):(50+days_num-day0_loc);
            for d = 1:days_num
                alldata.subjectdata(participant).data_aligned.day(days_loc(d)) = alldata.subjectdata(participant).allData.Day(d);
                alldata.subjectdata(participant).data_aligned.pos.sEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).sEst(d);
                alldata.subjectdata(participant).data_aligned.pos.muEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).muEst(d);
                alldata.subjectdata(participant).data_aligned.pos.vmuEst(days_loc(d)) = alldata.subjectdata(participant).learner(1).vmuEst(d);
        
                alldata.subjectdata(participant).data_aligned.neg.sEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).sEst(d);
                alldata.subjectdata(participant).data_aligned.neg.vmuEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).vmuEst(d);
                alldata.subjectdata(participant).data_aligned.neg.muEst(days_loc(d)) = alldata.subjectdata(participant).learner(2).muEst(d);
            end
        end
        
    end
elseif ~realign
    if to_learn ~= "mood"
        for participant = 1:size(alldata.subjectdata,2)
            alldata.subjectdata(participant).data_aligned.L5.sEst = alldata.subjectdata(participant).learner(1).sEst;
            alldata.subjectdata(participant).data_aligned.L5.vmuEst = alldata.subjectdata(participant).learner(1).vmuEst;
            alldata.subjectdata(participant).data_aligned.L5.muEst = alldata.subjectdata(participant).learner(1).muEst;
            alldata.subjectdata(participant).data_aligned.M10.sEst = alldata.subjectdata(participant).learner(2).sEst;
            alldata.subjectdata(participant).data_aligned.M10.vmuEst = alldata.subjectdata(participant).learner(2).vmuEst;
            alldata.subjectdata(participant).data_aligned.M10.muEst = alldata.subjectdata(participant).learner(2).muEst;
        end
    elseif to_learn == "mood"
            alldata.subjectdata(participant).data_aligned.pos.sEst = alldata.subjectdata(participant).learner(1).sEst;
            alldata.subjectdata(participant).data_aligned.pos.vmuEst = alldata.subjectdata(participant).learner(1).vmuEst;
            alldata.subjectdata(participant).data_aligned.pos.muEst = alldata.subjectdata(participant).learner(1).muEst;
            alldata.subjectdata(participant).data_aligned.neg.sEst = alldata.subjectdata(participant).learner(2).sEst;
            alldata.subjectdata(participant).data_aligned.neg.vmuEst = alldata.subjectdata(participant).learner(2).vmuEst;
            alldata.subjectdata(participant).data_aligned.neg.muEst = alldata.subjectdata(participant).learner(2).muEst;
    end
end

%% export L5/M10 to csv for HLM
participant_data = {};
all_data_out = {};
for participant = 1:size(alldata.subjectdata,2)
    participant_data(:,1:4) = repmat([alldata.subjectdata(participant).ID alldata.subjectdata(participant).Gender alldata.subjectdata(participant).Age alldata.subjectdata(participant).Allocation],[100 1]);
    participant_data(:,5) = num2cell(-49:50);
    if to_learn == "movement"
        participant_data(:,6:9) = [num2cell(alldata.subjectdata(participant).data_aligned.L5.mean) num2cell(alldata.subjectdata(participant).data_aligned.L5.sEst) num2cell(alldata.subjectdata(participant).data_aligned.L5.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.L5.muEst)];
        participant_data(:,10:13) = [num2cell(alldata.subjectdata(participant).data_aligned.M10.mean) num2cell(alldata.subjectdata(participant).data_aligned.M10.sEst) num2cell(alldata.subjectdata(participant).data_aligned.M10.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.M10.muEst)];
    elseif to_learn == "time"
        participant_data(:,6:9) = [num2cell(alldata.subjectdata(participant).data_aligned.Time5.mean) num2cell(alldata.subjectdata(participant).data_aligned.L5.sEst) num2cell(alldata.subjectdata(participant).data_aligned.L5.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.L5.muEst)];
        participant_data(:,10:13) = [num2cell(alldata.subjectdata(participant).data_aligned.Time10.mean) num2cell(alldata.subjectdata(participant).data_aligned.M10.sEst) num2cell(alldata.subjectdata(participant).data_aligned.M10.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.M10.muEst)];
    elseif to_learn == "mood"
        participant_data(:,6:9) = [num2cell(alldata.subjectdata(participant).data_aligned.pos.mean) num2cell(alldata.subjectdata(participant).data_aligned.pos.sEst) num2cell(alldata.subjectdata(participant).data_aligned.pos.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.pos.muEst)];
        participant_data(:,10:13) = [num2cell(alldata.subjectdata(participant).data_aligned.neg.mean) num2cell(alldata.subjectdata(participant).data_aligned.neg.sEst) num2cell(alldata.subjectdata(participant).data_aligned.neg.vmuEst) num2cell(alldata.subjectdata(participant).data_aligned.neg.muEst)];
    end
    all_data_out = [all_data_out;participant_data];
end

if to_learn == "mood"
     all_data_out = cell2table(all_data_out,"VariableNames",["ID" "Gender" "Age" "Allocation" "day" "pos_mean" "pos_s" "pos_vmu" "pos_mu" "neg_mean" "neg_s" "neg_vmu" "neg_mu"]);
else
    all_data_out = cell2table(all_data_out,"VariableNames",["ID" "Gender" "Age" "Allocation" "day" "L5_mean" "L5_s" "L5_vmu" "L5_mu" "M10_mean" "M10_s" "M10_vmu" "M10_mu"]);
end

writetable(all_data_out,[data_path+to_learn+"_data_HLM_-14toEnd_allparticipants_lowvmu.csv"]);

%% organise into one matrix for further analysis
if to_learn ~= "mood"
    L5 = nan(4,2,34,100); %param, group, participant, time
    M10 = nan(4,2,34,100); %param, group, participant, time

    for participant = 1:size(alldata.subjectdata,2)
        dayLength = length(alldata.subjectdata(participant).data_aligned.L5.sEst);
        if char(alldata.subjectdata(participant).Allocation) == 'A'
            groupid = 1;
        else
            groupid = 2;
        end
    
        L5(1,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.L5.sEst;
        L5(2,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.L5.vmuEst;
        L5(3,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.L5.muEst;

        M10(1,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.M10.sEst;
        M10(2,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.M10.vmuEst;
        M10(3,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.M10.muEst;
    
        if to_learn == "movement"
            if realign
                L5(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.L5.mean;
                M10(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.M10.mean;
                    %day 0 missing
                L5(4,groupid,participant,50) =  (L5(4,groupid,participant,49)+L5(4,groupid,participant,51))/2;
                M10(4,groupid,participant,50) =  (M10(4,groupid,participant,49)+M10(4,groupid,participant,51))/2;
            elseif ~realign
                L5(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.L5;
                M10(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.M10;
            end
        elseif to_learn == "time"
                if realign
                    L5(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.Time5.mean;
                    M10(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.Time10.mean;
                    %day 0 missing
                    L5(4,groupid,participant,50) =  (L5(4,groupid,participant,49)+L5(4,groupid,participant,51))/2;
                    M10(4,groupid,participant,50) =  (M10(4,groupid,participant,49)+M10(4,groupid,participant,51))/2;
                elseif ~realign
                    L5(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.Time5;
                    M10(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.Time10;
                end
        end

    end
elseif to_learn == "mood"
    pos = nan(4,2,34,100); %param, group, participant, time
    neg = nan(4,2,34,100); %param, group, participant, time

    for participant = 1:size(alldata.subjectdata,2)
        dayLength = length(alldata.subjectdata(participant).data_aligned.pos.sEst);
        if char(alldata.subjectdata(participant).Allocation) == 'A'
            groupid = 1;
        else
            groupid = 2;
        end
    
        pos(1,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.pos.sEst;
        pos(2,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.pos.vmuEst;
        pos(3,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.pos.muEst;

        neg(1,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.neg.sEst;
        neg(2,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.neg.vmuEst;
        neg(3,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.neg.muEst;
    
        
        if realign
                    pos(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.pos.mean;
                    neg(4,groupid,participant,(1:dayLength)) = alldata.subjectdata(participant).data_aligned.neg.mean;
                    %day 0 missing
                    pos(4,groupid,participant,50) =  (pos(4,groupid,participant,49)+pos(4,groupid,participant,51))/2;
                    neg(4,groupid,participant,50) =  (neg(4,groupid,participant,49)+neg(4,groupid,participant,51))/2;
                elseif ~realign
                    pos(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.pos;
                    neg(4,groupid,participant,(1:(dayLength-1))) = alldata.subjectdata(participant).allData.neg;
        end

    end
end

L5_sEst = squeeze(mean(L5(1,:,:,:),3,'omitnan'));
L5_vmuEst = squeeze(mean(L5(2,:,:,:),3,'omitnan'));
L5_muEst = squeeze(mean(L5(3,:,:,:),3,'omitnan'));
L5_mean = squeeze(mean(L5(4,:,:,:),3,'omitnan'));
M10_sEst = squeeze(mean(M10(1,:,:,:),3,'omitnan'));
M10_vmuEst = squeeze(mean(M10(2,:,:,:),3,'omitnan'));
M10_muEst = squeeze(mean(M10(3,:,:,:),3,'omitnan'));
M10_mean = squeeze(mean(M10(4,:,:,:),3,'omitnan'));

L5_ci = nan(4,2,100,2);
M10_ci = nan(4,2,100,2);
for i = 1:100
    %only operate on non-NaN values
    if ~isnan(L5_sEst(1,i))
        %upper
        L5_ci(:,:,i,1) = mean(L5(:,:,:,i),3,'omitnan') - 1.96*std(L5(:,:,:,i),0,3,'omitnan')./sqrt(sum(~isnan(L5(:,:,:,i)),3));
        M10_ci(:,:,i,1) = mean(M10(:,:,:,i),3,'omitnan') - 1.96*std(M10(:,:,:,i),0,3,'omitnan')./sqrt(sum(~isnan(M10(:,:,:,i)),3));
        %lower
        L5_ci(:,:,i,2) = mean(L5(:,:,:,i),3,'omitnan') + 1.96*std(L5(:,:,:,i),0,3,'omitnan')./sqrt(sum(~isnan(L5(:,:,:,i)),3));
        M10_ci(:,:,i,2) = mean(M10(:,:,:,i),3,'omitnan') + 1.96*std(M10(:,:,:,i),0,3,'omitnan')./sqrt(sum(~isnan(M10(:,:,:,i)),3));
    end
end


%%
if realign
    axis = (min(find(~isnan(L5_ci(1,1,:,1))))-50): (max(find(~isnan(L5_ci(1,1,:,1))))-50);
    axis_id = axis+50;
elseif ~realign
    axis = min(find(~isnan(L5_ci(1,1,:,1)))): max(find(~isnan(L5_ci(1,1,:,1))));
    axis_id = axis;
end

makefigure();

subplot(2,4,1);
p1 = plot(axis,L5_sEst(1,axis_id),'lineWidth',2); hold on;
p2 = plot(axis,L5_sEst(2,axis_id),'lineWidth',2); 
title('L5 noise');
legend('A','B')
ciplot(squeeze(L5_ci(1,1,axis_id,1)),squeeze(L5_ci(1,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(L5_ci(1,2,axis_id,1)),squeeze(L5_ci(1,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,2);
plot(axis,L5_vmuEst(1,axis_id),'lineWidth',2); hold on;
plot(axis,L5_vmuEst(2,axis_id),'lineWidth',2); 
title('L5 volatility');
legend('A','B')
ciplot(squeeze(L5_ci(2,1,axis_id,1)),squeeze(L5_ci(2,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(L5_ci(2,2,axis_id,1)),squeeze(L5_ci(2,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,3);
plot(axis,L5_muEst(1,axis_id),'lineWidth',2); hold on;
plot(axis,L5_muEst(2,axis_id),'lineWidth',2); 
title('L5 mu');
legend('A','B')
ciplot(squeeze(L5_ci(3,1,axis_id,1)),squeeze(L5_ci(3,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(L5_ci(3,2,axis_id,1)),squeeze(L5_ci(3,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,4);
plot(axis,L5_mean(1,axis_id),'lineWidth',2); hold on;
plot(axis,L5_mean(2,axis_id),'lineWidth',2); 
title('L5 mean');
legend('A','B')
ciplot(squeeze(L5_ci(4,1,axis_id(8:(length(axis_id)-1)),1)),squeeze(L5_ci(4,1,axis_id(8:(length(axis_id)-1)),2)),axis(8:(length(axis_id)-1)),p1.Color,0.3);
ciplot(squeeze(L5_ci(4,2,axis_id(1:(length(axis_id)-1)),1)),squeeze(L5_ci(4,2,axis_id(1:(length(axis_id)-1)),2)),axis(1:(length(axis_id)-1)),p2.Color,0.3);

subplot(2,4,5);
plot(axis,M10_sEst(1,axis_id),'lineWidth',2); hold on;
plot(axis,M10_sEst(2,axis_id),'lineWidth',2); 
title('M10 noise');
legend('A','B')
ciplot(squeeze(M10_ci(1,1,axis_id,1)),squeeze(M10_ci(1,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(M10_ci(1,2,axis_id,1)),squeeze(M10_ci(1,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,6);
plot(axis,M10_vmuEst(1,axis_id),'lineWidth',2); hold on;
plot(axis,M10_vmuEst(2,axis_id),'lineWidth',2); 
title('M10 volatility');
legend('A','B')
ciplot(squeeze(M10_ci(2,1,axis_id,1)),squeeze(M10_ci(2,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(M10_ci(2,2,axis_id,1)),squeeze(M10_ci(2,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,7);
plot(axis,M10_muEst(1,axis_id),'lineWidth',2); hold on;
plot(axis,M10_muEst(2,axis_id),'lineWidth',2); 
title('M10 mu');
legend('A','B')
ciplot(squeeze(M10_ci(3,1,axis_id,1)),squeeze(M10_ci(3,1,axis_id,2)),axis,p1.Color,0.3);
ciplot(squeeze(M10_ci(3,2,axis_id,1)),squeeze(M10_ci(3,2,axis_id,2)),axis,p2.Color,0.3);

subplot(2,4,8);
plot(axis,M10_mean(1,axis_id),'lineWidth',2); hold on;
plot(axis,M10_mean(2,axis_id),'lineWidth',2); 
title('M10 mean');
legend('A','B')
ciplot(squeeze(M10_ci(4,1,axis_id(8:(length(axis_id)-1)),1)),squeeze(M10_ci(4,1,axis_id(8:(length(axis_id)-1)),2)),axis(8:(length(axis_id)-1)),p1.Color,0.3);
ciplot(squeeze(M10_ci(4,2,axis_id(1:(length(axis_id)-1)),1)),squeeze(M10_ci(4,2,axis_id(1:(length(axis_id)-1)),2)),axis(1:(length(axis_id)-1)),p2.Color,0.3);

exportgraphics(gcf,[data_path+data + "_lowvmu.png"],'Resolution',300)

%% spaghetti plot
for participant = 1:size(alldata.subjectdata,2)
    line = plot(axis, squeeze(M10(4,1,participant,axis_id)),'b--','lineWidth',2); hold on;
    line.Color(4) = 0.6;
end

for participant = 1:size(alldata.subjectdata,2)
    line = plot(axis, squeeze(M10(4,2,participant,axis_id)),'r','lineWidth',2); hold on;
    line.Color(4) = 0.3;
end

set(gcf,'Position',[50 50 700 200])
exportgraphics(gcf,[data_path + data + "_spaghetti_lowvmu.png"],'Resolution',300)
%% t test prep
addpath('/Users/rh/Documents/GitHub/drumtrainer/analysis/functions/');
if data ~= "all"
    if realign
        axis = 50:100;
    else axis = 1:50;
    end
    
    L5_sEst_par = nan(2,20);
    for i = 1:2
        temp = mean(L5(1,i,:,axis),4,'omitnan');
        L5_sEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end
    
    L5_vmuEst_par = nan(2,20);
    for i = 1:2
        temp = mean(L5(2,i,:,axis),4,'omitnan');
        L5_vmuEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    L5_muEst_par = nan(2,20);
    for i = 1:2
        temp = mean(L5(3,i,:,axis),4,'omitnan');
        L5_muEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    L5_mean_par = nan(2,20);
    for i = 1:2
        temp = mean(L5(4,i,:,axis),4,'omitnan');
        L5_mean_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    M10_sEst_par = nan(2,20);
    for i = 1:2
        temp = mean(M10(1,i,:,axis),4,'omitnan');
        M10_sEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    M10_vmuEst_par = nan(2,20);
    for i = 1:2
        temp = mean(M10(2,i,:,axis),4,'omitnan');
        M10_vmuEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    M10_muEst_par = nan(2,20);
    for i = 1:2
        temp = mean(M10(3,i,:,axis),4,'omitnan');
        M10_muEst_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    M10_mean_par = nan(2,20);
    for i = 1:2
        temp = mean(M10(4,i,:,axis),4,'omitnan');
        M10_mean_par(i,(1:sum(~isnan(temp)))) = temp(~isnan(temp));
    end

    subplot(4,2,1);
    nbp = notBoxPlot(L5_sEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest(L5_sEst_par(1,:),L5_sEst_par(2,:));
    title(["L5 noise, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,2);
    nbp = notBoxPlot(M10_sEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(M10_sEst_par(1,:),M10_sEst_par(2,:),'Vartype','unequal')
    title(["M10 noise, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,3);
    nbp = notBoxPlot(L5_vmuEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(L5_vmuEst_par(1,:),L5_vmuEst_par(2,:),'Vartype','unequal')
    title(["L5 volatility, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,4);
    nbp = notBoxPlot(M10_vmuEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(M10_vmuEst_par(1,:),M10_vmuEst_par(2,:),'Vartype','unequal')
    title(["M10 volatility, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,5);
    nbp = notBoxPlot(L5_muEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(L5_muEst_par(1,:),L5_muEst_par(2,:),'Vartype','unequal')
    title(["L5 mu, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,6);
    nbp = notBoxPlot(M10_muEst_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(M10_muEst_par(1,:),M10_muEst_par(2,:),'Vartype','unequal')
    title(["M10 mu, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,7);
    nbp = notBoxPlot(L5_mean_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(L5_mean_par(1,:),L5_mean_par(2,:),'Vartype','unequal')
    title(["L5 mean, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    subplot(4,2,8);
    nbp = notBoxPlot(M10_mean_par');formatNBP(nbp);
    [H,P,CI,STATS] = ttest2(M10_mean_par(1,:),M10_mean_par(2,:),'Vartype','unequal')
    title(["M10 mean, t = "+ num2str(round(STATS.tstat,3)) + ", p = "+ char(num2str(round(P,3)))]);

    set(gcf,'Position',[50 50 700 700])

    exportgraphics(gcf,[data + "_boxplot.png"],'Resolution',300)

end
%% repeated measure ANOVA
if data == "all"
    axis_pre = 1:50;
    axis_post = 51:100;
    
    L5_sEst_par = nan(2,2,20); %group, pre/post, participant
    for i = 1:2 %group
        temp_pre = mean(squeeze(L5(1,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(L5(1,i,:,axis_post)),2,'omitnan');
        L5_sEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        L5_sEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    L5_sEst = reshape(L5_sEst_par,[4 20])';

    L5_vmuEst_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(L5(2,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(L5(2,i,:,axis_post)),2,'omitnan');
        L5_vmuEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        L5_vmuEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    L5_vmuEst = reshape(L5_vmuEst_par,[4 20])';

    L5_muEst_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(L5(3,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(L5(3,i,:,axis_post)),2,'omitnan');
        L5_muEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        L5_muEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    L5_muEst = reshape(L5_muEst_par,[4 20])';
    
    L5_mean_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(L5(4,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(L5(4,i,:,axis_post)),2,'omitnan');
        L5_mean_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        L5_mean_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    L5_mean = reshape(L5_mean_par,[4 20])';

    M10_sEst_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(M10(1,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(M10(1,i,:,axis_post)),2,'omitnan');
        M10_sEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        M10_sEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    M10_sEst = reshape(M10_sEst_par,[4 20])';

    M10_vmuEst_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(M10(2,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(M10(2,i,:,axis_post)),2,'omitnan');
        M10_vmuEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        M10_vmuEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    M10_vmuEst = reshape(M10_vmuEst_par,[4 20])';

    M10_muEst_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(M10(3,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(M10(3,i,:,axis_post)),2,'omitnan');
        M10_muEst_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        M10_muEst_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    M10_muEst = reshape(M10_muEst_par,[4 20])';

    M10_mean_par = nan(2,2,20);
    for i = 1:2
        temp_pre = mean(squeeze(M10(4,i,:,axis_pre)),2,'omitnan');
        temp_post = mean(squeeze(M10(4,i,:,axis_post)),2,'omitnan');
        M10_mean_par(i,1,(1:sum(~isnan(temp_pre)))) = temp_pre(~isnan(temp_pre));
        M10_mean_par(i,2,(1:sum(~isnan(temp_post)))) = temp_post(~isnan(temp_post));
    end
    M10_mean = reshape(M10_mean_par,[4 20])';

subplot(2,4,1);
nbp = notBoxPlot(L5_sEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("L5 noise");

subplot(2,4,2);
nbp = notBoxPlot(L5_vmuEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("L5 volatility");

subplot(2,4,3);
nbp = notBoxPlot(L5_muEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("L5 mu");

subplot(2,4,4);
nbp = notBoxPlot(L5_mean,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("L5 mean");

subplot(2,4,5);
nbp = notBoxPlot(M10_sEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("M10 noise");

subplot(2,4,6);
nbp = notBoxPlot(M10_vmuEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("M10 volatility");

subplot(2,4,7);
nbp = notBoxPlot(M10_muEst,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("M10 mu");

subplot(2,4,8);
nbp = notBoxPlot(M10_mean,[1,3,2,4]);%1"A pre",2"B pre", 3"A post", 4"B post"
set(gca,'xticklabel',{'A pre','A post','B pre','B post'});
formatNBP(nbp);
title("M10 mean");

    set(gcf,'Position',[50 50 1000 600])

%     exportgraphics(gcf,[data + "_boxplot.png"],'Resolution',300)
    
end

%% heatmap to see parameter estimation
participant =12;
a = alldata.subjectdata(participant).learner.vmuDist;
axis2 = alldata.subjectdata(participant).learner.vmulog;
axis1 = alldata.subjectdata(participant).allData.Day(1):(alldata.subjectdata(participant).allData.Day(end)+1);
heatmap(axis1,axis2, a');
