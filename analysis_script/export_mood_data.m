clear; clc;
addpath('/Users/rh/Desktop/BayesM10/');
data_path = '/Users/rh/Desktop/BayesM10/data/';

file_name = 'movement_data.csv';
movement_data = readtable([data_path file_name]);

kept_cols = ["ID","Day","pos","neg","mood_sum"];

mood_data = movement_data(:,kept_cols);

writetable(mood_data,[data_path+"data_mood_allparticipants.csv"]);