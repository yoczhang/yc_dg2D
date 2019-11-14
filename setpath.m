%% Set path to FEM  
%
% Add all subdirectories under pathstr to search path.
%
%

clear 
close all
clc
setpath_pwd = pwd;
addpath(genpath(setpath_pwd),'-begin');

saveFilename = [setpath_pwd,'/setpath_pwd/setpath_pwd.mat'];
save(saveFilename, 'setpath_pwd');
clear

% % % savepath;
