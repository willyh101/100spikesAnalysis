%%
% Recreates the key plots of the paper (averaged over ensembles), and 
% reproducing them using a cell-by-cell analysis 
%%
clear; close all; clc;

%%
% Adds all subfolders to the path
restoredefaultpath;
folder = fileparts(which('cellByCellAnalysis_GH.m')); 
addpath(genpath(folder));
rmpath(folder);

%%
load('./compressedData/cellTable220816.mat')

%% Cell conditions used in the functions
cellCond = cellTable.offTarget==0; 
cellCondTuned = cellTable.offTarget==0 & cellTable.visP<0.05 & cellTable.cellOSI > 0.25;
cellCondNonVis = cellTable.offTarget==0 & cellTable.visP>0.05;

%% Figure 2: min distance plot
fprintf('-------------------------\n')
fprintf('Creating Fig 2: min dist vs. dF/F plot\n')
Fig2(cellTable,cellCond)
fprintf('-------------------------\n')
% Fig2_cbc(cellTable,cellCond)


%% Figure 3: Effect of the spread of ensemble
fprintf('Creating Fig 3: ensemble spread vs. dF/F plot\n')
Fig3(cellTable,cellCond);
fprintf('-------------------------\n')

%% Figure 4: Iso vs. ortho
% Note: the cells used here are tuned 
fprintf('Creating Fig 4: Iso vs. Ortho\n')
Fig4(cellTable,cellCondTuned);
% Fig4_cbc(cellTable,cellCondTuned);
fprintf('-------------------------\n')

%% Figure 5: Tight co-tuned investigation
fprintf('Creating Fig 5: holo-squares\n')
Fig5a(cellTable,cellCond);
Fig5bc(cellTable,cellCondTuned,cellCondNonVis)
fprintf('-------------------------\n')
% Fig5_cbc(cellTable,cellCond);
% Fig5_cbc_dataFit(cellTable,cellCond);
