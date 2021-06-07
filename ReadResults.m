clear;
clc;
addpath('Results\BLOG_L\20\');
cd Results\BLOG_L\20
S = dir('*.txt');
N = {S.name};
cd ..
format short;

% AUCs = zeros(10,size(N,2));
% AUCPRs = zeros(10,size(N,2));
% Acc = zeros(10,size(N,2));
MacroF1s = zeros(8,size(N,2));
MicroF1s = zeros(8,size(N,2));
% 
for j = 1: size(N,2)
    fileID = fopen(N{j});
    formatSpec = '%s';
       % while ~feof(fileID)

        C = textscan(fileID,formatSpec,'Delimiter','\t', 'HeaderLines',0);
        %end
    fclose(fileID);

    for i =1:size(C{1},1)

        A = strsplit(C{1}{i});
%         AUCs(i,j) = str2double(A{3});
%         AUCPRs(i,j) = str2double(A{5});
%         Acc(i,j) = str2double(A{7});
         MacroF1s(i,j) = str2double(A{3});
         MicroF1s(i,j) = str2double(A{5});

    end
end
%---------------------------------------
% 
% % 
 %%rmpath('Results\CORA\20\');
 cd ..
 cd ..
 cd C:\Users\Secil\Desktop\BioInfo\DGI-master\PreProcess\Results\BoxPlotTool\code

 save('BLOG_L_20.mat','MacroF1s');
 %Come back
 cd C:\Users\Secil\Desktop\CIKM\LANE-master
 %  cd BoxPlotTool\code\
%  save('DDIPartition8AUC.mat','AUCs');
% 
% % 
