%
clear;
close all;
%addpath("C:\Users\kmiya\Documents\GitHub\J349546")
addpath("C:\Users\kmiya\Documents\GitHub\J349554")

% File path
filename = 'J349554.xlsx';

% Read the entire Excel file (modify sheet name or range if needed)
data = readtable(filename);

% Extract and rename variables with meaningful names
own97        = data.ER10035;
housevalue97 = data.ER10036/69.557;
moved97      = data.ER10072;

own99        = data.ER13040;
housevalue99 = data.ER13041/74.151;
moved99      = data.ER13077;
sold99 = data.ER15046;

own01        = data.ER17043;
housevalue01 = data.ER17044/80.99;
moved01      = data.ER17088;
sold01 = data.ER19242;

% (Optional) Combine into a single table
housingData = table( ...
   own97, housevalue97, moved97, ...
   own99, housevalue99, moved99, ...
   own01,housevalue01, moved01, ...
    sold99, price99, ...
    sold01, price01);

% Save to .mat or write to file if needed
%save('housingData.mat', 'housingData');

%% Filter

filterOwn = all([housingData.own97 == 1 , housingData.own99 == 1, housingData.own01 == 1],2);

housingDataOwn = housingData(filterOwn,:);
thr = .15;

%% 1999
filterSold99 = all([housingDataOwn.moved99 == 1,housingDataOwn.sold99 == 1],2);

dataTemp = housingDataOwn(filterSold99,:);


filterAdj99 = any([abs(dataTemp.housevalue99./dataTemp.housevalue97) > 1+thr,...
    abs(dataTemp.housevalue99./dataTemp.housevalue97) < 1-thr],2);

dataTemp = dataTemp(filterAdj99,:);

sizeAdj99 = dataTemp.housevalue99./dataTemp.housevalue97;

%% 2001
filterSold01 = all([housingDataOwn.moved01 == 1,housingDataOwn.sold01 == 1],2);

dataTemp = housingDataOwn(filterSold01,:);


filterAdj01 = any([abs(dataTemp.housevalue01./dataTemp.housevalue99) > 1+thr,...
    abs(dataTemp.housevalue01./dataTemp.housevalue99) < 1-thr],2);

dataTemp = dataTemp(filterAdj01,:);

sizeAdj01 = dataTemp.housevalue01./dataTemp.housevalue99;

%%
figure(1)
histogram([sizeAdj99;sizeAdj01],40,'BinEdges',[linspace(0,1-thr,20) linspace(1+thr,4,20)]);
title('Distribution of size of housing adjustment ')
xlabel('Size of Adjustment')
ylabel('Frequency')








