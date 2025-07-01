clear;
close all;

addpath("C:\Users\kmiya\Documents\GitHub\J349554");

% File path
filename = 'J349554.xlsx';

% Read all data
data = readtable(filename);

% Convert from nominal to numeric if needed
% data.ER13040 = double(data.ER13040);  % uncomment if variables are stored as categorical

% List of years and corresponding inflation adjustment factors
years = [1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023];

price_index = [61.833,	63.852	,66.355	,68.976,	70.715,	72.503,	74.360,	76.968	,79.712	,81.888,	83.320,	83.358,	84.478	,86.309	,88.308,	90.742,	93.536,	96.666,	100.000	,103.391,	106.915,	109.978	,112.744,	119.342,	128.362];
indx = [true false true false true false true false true false true false true false true false  true false  true false true false  true false true];
price_index = price_index(indx);

n = length(years);
own = cell(n,1);
value = cell(n,1);
moved = cell(n,1);
sold = cell(n,1);
age  = cell(n,1);

% Variable roots and suffixes for each type
own_vars   = ["ER13040", "ER17043", "ER21042", "ER25028", "ER36028", "ER42029", "ER47329", "ER53029", "ER60030", "ER66030", "ER72030", "ER78031", "ER82032"];
value_vars = ["ER13041", "ER17044", "ER21043", "ER25029", "ER36029", "ER42030", "ER47330", "ER53030", "ER60031", "ER66031", "ER72031", "ER78032", "ER82033"];
moved_vars = ["ER13077", "ER17088", "ER21117", "ER25098", "ER36103", "ER42132", "ER47440", "ER53140", "ER60155", "ER66156", "ER72156", "ER78158", "ER82141"];
sold_vars  = ["ER15046", "ER19242", "ER22637", "ER26618", "ER37636", "ER43627", "ER48972", "ER54734", "ER61845", "ER67899", "ER73927", "ER80049", "ER84019"];
age_vars   = ["ER13010", "ER17013", "ER21017", "ER25017", "ER36017", "ER42017", "ER47317", "ER53017", "ER60017", "ER66017", "ER72017", "ER78017", "ER82018"];

% Extract variables year by year
for it = 1:n
    yr = years(it);

    % Safely extract each variable
    own{it}   = data.(own_vars(it));
    value{it} = data.(value_vars(it)) / price_index(it);   % inflation-adjusted house value
    moved{it} = data.(moved_vars(it));
    sold{it}  = data.(sold_vars(it));
    age{it}   = data.(age_vars(it));
end

% Convert to a named table with year suffixes
housingData = table();
for it = 1:n
    suffix = num2str(years(it));
    housingData.(['own' suffix])   = own{it};
    housingData.(['value' suffix]) = value{it};
    housingData.(['moved' suffix]) = moved{it};
    housingData.(['sold' suffix])  = sold{it};
    housingData.(['age' suffix])  = sold{it};
end

% Optional: Save as .mat file
%save('psid_housing_clean.mat', 'housingData');


%% Filter
thr = 0.15;  % threshold for meaningful adjustment
sizeAdjAllYears = [];
dist_adj = cell(length(years) - 1,1);
prop_adj_up = zeros(length(years) - 1,1);
for it = 1:(length(years) - 1)
    year_t = num2str(years(it));
    year_tp1 = num2str(years(it+1));

    % Extract relevant variables
    own_t = housingData.(['own' year_t]);
    age_t = housingData.(['age' year_t]);
    own_tp1 = housingData.(['own' year_tp1]);
    moved_tp1 = housingData.(['moved' year_tp1]);
    sold_tp1 = housingData.(['sold' year_tp1]);
    val_t = housingData.(['value' year_t]);
    val_tp1 = housingData.(['value' year_tp1]);

    % Step 1: Owners in both years and younger than 65 in both years
    filterOwn = all([own_t == 1, own_tp1 == 1, age_t < 63], 2);
    dataOwn = table(own_t, own_tp1, moved_tp1, sold_tp1, val_t, val_tp1);
    dataOwn = dataOwn(filterOwn, :);

    % Step 2: Moved and sold in later year
    filterTransact = all([dataOwn.moved_tp1 == 1, dataOwn.sold_tp1 == 1], 2);
    dataAdj = dataOwn(filterTransact, :);

    % Step 3: Remove missing, zero, or negative values
    validVals = dataAdj.val_t > 0 & dataAdj.val_tp1 > 0 & ...
                ~isnan(dataAdj.val_t) & ~isnan(dataAdj.val_tp1);
    dataAdj = dataAdj(validVals, :);

    % Step 4: Compute relative adjustment
    sizeRatio = dataAdj.val_tp1 ./ dataAdj.val_t;

    % Step 5: Keep meaningful adjustments
    filterSize = abs(sizeRatio - 1) > thr;
    sizeAdj = sizeRatio(filterSize);

    % Step 6: Stack results
    sizeAdjAllYears = [sizeAdjAllYears; sizeAdj];
    
    dist_adj{it} = sizeAdj;
    prop_adj_up(it) = mean(sizeAdj>1);
end


%%
numBins =40;
figure(1)
histogram(sizeAdjAllYears,numBins,'BinEdges',[linspace(0.1,1-thr,numBins/2) linspace(1+thr,4,numBins/2)]);
title('Distribution of size of housing adjustment ')
xlabel('Size of Adjustment')
ylabel('Frequency')

figure(2)
for it = 1:(length(years) - 1)
    subplot(3,4,it)
    histogram(dist_adj{it},numBins,'BinEdges',[linspace(0.1,1-thr,numBins/2) linspace(1+thr,4,numBins/2)]);
    title(['Year ',num2str(years(it+1))])
end
sgtitle('Distribution of adjustment across years')






