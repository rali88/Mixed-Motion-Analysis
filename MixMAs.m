function MixMAs
%% Main script that decides what to do based on choices made by MixMAsChoicesDialog
[choice1, choice2] = MixMAsChoicesDialog;
%% Based on 1-D choices actions to do
if choice1 == 1
    if choice2 == 1
        [SimulationResults_1D, X_1D, Y_1D] = MixMSi1D;
        clear choice1; clear choice2;
        save('SimulationResults_1D.mat')
    elseif choice2 == 2
        [CutoffDisplacement_1D,DynamicMissedAtCutoffDisplacement_1D,StaticMissedAtCutoffDisplacement_1D] = FindDynamicStaticCutoff1D;
        clear choice1; clear choice2;
        save('DynamicStaticCutoff_1D.mat')
    elseif choice2 == 3
        [CutoffAngle_1D,DiffusiveMissedAtCutoffAngle_1D,DirectedMissedAtCutoffAngle_1D] = FindDirectedDiffusiveCutoff1D;
        clear choice1; clear choice2;
        save('DirectedDiffusiveCutoff_1D.mat')
    elseif choice2 == 4
        Results_1D = MixMAsMain;
        clear choice1; clear choice2;
        save('AnalysisResults_1D.mat')
    end
%% Based on 2-D choices actions to do    
else
    if choice2 == 1
        [SimulationResults_2D, X_2D, Y_2D] = MixMSi;
        clear choice1; clear choice2;
        save('SimulationResults_2D.mat')
    elseif choice2 == 2
        [CutoffDisplacement_2D,DynamicMissedAtCutoffDisplacement_2D,StaticMissedAtCutoffDisplacement_2D] = FindDynamicStaticCutoff;
        clear choice1; clear choice2;
        save('DynamicStaticCutoff_2D.mat')
    elseif choice2 == 3
        [CutoffAngle_2D,DiffusiveMissedAtCutoffAngle_2D,DirectedMissedAtCutoffAngle_2D] = FindDirectedDiffusiveCutoff;
        clear choice1; clear choice2;
        save('DirectedDiffusiveCutoff_2D.mat')
    elseif choice2 == 4
        Results_2D = MixMAsMain;
        save('AnalysisResults_2D.mat')
    end
end
end
%%
function [choice1, choice2] = MixMAsChoicesDialog
%% Gui for MixMAs

%% Creating dialog and buttons for initial selection
d1 = dialog('Position',[300 300 250 190],'Name','Select dimensionality');

btn1 = uicontrol('Parent',d1,...
    'Position',[50 50 150 50],...
    'String','1-D Motion',...
    'Callback',@btn1_callback );
btn2 = uicontrol('Parent',d1,...
    'Position',[50 120 150 50],...
    'String','2-D Motion',...
    'Callback',@btn2_callback);

uiwait(d1)

    function btn1_callback(btn1,callbackdata)
        if  btn1.Value == 1;
            choice1 = 1;
            close(d1)
        end
    end
    function btn2_callback(btn2,callbackdata)
        if  btn2.Value == 1;
            choice1 = 2;
            close(d1)
        end
    end

%% Creating dialog and buttons for further selection

d2 = dialog('Position',[300 300 250 360],'Name','What to do');

Btn1 = uicontrol('Parent',d2,...
    'Position',[50 50 150 50],...
    'String','Simulate motion',...
    'Callback',@Btn1_callback );
Btn2 = uicontrol('Parent',d2,...
    'Position',[50 120 150 50],...
    'String','Find dynamic/static cutoff',...
    'Callback',@Btn2_callback);
Btn3 = uicontrol('Parent',d2,...
    'Position',[50 190 150 50],...
    'String','Find diffusive/directed cutoff',...
    'Callback',@Btn3_callback);
Btn4 = uicontrol('Parent',d2,...
    'Position',[50 260 150 50],...
    'String','Analyze coordinates file',...
    'Callback',@Btn4_callback);
uiwait(d2)

    function Btn1_callback(Btn1,callbackdata)
        if  Btn1.Value == 1;
            choice2 = 1;
            close(d2)
        end
    end
    function Btn2_callback(Btn2,callbackdata)
        if  Btn2.Value == 1;
            choice2 = 2;
            close(d2)
        end
    end
    function Btn3_callback(Btn3,callbackdata)
        if  Btn3.Value == 1;
            choice2 = 3;
            close(d2)
        end
    end
    function Btn4_callback(Btn4,callbackdata)
        if  Btn4.Value == 1;
            choice2 = 4;
            close(d2)
        end
    end
end
%%
function Results = MixMAsMain
%% Analyzes motion of particles with multiple components of motion
% Main function
%% Requesting user input
whattodo = input('To analyze excel track file press 1, To fit MSD to time press 2, To exit press 3 ');
if whattodo == 1
    %% Making vector for alphabets and numbers
    alphabets = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'];
    numbers = [1:length(alphabets)];
    %% selecting excel file and getting column info from user
    [filetoread, path] = uigetfile({'*.xls;*.xlsx;.csv', 'Required files (*.xls,*.xlsx,*.csv)'},'Select Excel File');
    tid = input('Enter TID Column title  ','s');
    it = find(alphabets == tid);
    tidno = numbers(it);
    xcord = input('Enter X coordinate Column Title  ','s');
    ix = find(alphabets == xcord);
    xcordno = numbers(ix);
    ycord = input('Enter Y coordinate Column Title  ','s');
    iy = find(alphabets == ycord);
    ycordno = numbers(iy);
    tcord = input('Enter T coordinate Column Title  ','s');
    ip = find(alphabets == tcord);
    tcordno = numbers(ip);
    advanced_options = input('Enter 1 for advanced options (not recommended) else enter 0 ');
    if advanced_options == 0
        minint = 2;
        sizeofslidingwindow = 1;
    else
        minint = input('Enter minimum number of intervals that a particle needs to have ');
        sizeofslidingwindow = input('Enter the size of the sliding window  ');
    end
    noofconsecchanges = 2;
    diffusioncutoff = input('enter dynamic cutoff in microns  ');
    directedcutoff = input('enter direction cutoff in std of degrees  ');
    timeinterval = input('enter time between frames in seconds ');
    microncheck = input('enter 1 if coordinates are in microns, else enter 0 ');
    if microncheck == 1
        micronfactor = 1;
    else
        micronfactor = input('enter pixel size (microns/pixel)');
    end
    clc
    %% Reading excel track file
    [TID, X, Y, T] = ReadTrackFile(filetoread, path, tidno,xcordno,ycordno,minint,tcordno);
    %% Finding points of state changes
    [state2, changepoint] = FindChangePoints(X, Y,sizeofslidingwindow,diffusioncutoff,directedcutoff);
    %% Calculating output parameters
    [diffusivedwelltime, staticdwelltime, directeddwelltime, numberdiffusive, numberstatic,...
        numberdirected,numberdiffusivetostatic, numberdiffusivetodirected, numberstatictodiffusive,...
        numberstatictodirected, numberdirectedtodiffusive, numberdirectedtostatic, dDelta,sDelta,...
        rDelta,dMSD,sMSD,rMSD...
        ] = CalculateParameters(state2, changepoint, X, Y, timeinterval, micronfactor);
    save(['DiffusiveMeanSquaredDisplacement',filetoread,'.mat'],'dMSD')
    save(['StaticMeanSquaredDisplacement',filetoread,'.mat'],'sMSD')
    save(['DirectedMeanSquaredDisplacement',filetoread,'.mat'],'rMSD')
    Results = struct('DiffusiveDwellTimeInSec',diffusivedwelltime,'StaticDwellTimeInSec', staticdwelltime,...
        'DirectedDwellTimeInSec', directeddwelltime, 'NumberOfDiffusiveEvents', numberdiffusive,...
        'NumberOfStaticEvents', numberstatic, 'NumberOfDirectedEvents', numberdirected,'DiffusiveToStaticNumber',numberdiffusivetostatic,...
        'DiffusiveToDirectedNumber',numberdiffusivetodirected,'StaticToDiffusiveNumber', numberstatictodiffusive,...
        'StaticToDirectedNumber', numberstatictodirected, 'DirectedToDiffusiveNumber', numberdirectedtodiffusive,...
        'DirectedToStaticNumber', numberdirectedtostatic, 'DiffusiveMeanSquaredDisplacement', {dMSD},...
        'StaticMeanSquaredDisplacement', {sMSD},'DirectedMeanSquaredDisplacement', {rMSD},'StateValues', {state2}, 'ChangePoints', {changepoint});
elseif whattodo == 2
    %% Fitting MSD to time
    [MeanSquaredDisplacement, StdMeanSquaredDisplacement,TimeInterval, DiffusionCoefficient, Alpha, Rsquared,fitresult,gof] = FitMsd2Time;
    Results = struct('MeanSquaredDisplacement',MeanSquaredDisplacement, 'StandardErrorOfMSD',StdMeanSquaredDisplacement,'TimeInterval',...
        TimeInterval, 'DiffusionCoefficientInMicorns', DiffusionCoefficient, 'AlphaValue', Alpha,'RSquared', Rsquared,...
        'ExtendedResults',fitresult,'GoodnessOfFit',gof);
elseif whattodo == 3
    Results = [];
end
end
%%
function [TID, X, Y, T] = ReadTrackFile(filetoread, path, tidno,xcordno,ycordno,minint,tcordno)
%% Read excel file containing tracking information
%nested function 1
%% Getting user input
if nargin == 0
    %% Making vector for alphabets and numbers
    alphabets = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'];
    numbers = [1:length(alphabets)];
    %% Reading excel file and sorting data in vectors
    [filetoread, path] = uigetfile('*','Select Excel File');
    tid = input('Enter TID Column title  ','s');
    it = find(alphabets == tid);
    tidno = numbers(it);
    xcord = input('Enter X coordinate Column Title  ','s');
    ix = find(alphabets == xcord);
    xcordno = numbers(ix);
    ycord = input('Enter Y coordinate Column Title  ','s');
    iy = find(alphabets == ycord);
    ycordno = numbers(iy);
    tcord = input('Enter T coordinate Column Title  ','s');
    ip = find(alphabets == tcord);
    tcordno = numbers(ip);
    minint = input('Enter minimum number of intervals  ');
    clc
end
Data = xlsread([path,filetoread]); % reading excel file
tid = Data(:,tidno); tid = tid'; % getting id no for particles
xcord = Data(:,xcordno); xcord = xcord'; % getting x-coordinates
ycord = Data(:,ycordno); ycord = ycord'; % getting y-coordinates
% adding exception, if there is not t coordinate 
if size(Data,2) > 3
    tcord = Data(:,tcordno); tcord = tcord'; % getting t coordinates
else
    tcord = zeros(1,length(xcord));
end
%% Sorting data into cells according to the tid no
uniquetids = unique(tid);
count = 0;
for k = 1 : length(uniquetids)
    locationoftid = find(tid == uniquetids(k));
    if length(locationoftid) >= minint + 1
        count = count + 1;
        TID{count} = tid(locationoftid);
        X{count} = xcord(locationoftid);
        Y{count} = ycord(locationoftid);
        T{count} = tcord(locationoftid);
    end
end
end
%%
function [stateF, changepointF] = FindChangePoints(X,Y,SlidingWindowSize,dynamiccutoff,directedcutoff)
%% Find points of transition and their respective states
% nested function 2
%% Calculating displacement and direction for every interval
% getting user input
if nargin == 2
    SlidingWindowSize = input('enter size of slideing window  ');
    dynamiccutoff = input('enter dynamic cutoff in microns  ');
    directedcutoff = input('enter direction cutoff in std of degrees  ');
end
%% Generating direction and angle values vector
for k = 1 : length(X)
    for j = 1 : length(X{k}) - 1
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand(rise/run);
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end
%% Running sliding window
for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
%% Generating state vectors
for k = 1 : length(displace2)
    state{k} = [];
    for j = 1 : length(displace2{k})
        if displace2{k}(j) < dynamiccutoff
            state{k}(j) = 2;
        else
            if direction2{k}(j) < directedcutoff
                state{k}(j) = 3;
            else
                state{k}(j) = 1;
            end
        end
    end
end

clear displace2; clear direction2;
%% Combining same state values and generating a final state and transition points vector
for count = 1 : length(state)
    diffusive = find(state{count} == 1);
    static = find(state{count} == 2);
    directed = find(state{count} == 3);
    a = 1;
    b = 2;
    starts = [];
    ends = [];
    while a <= length(static) - 2
        if ((static(a) + 1) == (static(a+1))) && ((static(a + 1) + 1) == (static(a + 2)))
            starts = [starts, a];
            while (a <= length(static) - 2) && ((static(a + 1) + 1) == (static(a + 2)))
                a = a + 1; b = b + 1;
            end
            ends = [ends, b];
            a = b; b = b + 1;
        end
        a = a + 1; b = b + 1;
    end
    
    a = 1;
    b = 2;
    startd = [];
    endd = [];
    while a <= length(diffusive) - 2
        if ((diffusive(a) + 1) == (diffusive(a+1))) && ((diffusive(a + 1) + 1) == (diffusive(a + 2)))
            startd = [startd, a];
            while (a <= length(diffusive) - 2) && ((diffusive(a + 1) + 1) == (diffusive(a + 2)))
                a = a + 1; b = b + 1;
            end
            endd = [endd, b];
            a = b; b = b + 1;
        end
        a = a + 1; b = b + 1;
    end
    a = 1;
    b = 2;
    startr = [];
    endr = [];
    while a <= length(directed) - 2
        if ((directed(a) + 1) == (directed(a+1))) && ((directed(a + 1) + 1) == (directed(a + 2)))
            startr = [startr, a];
            while (a <= length(directed) - 2) && ((directed(a + 1) + 1) == (directed(a + 2)))
                a = a + 1; b = b + 1;
            end
            endr = [endr, b];
            a = b; b = b + 1;
        end
        a = a + 1; b = b + 1;
    end
    starts2 = static(starts);
    ends2 = static(ends);
    startd2 = diffusive(startd);
    endd2 = diffusive(endd);
    startr2 = directed(startr);
    endr2 = directed(endr);
    combined = sort([starts2 ends2 startd2 endd2 startr2 endr2]);
    statem = [];
    a = 1;
    while a <= length(combined) - 1
        for j = 1 : length(starts2)
            if combined(a) == starts2(j);
                statem = [statem, 2, 0];
                break
            end
        end
        for j = 1 : length(startd2)
            if combined(a) == startd2(j);
                statem = [statem, 1, 0];
                break
            end
        end
        for j = 1 : length(startr2)
            if combined(a) == startr2(j);
                statem = [statem, 3, 0];
                break
            end
        end
        a = a + 2;
    end
    statem = statem(1:length(statem) - 1);
    
    combined2 = [1];
    a = 2; b = 3;
    while b <= length(combined)
        combined2 = [combined2, round((combined(b) - combined(a))/2) + combined(a)];
        a = a + 2; b = b + 2;
    end
    statem(find(statem == 0)) = [];
    combined2 = [combined2, length(X{count})];
    a = 1;
    consec = [];
    while a <= length(statem) - 1
        if statem(a) == statem(a + 1);
            consec = [consec, a];
            while (a <= length(statem) - 1) && (statem(a) == statem(a + 1))
                a = a + 1;
            end
            consec = [consec, a];
        end
        a = a + 1;
    end
    oddnum = 1:2:length(consec) - 1;
    evennum = 2:2:length(consec);
    combined2remove = consec(evennum);
    statem2remove = consec(oddnum);
    combined2(combined2remove) = []; statem(statem2remove) = [];
    stateF{count} = statem; changepointF{count} = combined2;
end
end
%%
function [diffusivedwelltime, staticdwelltime, directeddwelltime, numberdiffusive, numberstatic,...
    numberdirected,numberdiffusivetostatic, numberdiffusivetodirected, numberstatictodiffusive,...
    numberstatictodirected, numberdirectedtodiffusive, numberdirectedtostatic, dDelta, sDelta,rDelta...
    ,dMSD, sMSD,rMSD...
    ] = CalculateParameters(state2, changepoint, X, Y, timeinterval, micronfactor)
%% Calculate required parameters based on transition points and state values
% nested function 3
%% getting user input
if nargin == 4
    timeinterval = input('enter time between frames in seconds ');
    microncheck = input('enter 1 if coordinates are in microns, else enter 0 ');
    if microncheck == 1
        micronfactor = 1;
    else
        micronfactor = input('enter pixel size (microns/pixel)');
    end
end
%% Initializing variables
diffusivedwelltime = [];
staticdwelltime = [];
directeddwelltime = [];
numberdiffusivetostatic = 0;
numberdiffusivetodirected = 0;
numberstatictodiffusive = 0;
numberstatictodirected = 0;
numberdirectedtodiffusive = 0;
numberdirectedtostatic = 0;
numberdiffusive = 0;
numberstatic = 0;
numberdirected = 0;
%% calculating transitions
for k = 1 : length(state2)
    for j = 1 : length(state2{k}) - 1
        if (state2{k}(j) == 1) && (state2{k}(j + 1) == 2)
            numberdiffusivetostatic = numberdiffusivetostatic + 1;
        elseif (state2{k}(j) == 1) && (state2{k}(j + 1) == 3)
            numberdiffusivetodirected = numberdiffusivetodirected + 1;
        elseif (state2{k}(j) == 2) && (state2{k}(j + 1) == 1)
            numberstatictodiffusive = numberstatictodiffusive + 1;
        elseif (state2{k}(j) == 2) && (state2{k}(j + 1) == 3)
            numberstatictodirected = numberstatictodirected + 1;
        elseif (state2{k}(j) == 3) && (state2{k}(j + 1) == 1)
            numberdirectedtodiffusive = numberdirectedtodiffusive + 1;
        elseif (state2{k}(j) == 3) && (state2{k}(j + 1) == 2)
            numberdirectedtostatic = numberdirectedtostatic + 1;
        end
    end
end
%% calculating no of particles for eachstate, dwelltimes and msd
for k = 1 : length(state2)
    for j = 1 : length(state2{k})
        if state2{k}(j) == 1
            numberdiffusive = numberdiffusive + 1;
            diffusivedwelltime = [diffusivedwelltime, (changepoint{k}(j+1) - changepoint{k}(j))*timeinterval];
            [dDelta{numberdiffusive}, dMSD{numberdiffusive}] = msdisplace(X{k}(changepoint{k}(j):changepoint{k}(j+1)),Y{k}(changepoint{k}(j):changepoint{k}(j+1)),micronfactor);
        elseif state2{k}(j) == 2
            numberstatic = numberstatic + 1;
            staticdwelltime = [staticdwelltime, (changepoint{k}(j+1) - changepoint{k}(j))*timeinterval];
            [sDelta{numberstatic}, sMSD{numberstatic}] = msdisplace(X{k}(changepoint{k}(j):changepoint{k}(j+1)),Y{k}(changepoint{k}(j):changepoint{k}(j+1)),micronfactor);
        elseif state2{k}(j) == 3
            numberdirected = numberdirected + 1;
            directeddwelltime = [directeddwelltime, (changepoint{k}(j+1) - changepoint{k}(j))*timeinterval];
            [rDelta{numberdirected}, rMSD{numberdirected}] = msdisplace(X{k}(changepoint{k}(j):changepoint{k}(j+1)),Y{k}(changepoint{k}(j):changepoint{k}(j+1)),micronfactor);
        end
    end
end
% Making sure that output variables are assigned
if numberdiffusive == 0
    dDelta = []; dMSD = [];
end
if numberstatic == 0
    sDelta = []; sMSD = [];
end
if numberdirected == 0
    rDelta = []; rMSD = [];
end
end
%%
function [delta, msd] = msdisplace(x,y,micronfactor)
%% This function calculates the mean squared displacement of a particle
% nested function 4
no_of_points = length(x);
no_of_deltas = length(x) - 1;
if no_of_deltas < 1
    displace2{1} = [];
end
for i = 1 : no_of_deltas
    for j = 1 : no_of_points - i
        displace2{i}(j) = (((x(j) - x(j+i))^2) + ((y(j) - y(j+i))^2))*micronfactor;
    end
end

for k = 1 : length(displace2)
    msd_sum = 0;
    for u = 1 : length(displace2{k})
        msd_sum = msd_sum + displace2{k}(u);
    end
    msd(k) = (msd_sum/(length(displace2{k})));
end
if no_of_deltas < 1
    delta = [];
end
for m = 1 : no_of_deltas
    delta(m) = m;
end
end
%%
function [MeanSquaredDisplacement, StdMeanSquaredDisplacement,TimeInterval, DiffusionCoefficient, Alpha, Rsquared,fitresult,gof] = FitMsd2Time
%% This function takes Mean Squared Displacement from MixMAs and fits it to time interval
% and calculates D and alpha value
% nested function 5
%% Get MSD.mat file generated by MixMAs
[filestoread, path] = uigetfile('.mat','Select one or more MSD data files','MultiSelect','on');
%% Checking if multiple files were selected
if iscell(filestoread) == 1
    numberoffiles = length(filestoread);
    lengthoffiles = 0;
    for mk = 1 : numberoffiles
        importdata([path,filestoread{mk}])
        lengthoffiles = lengthoffiles + (length(filestoread{mk}));
    end
    msdD_C = cell(1,lengthoffiles);
    for mj = 1 : numberoffiles
        msdD_C = cat(2,msdD_C,importdata([path,filestoread{mj}]));
    end
else
    msdD_C = importdata([path,filestoread]);
end
%% Generating length vector for all MSD cells
len_vec = [];
for i = 1 : length(msdD_C)
    len_vec = [len_vec length(msdD_C{i})];
end
%% Selecting only those values of MSD that are less than the minimum length
for j = 1 : min(len_vec)
    sum_msdD = [];
    for i = 1 : length(msdD_C)
        if len_vec(i) >= j
            sum_msdD = [sum_msdD msdD_C{i}(j)];
        end
    end
    delta_msdD{j} = sum_msdD;
end
%% Caluclating mean and std of MSD
for k = 1 : length(delta_msdD)
    temp_vec = delta_msdD{k};
    mean_msdD(k) = mean(temp_vec);
    std_msdD(k) = (std(temp_vec))/(sqrt(length(temp_vec)));
    time_interval(k) = k;
end
%% Removing MSD which only has a single particle
for k = 1 : length(delta_msdD)
    temp_vec = delta_msdD{k};
    toR = [];
    if length(temp_vec) == 1
        toR = [toR, k];
    end
    mean_msdD(toR) = []; std_msdD(toR) = []; time_interval(toR) = [];
end
%% Fitting MSD to time interval
MeanSquaredDisplacement = mean_msdD; StdMeanSquaredDisplacement = std_msdD; TimeInterval = time_interval;
yTofit = MeanSquaredDisplacement; xTofit = TimeInterval; yStd = StdMeanSquaredDisplacement;
initialD = num2str((((yTofit(1))/(xTofit(1)))));
prompt = {'Initial D to fit:','Initial alpha to fit:','Lower Bound For D:','Upper Bound For D:','Lower Bound For alpha:','Upper Bound For alpha:'};
dlg_title = 'Enter initial values to fit';
num_lines = 1;
defaultans = {initialD,'','-inf','inf','0','2'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
[fitresult, gof] = FitMSD2T(xTofit, yTofit, yStd, str2num(answer{1}), str2num(answer{2}),...
    [str2num(answer{3}),str2num(answer{5})], [str2num(answer{4}),str2num(answer{6})]);
DiffusionCoefficient = fitresult.a; Alpha = fitresult.b; Rsquared = gof.rsquare;
end
%%
function [fitresult, gof] = FitMSD2T(Dt, MSD, semMSD, initialD, initialAlpha, lowB, UpB)
%% Fits MSD to time interval
% nested function 6
weights = 1./semMSD;
[xData, yData, weights] = prepareCurveData( Dt, MSD, weights );

% Set up fittype and options.
ft = fittype( 'a*(x^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [initialD initialAlpha];
opts.Weights = weights;
opts.Lower = lowB;
opts.Upper = UpB;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
plot( fitresult, xData, yData );
hold on
h = errorbar(Dt,MSD,semMSD);
h.LineStyle = 'none';
% Label axes
xlabel DeltaT
ylabel MSD
grid off
hold off
end
%%
function [SimulationResults, X, Y] = MixMSi
%% Simulates particle motion in 2D
%% getting user input
NameofFileToWrite = input('Name of Excel file to write ', 's');
NoofParticlesToSimulate = input('Number of particles to simulate ');
DiffusiveStepSize = input('Dynamic step size in um ');
StaticStepSize = 0;
DirectedStepSize = input('Directed step size in um ');
TrackingError = input('Tracking error in um ');
TimeStep = input('Time between frames in s ');
% Initializing variables
StaticRandomAngles = [0:360];
DiffusiveDwellTime  = [];
StaticDwellTime = [];
DirectedDwellTime = [];
D2S = 0;
D2R = 0;
S2D = 0;
S2R = 0;
R2D = 0;
R2S = 0;
%% Loop over particles
for j = 1 : NoofParticlesToSimulate
    clc
    disp(['Simulation Progress ',num2str((j/NoofParticlesToSimulate)*100),' %'])
    DirectedAngles = randsample(StaticRandomAngles,1);
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 0; ID{j}(1) = j;
    % directed
    for k = 2 : 11
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize*cosd(DirectedAngles) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + DirectedStepSize*sind(DirectedAngles) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    % static
    for k = 12 : 21
        ErrorDirection = randsample(StaticRandomAngles,1);
        MoveDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + StaticStepSize*cosd(MoveDirection) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + StaticStepSize*sind(MoveDirection) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    R2S = R2S + 1;
    StaticDwellTime = [StaticDwellTime ,TimeStep*10];
    % diffusive
    for k = 22 : 31
        ErrorDirection = randsample(StaticRandomAngles,1);
        MoveDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DiffusiveStepSize*cosd(MoveDirection) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + DiffusiveStepSize*sind(MoveDirection) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    S2D = S2D + 1;
    DiffusiveDwellTime = [DiffusiveDwellTime ,TimeStep*10];
    % directed
    for k = 32 : 41
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize*cosd(DirectedAngles) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + DirectedStepSize*sind(DirectedAngles) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    D2R = D2R + 1;
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    % diffucive
    for k = 42 : 51
        ErrorDirection = randsample(StaticRandomAngles,1);
        MoveDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DiffusiveStepSize*cosd(MoveDirection) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + DiffusiveStepSize*sind(MoveDirection) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DiffusiveDwellTime = [DiffusiveDwellTime ,TimeStep*10];
    R2D = R2D + 1;
    % static
    for k = 52 : 61
        ErrorDirection = randsample(StaticRandomAngles,1);
        MoveDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + StaticStepSize*cosd(MoveDirection) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + StaticStepSize*sind(MoveDirection) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    D2S = D2S + 1;
    StaticDwellTime = [StaticDwellTime ,TimeStep*10];
    for k = 62 : 71
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (DirectedStepSize*cosd(DirectedAngles)) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (DirectedStepSize*sind(DirectedAngles)) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    S2R = S2R + 1;
end
% Storing results in a structure
SimulationResults = struct('DiffusiveDwellTimeInSec',DiffusiveDwellTime,'StaticDwellTimeInSec', StaticDwellTime,...
    'DirectedDwellTimeInSec', DirectedDwellTime, 'NumberOfDiffusiveEvents', length(DiffusiveDwellTime),...
    'NumberOfStaticEvents', length(StaticDwellTime), 'NumberOfDirectedEvents', length(DirectedDwellTime),...
    'DiffusiveToStaticNumber',D2S,'DiffusiveToDirectedNumber',D2R,'StaticToDiffusiveNumber', S2D,...
    'StaticToDirectedNumber', S2R, 'DirectedToDiffusiveNumber', R2D,'DirectedToStaticNumber', R2S,...
    'DiffusionCoefficient', ((DiffusiveStepSize)^2)/4);
% writing to an excel file
disp(['Writing Coordinate Data To Excel File'])

XX = []; YY = []; IDD = []; TT = [];

for kk = 1 : NoofParticlesToSimulate
    XX = [XX, X{kk}]; YY = [YY, Y{kk}]; IDD = [IDD, ID{kk}]; TT = [TT, T{kk}];
end

CoordinateData(:,1) = IDD';
CoordinateData(:,2) = XX';
CoordinateData(:,3) = YY';
CoordinateData(:,4) = TT';
xlswrite([NameofFileToWrite,'.xlsx'],CoordinateData)

clc
disp(['Done!'])
end
%%
function [SimulationResults, X, Y] = MixMSi1D
%% Simulates particle motion in 1D
% For comments see MixMSi
NameofFileToWrite = input('Name of Excel file to write ', 's');
NoofParticlesToSimulate = input('Number of particles to simulate ');
DiffusiveStepSize = input('Dynamic step size in um ');
DirectedStepSize = input('Directed step size in um ');
TrackingError = input('Tracking error in um ');
TimeStep = input('Time between frames in s ');

StaticRandomAngles = [0:360];
DiffusiveDwellTime  = [];
StaticDwellTime = [];
DirectedDwellTime = [];
D2S = 0;
D2R = 0;
S2D = 0;
S2R = 0;
R2D = 0;
R2S = 0;

for j = 1 : NoofParticlesToSimulate
    clc
    disp(['Simulation Progress ',num2str((j/NoofParticlesToSimulate)*100),' %'])
    DirectedAngles = randsample(StaticRandomAngles,1);
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 0; ID{j}(1) = j;
    for k = 2 : 11
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    for k = 12 : 21
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    R2S = R2S + 1;
    StaticDwellTime = [StaticDwellTime ,TimeStep*10];
    for k = 22 : 31
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (DiffusiveStepSize*(randsample([1,-1],1))) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    S2D = S2D + 1;
    DiffusiveDwellTime = [DiffusiveDwellTime ,TimeStep*10];
    for k = 32 : 41
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    D2R = D2R + 1;
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    for k = 42 : 51
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (DiffusiveStepSize*(randsample([1,-1],1))) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DiffusiveDwellTime = [DiffusiveDwellTime ,TimeStep*10];
    R2D = R2D + 1;
    for k = 52 : 61
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    D2S = D2S + 1;
    StaticDwellTime = [StaticDwellTime ,TimeStep*10];
    for k = 62 : 71
        ErrorDirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize + (TrackingError*cosd(ErrorDirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(ErrorDirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
    DirectedDwellTime = [DirectedDwellTime ,TimeStep*10];
    S2R = S2R + 1;
end

SimulationResults = struct('DiffusiveDwellTimeInSec',DiffusiveDwellTime,'StaticDwellTimeInSec', StaticDwellTime,...
    'DirectedDwellTimeInSec', DirectedDwellTime, 'NumberOfDiffusiveEvents', length(DiffusiveDwellTime),...
    'NumberOfStaticEvents', length(StaticDwellTime), 'NumberOfDirectedEvents', length(DirectedDwellTime),...
    'DiffusiveToStaticNumber',D2S,'DiffusiveToDirectedNumber',D2R,'StaticToDiffusiveNumber', S2D,...
    'StaticToDirectedNumber', S2R, 'DirectedToDiffusiveNumber', R2D,'DirectedToStaticNumber', R2S,...
    'DiffusionCoefficient', ((DiffusiveStepSize)^2)/4);

disp(['Writing Coordinate Data To Excel File'])

XX = []; YY = []; IDD = []; TT = [];

for kk = 1 : NoofParticlesToSimulate
    XX = [XX, X{kk}]; YY = [YY, Y{kk}]; IDD = [IDD, ID{kk}]; TT = [TT, T{kk}];
end

CoordinateData(:,1) = IDD'; 
CoordinateData(:,2) = XX'; 
CoordinateData(:,3) = YY'; 
CoordinateData(:,4) = TT';
xlswrite([NameofFileToWrite,'.xlsx'],CoordinateData)

clc
disp(['Done!'])
end
%%
function  [CutoffAngle,DiffusiveMissedAtCutoffAngle,DirectedMissedAtCutoffAngle] = FindDirectedDiffusiveCutoff
%% Find cutoff value for directed and diffusive states for 2-D motion
% getting user input
NoofParticlesToSimulate = input('Number of particles to simulate ');
NoOfIters = input('Number of times to run the simulation ');
DiffusiveStepSize = input('Diffusive step size in um ');
DirectedStepSize = input('Directed step size in um ');
TrackingError = input('Localization error in um ');
checkwindow = input('Use default window size enter 1, otherwise enter 0 ');
if checkwindow == 1
    SlidingWindowSize = 1;
else
    SlidingWindowSize = input('Size of sliding window ');
end
h = waitbar(0,'0 % Progress...');
%% loop over number of time for simulation to run
for j = 1 : NoOfIters
    waitbar(j/NoOfIters,h,[num2str((j/NoOfIters)*100), ' % progress...'])
    [X,Y, displaceR, directionR] = MixMSiDirected(DirectedStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate); % simulate directed motion
    [X,Y, displaceD, directionD] = MixMSiDiffusive(DiffusiveStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate); % simulate diffucive motion
    stddirectedangles = []; stddiffusiveangles = [];
    %% Calculate standard deviation of angles
    for k = 1 : length(directionR)
        stddirectedangles = [stddirectedangles, directionR{k}];
    end
    for k = 1 : length(directionD)
        stddiffusiveangles = [stddiffusiveangles,directionD{k}];
    end
    
    directederror{j} = []; diffusiveerror{j} = []; cutoffangle{j} = [];
    %% Calculate missed events for different cutoff values
    runfor = 180;
    for k = 0 : 0.1 : runfor
        directederror{j} = [directederror{j},(length(find(stddirectedangles>k))/(length(stddirectedangles)))*100];
        diffusiveerror{j} = [diffusiveerror{j},(length(find(stddiffusiveangles<k))/(length(stddiffusiveangles)))*100];
        cutoffangle{j} = [cutoffangle{j}, k];
    end
    OptimizedCutoffAngle = (find(directederror{j} == diffusiveerror{j}));
    
    if isempty(OptimizedCutoffAngle) == 1
        for q = 1 : length(cutoffangle{j})
            errorDiff(q) = abs(directederror{j}(q) - diffusiveerror{j}(q));
        end
        OptimizedCutoffAngle = (find(errorDiff == min(errorDiff)));
    end
    if length(OptimizedCutoffAngle) > 1
        OptimizedCutoffAngle = round(mean(OptimizedCutoffAngle));
    end
    DirectedMissed(j) = directederror{j}(OptimizedCutoffAngle);
    DiffusiveMissed(j) = diffusiveerror{j}(OptimizedCutoffAngle);
    OptimizedCutoffAngle = OptimizedCutoffAngle/10;
    OptimizedCutoffAngle2(j) = OptimizedCutoffAngle;
end
CutoffAngle = mean(OptimizedCutoffAngle2); % Mean at which missed values are minimized
DirectedMissedAtCutoffAngle = mean(DirectedMissed); % Mean directed missed
DiffusiveMissedAtCutoffAngle = mean(DiffusiveMissed); % Mean diffusive missed
close(h)
disp('Generating plots')

cutoffangleplot = []; meandiffusiveerrorplot = []; meandirectederrorplot = [];
semdiffusiveerrorplot = []; semdirectederrorplot = [];

%% Plotting graphs
for kk = 1 : length(cutoffangle{1})
    cutoffangleplot1 = []; diffusiveerrorplot1 = []; directederrorplot1 = [];
    for jj = 1 : length(cutoffangle)
        cutoffangleplot1 = [cutoffangleplot1, cutoffangle{jj}(kk)];
        diffusiveerrorplot1 = [diffusiveerrorplot1, diffusiveerror{jj}(kk)];
        directederrorplot1 = [directederrorplot1, directederror{jj}(kk)];
    end
    cutoffangleplot = [cutoffangleplot, mean(cutoffangleplot1)];
    meandiffusiveerrorplot = [meandiffusiveerrorplot, mean(diffusiveerrorplot1)];
    meandirectederrorplot = [meandirectederrorplot, mean(directederrorplot1)];
    semdiffusiveerrorplot = [semdiffusiveerrorplot, (std(diffusiveerrorplot1))/(sqrt(length(diffusiveerrorplot1)))];
    semdirectederrorplot = [semdirectederrorplot, (std(directederrorplot1))/(sqrt(length(directederrorplot1)))];
end
figure
errorbar(cutoffangleplot, meandiffusiveerrorplot, semdiffusiveerrorplot)
hold on
title('Cutoff angle vs error 2D')
errorbar(cutoffangleplot, meandirectederrorplot, semdirectederrorplot)
hold off
figure
hold on
title('Cutoff angle 2D')
histogram(OptimizedCutoffAngle2);
hold off
figure
hold on
title('Directed missed 2D')
histogram(DirectedMissed);
hold off
figure
hold on
title('Diffusive missed 2D')
histogram(DiffusiveMissed);
hold off
disp('Done!')
end
%%
function [X,Y, displace2, direction2] = MixMSiDirected(DirectedStepSize, TrackingError, SlidingWindowSize,NoofParticlesToSimulate)
%% Simulating directed motion
if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

if nargin == 3
    NoofParticlesToSimulate = 200;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    DirectedAngles = randsample(StaticRandomAngles,1);
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 0; ID{j}(1) = j;
    for k = 2 : 36
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize*cosd(DirectedAngles) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + DirectedStepSize*sind(DirectedAngles) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
end
%% Calculating displacement and angles
% nested function 1
for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand(rise/run);
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end
%% Calculating displacement values and angle values for sliding windows
for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end
%%
function [X,Y, displace2, direction2] = MixMSiDiffusive(DiffusiveStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate)
%% Simulatng diffusive motion
% nested function 2, see nested function 1 for comments
if nargin == 3
    NoofParticlesToSimulate = 200;
end

if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        DiffusiveAngles = randsample(StaticRandomAngles,1);
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DiffusiveStepSize*cosd(DiffusiveAngles) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + DiffusiveStepSize*sind(DiffusiveAngles) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end

end
%%
function  [CutoffAngle,DiffusiveMissedAtCutoffAngle,DirectedMissedAtCutoffAngle] = FindDirectedDiffusiveCutoff1D
%% Find cutoff value for directed and diffusive states for 1-D motion
% for comments see FindDirectedDiffusiveCutoff
NoofParticlesToSimulate = input('Number of particles to simulate ');
NoOfIters = input('Number of times to run the simulation ');
DiffusiveStepSize = input('Diffusive step size in um ');
DirectedStepSize = input('Directed step size in um ');
TrackingError = input('Localization error in um ');
checkwindow = input('Use default window size enter 1, otherwise enter 0 ');
if checkwindow == 1
    SlidingWindowSize = 1;
else
    SlidingWindowSize = input('Size of sliding window ');
end

h = waitbar(0,'0 % Progress...');
for j = 1 : NoOfIters
    waitbar(j/NoOfIters,h,[num2str((j/NoOfIters)*100), ' % progress...'])
    [X,Y, displaceR, directionR] = MixMSiDirected1D(DirectedStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate);
    [X,Y, displaceD, directionD] = MixMSiDiffusive1D(DiffusiveStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate);
    stddirectedangles = []; stddiffusiveangles = [];
    for k = 1 : length(directionR)
        stddirectedangles = [stddirectedangles, directionR{k}];
    end
    for k = 1 : length(directionD)
        stddiffusiveangles = [stddiffusiveangles,directionD{k}];
    end
    
    directederror{j} = []; diffusiveerror{j} = []; cutoffangle{j} = [];
    runfor = 180; %max([max(stddiffusiveangles), max(stddirectedangles)]);
    for k = 0 : 0.1 : runfor
        directederror{j} = [directederror{j},(length(find(stddirectedangles>k))/(length(stddirectedangles)))*100];
        diffusiveerror{j} = [diffusiveerror{j},(length(find(stddiffusiveangles<k))/(length(stddiffusiveangles)))*100];
        cutoffangle{j} = [cutoffangle{j}, k];
    end
    OptimizedCutoffAngle = (find(directederror{j} == diffusiveerror{j}));
    
    if isempty(OptimizedCutoffAngle) == 1
        for q = 1 : length(cutoffangle{j})
            errorDiff(q) = abs(directederror{j}(q) - diffusiveerror{j}(q));
        end
        OptimizedCutoffAngle = (find(errorDiff == min(errorDiff)));
    end
    if length(OptimizedCutoffAngle) > 1
        OptimizedCutoffAngle = round(mean(OptimizedCutoffAngle));
    end
    DirectedMissed(j) = directederror{j}(OptimizedCutoffAngle);
    DiffusiveMissed(j) = diffusiveerror{j}(OptimizedCutoffAngle);
    OptimizedCutoffAngle = OptimizedCutoffAngle/10;
    OptimizedCutoffAngle2(j) = OptimizedCutoffAngle;
end
CutoffAngle = median(OptimizedCutoffAngle2);
DirectedMissedAtCutoffAngle = median(DirectedMissed);
DiffusiveMissedAtCutoffAngle = median(DiffusiveMissed);
close(h)
disp('Generating plots')

cutoffangleplot = []; meandiffusiveerrorplot = []; meandirectederrorplot = [];
semdiffusiveerrorplot = []; semdirectederrorplot = [];

for kk = 1 : length(cutoffangle{1})
    cutoffangleplot1 = []; diffusiveerrorplot1 = []; directederrorplot1 = [];
    for jj = 1 : length(cutoffangle)
        cutoffangleplot1 = [cutoffangleplot1, cutoffangle{jj}(kk)];
        diffusiveerrorplot1 = [diffusiveerrorplot1, diffusiveerror{jj}(kk)];
        directederrorplot1 = [directederrorplot1, directederror{jj}(kk)];
    end
    cutoffangleplot = [cutoffangleplot, mean(cutoffangleplot1)];
    meandiffusiveerrorplot = [meandiffusiveerrorplot, mean(diffusiveerrorplot1)];
    meandirectederrorplot = [meandirectederrorplot, mean(directederrorplot1)];
    semdiffusiveerrorplot = [semdiffusiveerrorplot, (std(diffusiveerrorplot1))/(sqrt(length(diffusiveerrorplot1)))];
    semdirectederrorplot = [semdirectederrorplot, (std(directederrorplot1))/(sqrt(length(directederrorplot1)))];
end
figure
errorbar(cutoffangleplot, meandiffusiveerrorplot, semdiffusiveerrorplot)
hold on
title('Cutoff angle vs error 1D')
errorbar(cutoffangleplot, meandirectederrorplot, semdirectederrorplot)
hold off
figure
hold on
title('Cutoff angle 1D')
histogram(OptimizedCutoffAngle2);
hold off
figure
hold on
title('Directed missed 1D')
histogram(DirectedMissed);
hold off
figure
hold on
title('Diffusive missed 1D')
histogram(DiffusiveMissed);
hold off
disp('Done!')
end
function [X,Y, displace2, direction2] = MixMSiDirected1D(DirectedStepSize, TrackingError, SlidingWindowSize,NoofParticlesToSimulate)

if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

if nargin == 3
    NoofParticlesToSimulate = 200;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    DirectedAngles = randsample(StaticRandomAngles,1);
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 0; ID{j}(1) = j;
    for k = 2 : 36
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DirectedStepSize + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = j;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand(rise/run);
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end

function [X,Y, displace2, direction2] = MixMSiDiffusive1D(DiffusiveStepSize, TrackingError, SlidingWindowSize, NoofParticlesToSimulate)

if nargin == 3
    NoofParticlesToSimulate = 200;
end

if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        DiffusiveAngles = randsample(StaticRandomAngles,1);
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (DiffusiveStepSize*(randsample(([-1,1]),1))) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end

end
%%
function [CutoffDisplacement,DynamicMissedAtCutoffDisplacement,StaticMissedAtCutoffDisplacement] = FindDynamicStaticCutoff
%% Find cutoff value for dynamic and static states for 2-D motion
% for comments see FindDirectedDiffusiveCutoff
NoofParticlesToSimulate = input('Number of particles to simulate ');
NoOfIters = input('Number of times to run the simulation ');
DynamicStepSize = input('Dynamic step size in um ');
TrackingError = input('Localization error in um ');
checkwindow = input('Use default window size enter 1, otherwise enter 0 ');
if checkwindow == 1
    SlidingWindowSize = 1;
else
    SlidingWindowSize = input('Size of sliding window ');
end

h = waitbar(0,'0 % Progress...');

for j = 1 : NoOfIters
    waitbar(j/NoOfIters,h,[num2str((j/NoOfIters)*100), ' % progress...'])
    [XD,YD, displaceD, directionD] = MixMSiDynamic(DynamicStepSize, TrackingError, SlidingWindowSize,...
        NoofParticlesToSimulate);
    [XS,YS, displaceS, directionS] = MixMSiStatic(TrackingError, SlidingWindowSize, NoofParticlesToSimulate);
    maxstaticdisplace = []; maxdynamicdisplace = [];
    for k = 1 : length(displaceD)
        maxdynamicdisplace = [maxdynamicdisplace, displaceD{k}];
    end
    for k = 1 : length(displaceS)
        maxstaticdisplace = [maxstaticdisplace, displaceS{k}];
    end
    
    staticerror{j} = []; dynamicerror{j} = []; cutoffdisplace{j} = [];
    runfor = 2; %max([max(maxstaticdisplace), max(maxdynamicdisplace)]);
    for k = 0 : 0.01 : runfor
        staticerror{j} = [staticerror{j},(length(find(maxstaticdisplace>k))/(length(maxstaticdisplace)))*100];
        dynamicerror{j} = [dynamicerror{j},(length(find(maxdynamicdisplace<=k))/(length(maxdynamicdisplace)))...
            *100];
        cutoffdisplace{j} = [cutoffdisplace{j}, k];
    end
    OptimizedCutoffdisplace = (find(staticerror{j} == dynamicerror{j}));
    
    if isempty(OptimizedCutoffdisplace) == 1
        for q = 1 : length(cutoffdisplace{j})
            errorDiff(q) = abs(staticerror{j}(q) - dynamicerror{j}(q));
        end
        OptimizedCutoffdisplace = (find(errorDiff == min(errorDiff)));
    end
    if length(OptimizedCutoffdisplace) > 1
        OptimizedCutoffdisplace = round(mean(OptimizedCutoffdisplace));
    end
    staticMissed(j) = staticerror{j}(OptimizedCutoffdisplace);
    DynamicMissed(j) = dynamicerror{j}(OptimizedCutoffdisplace);
    OptimizedCutoffdisplace = OptimizedCutoffdisplace/100;
    OptimizedCutoffdisplace2(j) = OptimizedCutoffdisplace;
end
CutoffDisplacement = median(OptimizedCutoffdisplace2);
StaticMissedAtCutoffDisplacement = median(staticMissed);
DynamicMissedAtCutoffDisplacement = median(DynamicMissed);
close(h)
disp('Generating plots')

cutoffdisplaceplot = []; meandynamicerrorplot = []; meanstaticerrorplot = [];
semdynamicerrorplot = []; semstaticerrorplot = [];

for kk = 1 : length(cutoffdisplace{1})
    cutoffdisplaceplot1 = []; dynamicerrorplot1 = []; staticerrorplot1 = [];
    for jj = 1 : length(cutoffdisplace)
        cutoffdisplaceplot1 = [cutoffdisplaceplot1, cutoffdisplace{jj}(kk)];
        dynamicerrorplot1 = [dynamicerrorplot1, dynamicerror{jj}(kk)];
        staticerrorplot1 = [staticerrorplot1, staticerror{jj}(kk)];
    end
    cutoffdisplaceplot = [cutoffdisplaceplot, mean(cutoffdisplaceplot1)];
    meandynamicerrorplot = [meandynamicerrorplot, mean(dynamicerrorplot1)];
    meanstaticerrorplot = [meanstaticerrorplot, mean(staticerrorplot1)];
    semdynamicerrorplot = [semdynamicerrorplot, (std(dynamicerrorplot1))/(sqrt(length(dynamicerrorplot1)))];
    semstaticerrorplot = [semstaticerrorplot, (std(staticerrorplot1))/(sqrt(length(staticerrorplot1)))];
end
figure
errorbar(cutoffdisplaceplot, meandynamicerrorplot, semdynamicerrorplot)
hold on
title('Cutoff displacement vs error 2D')
errorbar(cutoffdisplaceplot, meanstaticerrorplot, semstaticerrorplot)
hold off
figure
hold on
title('Cutoff displacement 2D')
histogram(OptimizedCutoffdisplace2);
hold off
figure
hold on
title('Dynamic missed 2D')
histogram(DynamicMissed);
hold off
figure
hold on
title('Static missed 2D')
histogram(staticMissed);
hold off
disp('Done!')
end

function [X,Y, displace2, direction2] = MixMSiDynamic(DiffusiveStepSize, TrackingError, SlidingWindowSize,...
    NoofParticlesToSimulate)

if nargin == 3
    NoofParticlesToSimulate = 200;
end

if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = round((rand*100) + 1); Y{j}(1) = round((rand*100) + 1); T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        DiffusiveAngles = randsample(StaticRandomAngles,1);
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + DiffusiveStepSize*cosd(DiffusiveAngles) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + DiffusiveStepSize*sind(DiffusiveAngles) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end

function [X,Y, displace2, direction2] = MixMSiStatic(TrackingError, SlidingWindowSize,...
    NoofParticlesToSimulate)

if nargin == 2
    NoofParticlesToSimulate = 200;
end

if nargin == 1
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = round((rand*100) + 1); Y{j}(1) = round((rand*100) + 1); T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end
%%
function [CutoffDisplacement,DynamicMissedAtCutoffDisplacement,StaticMissedAtCutoffDisplacement] = FindDynamicStaticCutoff1D
%% Find cutoff value for dynamic and static states for 1-D motion
% for comments see FindDirectedDiffusiveCutoff
NoofParticlesToSimulate = input('Number of particles to simulate ');
NoOfIters = input('Number of times to run the simulation ');
DynamicStepSize = input('Dynamic step size in um ');
TrackingError = input('Localization error in um ');
checkwindow = input('Use default window size enter 1, otherwise enter 0 ');
if checkwindow == 1
    SlidingWindowSize = 1;
else
    SlidingWindowSize = input('Size of sliding window ');
end
h = waitbar(0,'0 % Progress...');
for j = 1 : NoOfIters
    waitbar(j/NoOfIters,h,[num2str((j/NoOfIters)*100), ' % progress...'])
    [XD,YD, displaceD, directionD] = MixMSiDynamic1D(DynamicStepSize, TrackingError, SlidingWindowSize,...
        NoofParticlesToSimulate);
    [XS,YS, displaceS, directionS] = MixMSiStatic1D(TrackingError, SlidingWindowSize, NoofParticlesToSimulate);
    maxstaticdisplace = []; maxdynamicdisplace = [];
    for k = 1 : length(displaceD)
        maxdynamicdisplace = [maxdynamicdisplace, displaceD{k}];
    end
    for k = 1 : length(displaceS)
        maxstaticdisplace = [maxstaticdisplace, displaceS{k}];
    end
    
    staticerror{j} = []; dynamicerror{j} = []; cutoffdisplace{j} = [];
    runfor = 2; %max([max(maxstaticdisplace), max(maxdynamicdisplace)]);
    for k = 0 : 0.01 : runfor
        staticerror{j} = [staticerror{j},(length(find(maxstaticdisplace>k))/(length(maxstaticdisplace)))*100];
        dynamicerror{j} = [dynamicerror{j},(length(find(maxdynamicdisplace<=k))/(length(maxdynamicdisplace)))...
            *100];
        cutoffdisplace{j} = [cutoffdisplace{j}, k];
    end
    OptimizedCutoffdisplace = (find(staticerror{j} == dynamicerror{j}));
    
    if isempty(OptimizedCutoffdisplace) == 1
        for q = 1 : length(cutoffdisplace{j})
            errorDiff(q) = abs(staticerror{j}(q) - dynamicerror{j}(q));
        end
        OptimizedCutoffdisplace = (find(errorDiff == min(errorDiff)));
    end
    if length(OptimizedCutoffdisplace) > 1
        OptimizedCutoffdisplace = round(mean(OptimizedCutoffdisplace));
    end
    staticMissed(j) = staticerror{j}(OptimizedCutoffdisplace);
    DynamicMissed(j) = dynamicerror{j}(OptimizedCutoffdisplace);
    OptimizedCutoffdisplace = OptimizedCutoffdisplace/100;
    OptimizedCutoffdisplace2(j) = OptimizedCutoffdisplace;
end
CutoffDisplacement = median(OptimizedCutoffdisplace2);
StaticMissedAtCutoffDisplacement = median(staticMissed);
DynamicMissedAtCutoffDisplacement = median(DynamicMissed);
close(h)
disp('Generating plots')

cutoffdisplaceplot = []; meandynamicerrorplot = []; meanstaticerrorplot = [];
semdynamicerrorplot = []; semstaticerrorplot = [];

for kk = 1 : length(cutoffdisplace{1})
    cutoffdisplaceplot1 = []; dynamicerrorplot1 = []; staticerrorplot1 = [];
    for jj = 1 : length(cutoffdisplace)
        cutoffdisplaceplot1 = [cutoffdisplaceplot1, cutoffdisplace{jj}(kk)];
        dynamicerrorplot1 = [dynamicerrorplot1, dynamicerror{jj}(kk)];
        staticerrorplot1 = [staticerrorplot1, staticerror{jj}(kk)];
    end
    cutoffdisplaceplot = [cutoffdisplaceplot, mean(cutoffdisplaceplot1)];
    meandynamicerrorplot = [meandynamicerrorplot, mean(dynamicerrorplot1)];
    meanstaticerrorplot = [meanstaticerrorplot, mean(staticerrorplot1)];
    semdynamicerrorplot = [semdynamicerrorplot, (std(dynamicerrorplot1))/(sqrt(length(dynamicerrorplot1)))];
    semstaticerrorplot = [semstaticerrorplot, (std(staticerrorplot1))/(sqrt(length(staticerrorplot1)))];
end
figure
errorbar(cutoffdisplaceplot, meandynamicerrorplot, semdynamicerrorplot)
hold on
title('Cutoff displacement vs error 1D')
errorbar(cutoffdisplaceplot, meanstaticerrorplot, semstaticerrorplot)
hold off
figure
hold on
title('Cutoff displacement 1D')
histogram(OptimizedCutoffdisplace2);
hold off
figure
hold on
title('Dynamic missed 1D')
histogram(DynamicMissed);
hold off
figure
hold on
title('Static missed 1D')
histogram(staticMissed);
hold off
disp('Done!')

end

function [X,Y, displace2, direction2] = MixMSiDynamic1D(DiffusiveStepSize, TrackingError, SlidingWindowSize,...
    NoofParticlesToSimulate)

if nargin == 3
    NoofParticlesToSimulate = 200;
end

if nargin == 2
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (DiffusiveStepSize*(randsample([1,-1],1))) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end

function [X,Y, displace2, direction2] = MixMSiStatic1D(TrackingError, SlidingWindowSize,...
    NoofParticlesToSimulate)

if nargin == 2
    NoofParticlesToSimulate = 200;
end

if nargin == 1
    NoofParticlesToSimulate = 200;
    SlidingWindowSize = 3;
end

TimeStep = 0.1;
StaticRandomAngles = [0:360];

for j = 1 : NoofParticlesToSimulate
    X{j}(1) = 10; Y{j}(1) = 10; T{j}(1) = 1; ID{j}(1) = 1;
    for k = 2 : 36
        randomdirection = randsample(StaticRandomAngles,1);
        X{j}(k) = X{j}(k-1) + (TrackingError*cosd(randomdirection));
        Y{j}(k) = Y{j}(k-1) + (TrackingError*sind(randomdirection));
        T{j}(k) = T{j}(k-1) + TimeStep;
        ID{j}(k) = 1;
    end
end

for k = 1 : 200
    for j = 1 : 35
        displace{k}(j) = sqrt(((X{k}(j+1) - X{k}(j))^2) + ((Y{k}(j+1) - Y{k}(j))^2));
        rise = Y{k}(j + 1) - Y{k}(j); run = X{k}(j + 1) - X{k}(j);
        if rise > 0 && run > 0
            angle1 = atand(rise/run);
        elseif rise > 0 && run < 0
            angle1 = 180 - abs(atand(rise/run));
        elseif rise < 0 && run < 0
            angle1 = 180 + abs(atand(rise/run));
        elseif rise < 0 && run > 0
            angle1 = 360 - abs(atand(rise/run));
        elseif (rise == 0 && run > 0) || (rise > 0 && run == 0)
            angle1 = atand((rise/run));
        elseif (rise == 0 && run < 0)
            angle1 = 180;
        elseif (rise < 0 && run == 0)
            angle1 = 270;
        end
        if angle1 > 180
            angle1 = 360 - angle1;
        end
        direction{k}(j) = angle1;
    end
end

for k = 1 : length(displace)
    runfor = length(displace{k}) - SlidingWindowSize;
    for j = 1 : runfor
        displace2{k}(j) = max(displace{k}(j:j+SlidingWindowSize));
        direction2{k}(j) = std(direction{k}(j:j+SlidingWindowSize));
    end
end
end
%%
% The end!!