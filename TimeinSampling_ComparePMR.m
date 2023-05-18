%% Plot  Figure 6 -  PMR variation across sampling times
clear all; close all; clc;

%% Determine PMR across sampling times

% PMR
PMRvec = [6 2:2:32];


% shift of beta distribution peak
ICshiftvec = linspace(0,1,9);
ICshiftvec(end) = [];
% level of initial synchrony (~9 hrs for 99% bursting)
ICshape = 100;

% Sequestered mortality rate
museq = 0;
% Non-sequestered mortality rate
mu = 0;
% Cycle length (hrs)
cycle_len = 48;

winvec = [6:6:360]/6; % minutes shifted

maxPMRobs = NaN(length(ICshiftvec),length(PMRvec),length(winvec));
maxPMRobs_early = NaN(length(ICshiftvec),length(PMRvec),length(winvec));
maxPMRobs_late = NaN(length(ICshiftvec),length(PMRvec),length(winvec));
PMRreg = NaN(length(ICshiftvec),length(PMRvec),length(winvec));
PMRreg_early = NaN(length(ICshiftvec),length(PMRvec),length(winvec));
PMRreg_late = NaN(length(ICshiftvec),length(PMRvec),length(winvec));

tic
for r1 = 1:length(PMRvec)
    R = PMRvec(r1);
    for j1 = 1:length(ICshiftvec)        
        
        ICshift = mod(0.5+ICshiftvec(j1),1);
        
        % Run modeling code
        [t, tot, totSeq, totNonSeq, p] = RunRBC(mu,museq,cycle_len,ICshift,ICshape,R);
        
        % Determine observed PMR by day
        PMRobs = NaN(length(totNonSeq)-480,1);
        for j = 1:length(totNonSeq)-480
            PMRobs(j) = totNonSeq(j+480)/totNonSeq(j);
        end
        
        for w1 = 1:length(winvec)
            window = winvec(w1);
            
            % Find time points of first 8 days
            td = NaN(1,8);
            td_early = NaN(1,8);
            td_late = NaN(1,8);
            for days = 1:8
                td(days) = find(t/24==days);
                td_early(days) = td(days)-window;
                td_late(days) = td(days)+window;
            end
            
            % Log-linear regression including data from first 7 days
            t1 = t(td(1:7))/24;
            y1 = log10(totNonSeq(td(1:7),:));
            pfit = polyfit(t1,y1,1);
            b = pfit(1);
            
            t1_early = t(td_early(1:7))/24;
            y1_early = log10(totNonSeq(td_early(1:7),:));
            pfit_early = polyfit(t1_early,y1_early,1);
            b_early = pfit_early(1);
            
            t1_late = t(td_late(1:7))/24;
            y1_late = log10(totNonSeq(td_late(1:7),:));
            pfit_late = polyfit(t1_late,y1_late,1);
            b_late = pfit_late(1);
            
            % PMR_reg, from log-linear regression
            PMRreg(j1,r1,w1) = 10^(2*b);
            PMRreg_early(j1,r1,w1) = 10^(2*b_early);
            PMRreg_late(j1,r1,w1) = 10^(2*b_late);
            
            % Max. PMR observed
            maxPMRobs(j1,r1,w1) = max(PMRobs(td(1:5)));
            maxPMRobs_early(j1,r1,w1) = max(PMRobs(td_early(1:5)));
            maxPMRobs_late(j1,r1,w1) = max(PMRobs(td_late(1:5)));
            
        end
    end
    toc
    save('Output_ComparePMR.mat')
    
end
%% Output Data

clear; close all; clc

load('Output_ComparePMR.mat')

%% PMR from Regression - Determine window of comparison

% From PMRreg and maxPMR
% dim1 = shift
% dim2 = PMR
% dim3 = window

shifts = length(ICshiftvec);
windows = length(winvec);

pmin = NaN(shifts,windows);
pmax = NaN(shifts,windows);
pmin1 = NaN(shifts,length(PMRvec),windows);
pmax1 = NaN(shifts,length(PMRvec),windows);


for shift = 1:shifts
    for win = 1:windows
        early_min = min(squeeze(PMRreg_early(shift,1,1:win)));
        early_max = max(squeeze(PMRreg_early(shift,1,1:win)));
        late_min = min(squeeze(PMRreg_late(shift,1,1:win)));
        late_max = max(squeeze(PMRreg_late(shift,1,1:win)));
        pmin(shift,win) = min([early_min early_max late_min late_max]);
        pmax(shift,win) = max([early_min early_max late_min late_max]);
        
        for pmr = 2:length(PMRvec)
            early_min1 = min(squeeze(PMRreg_early(shift,pmr,1:win)));
            early_max1 = max(squeeze(PMRreg_early(shift,pmr,1:win)));
            late_min1 = min(squeeze(PMRreg_late(shift,pmr,1:win)));
            late_max1 = max(squeeze(PMRreg_late(shift,pmr,1:win)));
            
            pmin1(shift,pmr,win) = min([early_min1 early_max1 late_min1 late_max1]);
            pmax1(shift,pmr,win) = max([early_min1 early_max1 late_min1 late_max1]);
        end
    end
end

out = NaN(shifts,length(PMRvec));
for j = 2:length(PMRvec)
    x = squeeze(pmax1(:,j,:));
    y = squeeze(pmin1(:,j,:));
    z = ((x-pmin)>0 & (y-pmax)<0);
    for shift = 1:shifts
        id = find(z(shift,:)==1,1,'first');
        if ~isempty(id)
            out(shift,j) = winvec(id); % first index with overlap
        end
    end
end

%% PMR from regression - last window with no overlap (can be distinguished)
out2 = out-1; % last index with no overlap
out2(out2>(360/12)) = NaN; % if index is above 6 hour window
out3reg = fliplr(out2(1:end,2:end))'*12/60; % changes to hours

%%  PMR from regression (compare true PMR of 6 and 12) - 12 minute window

minPMRreg1 = NaN(1,size(PMRreg_early,1));
maxPMRreg1 = NaN(1,size(PMRreg_early,1));
minPMRreg7 = NaN(1,size(PMRreg_early,1));
maxPMRreg7 = NaN(1,size(PMRreg_early,1));

for j = 1:size(PMRreg_early,1)
    minPMRreg1(j) = min(min(PMRreg_early(j,1,1)),min(PMRreg_late(j,1,1)));
    maxPMRreg1(j) = max(max(PMRreg_early(j,1,1)),max(PMRreg_late(j,1,1)));
    minPMRreg7(j) = min(min(PMRreg_early(j,7,1)),min(PMRreg_late(j,7,1)));
    maxPMRreg7(j) = max(max(PMRreg_early(j,7,1)),max(PMRreg_late(j,7,1)));    
end

outMaxPMRreg12m = [minPMRreg1(1:end)' maxPMRreg1(1:end)' minPMRreg7(1:end)' maxPMRreg7(1:end)'];

%%  PMR from regression (compare true PMR of 6 and 12) - 6 hr window

minPMRreg1 = NaN(1,size(PMRreg_early,1));
maxPMRreg1 = NaN(1,size(PMRreg_early,1));
minPMRreg7 = NaN(1,size(PMRreg_early,1));
maxPMRreg7 = NaN(1,size(PMRreg_early,1));

for j = 1:size(PMRreg_early,1)
    minPMRreg1(j) = min(min(PMRreg_early(j,1,1:30)),min(PMRreg_late(j,1,1:30)));
    maxPMRreg1(j) = max(max(PMRreg_early(j,1,1:30)),max(PMRreg_late(j,1,1:30)));
    minPMRreg7(j) = min(min(PMRreg_early(j,7,1:30)),min(PMRreg_late(j,7,1:30)));
    maxPMRreg7(j) = max(max(PMRreg_early(j,7,1:30)),max(PMRreg_late(j,7,1:30)));
end

outMaxPMRreg6h = [minPMRreg1(1:end)' maxPMRreg1(1:end)' minPMRreg7(1:end)' maxPMRreg7(1:end)'];


%% Max observed PMR - Determine window of comparison

% dim1 = shift
% dim2 = PMR
% dim3 = window

pmin = NaN(shifts,windows);
pmax = NaN(shifts,windows);
pmin1 = NaN(shifts,length(PMRvec),windows);
pmax1 = NaN(shifts,length(PMRvec),windows);


for shift = 1:shifts
    for win = 1:windows
        early_min = min(squeeze(maxPMRobs_early(shift,1,1:win)));
        early_max = max(squeeze(maxPMRobs_early(shift,1,1:win)));
        late_min = min(squeeze(maxPMRobs_late(shift,1,1:win)));
        late_max = max(squeeze(maxPMRobs_late(shift,1,1:win)));
        pmin(shift,win) = min([early_min early_max late_min late_max]);
        pmax(shift,win) = max([early_min early_max late_min late_max]);
        
        for pmr = 2:length(PMRvec)
            early_min1 = min(squeeze(maxPMRobs_early(shift,pmr,1:win)));
            early_max1 = max(squeeze(maxPMRobs_early(shift,pmr,1:win)));
            late_min1 = min(squeeze(maxPMRobs_late(shift,pmr,1:win)));
            late_max1 = max(squeeze(maxPMRobs_late(shift,pmr,1:win)));
            
            pmin1(shift,pmr,win) = min([early_min1 early_max1 late_min1 late_max1]);
            pmax1(shift,pmr,win) = max([early_min1 early_max1 late_min1 late_max1]);
        end
    end
end

out = NaN(shifts,length(PMRvec));
for j = 2:length(PMRvec)
    x = squeeze(pmax1(:,j,:));
    y = squeeze(pmin1(:,j,:));
    z = ((x-pmin)>0 & (y-pmax)<0);
    for shift = 1:shifts
        id = find(z(shift,:)==1,1,'first');
        if ~isempty(id)
            out(shift,j) = winvec(id); % first index with overlap
        end
    end
end

%% Max observed PMR - last window with no overlap (can be distinguished)
out2 = out-1; % last index with no overlap
out2(out2>(360/12)) = NaN; % if index is above 6 hour window
out3 = fliplr(out2(1:end,2:end))'*12/60; % change to hours

%% Max observed PMR (compare true PMR of 6 and 12) - 12 minute window

minPMR1 = NaN(1,size(PMRreg_early,1));
maxPMR1 = NaN(1,size(PMRreg_early,1));
minPMR7 = NaN(1,size(PMRreg_early,1));
maxPMR7 = NaN(1,size(PMRreg_early,1));

for j = 1:size(maxPMRobs_early,1)
    minPMR1(j) = min(min(maxPMRobs_early(j,1,1)),min(maxPMRobs_late(j,1,1)));
    maxPMR1(j) = max(max(maxPMRobs_early(j,1,1)),max(maxPMRobs_late(j,1,1)));
    minPMR7(j) = min(min(maxPMRobs_early(j,7,1)),min(maxPMRobs_late(j,7,1)));
    maxPMR7(j) = max(max(maxPMRobs_early(j,7,1)),max(maxPMRobs_late(j,7,1)));
end
outMaxPMR12m = [minPMR1(1:end)' maxPMR1(1:end)' minPMR7(1:end)' maxPMR7(1:end)'];

%% Max observed PMR (compare true PMR of 6 and 12) - 6 hr window

minPMR1 = NaN(1,size(PMRreg_early,1));
maxPMR1 = NaN(1,size(PMRreg_early,1));
minPMR7 = NaN(1,size(PMRreg_early,1));
maxPMR7 = NaN(1,size(PMRreg_early,1));

for j = 1:size(maxPMRobs_early,1)
    minPMR1(j) = min(min(maxPMRobs_early(j,1,1:30)),min(maxPMRobs_late(j,1,1:30)));
    maxPMR1(j) = max(max(maxPMRobs_early(j,1,1:30)),max(maxPMRobs_late(j,1,1:30)));
    minPMR7(j) = min(min(maxPMRobs_early(j,7,1:30)),min(maxPMRobs_late(j,7,1:30)));
    maxPMR7(j) = max(max(maxPMRobs_early(j,7,1:30)),max(maxPMRobs_late(j,7,1:30)));    
end

outMaxPMR6h = [minPMR1(1:end)' maxPMR1(1:end)' minPMR7(1:end)' maxPMR7(1:end)'];


%% Output data to spreadsheets

%
writematrix(out3,'PMRHeatmapMaxObs.csv')
writematrix(out3reg,'PMRHeatmapFromReg.csv')
writematrix(outMaxPMRreg6h,'PMRFromReg6h.csv')
writematrix(outMaxPMRreg12m,'PMRFromReg12m.csv')
writematrix(outMaxPMR6h,'PMRMaxObs6h.csv')
writematrix(outMaxPMR12m,'PMRMaxObs12m.csv')
%}
