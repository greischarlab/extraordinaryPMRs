%% Plot  Figure S4 - Vary True PMR

clear all; close all; clc;

%% Determine PMR across sampling times

% shift of beta distribution peak
ICshiftvec = linspace(0,1,49);
% level of initial synchrony (~9 hrs for 99% bursting)
ICshape = 100;

% Sequestered mortality rate
museq = 0;
% Non-sequestered mortality rate
mu = 0;
% Cycle length (hrs)
cycle_len = 48;

% True PMRs considered
PMRvec = [0.1 0.25 0.5 0.75 0.9];


maxPMRobs = NaN(length(ICshiftvec),length(PMRvec));
PMRreg = NaN(length(ICshiftvec),length(PMRvec));
tic
for r1 = 1:length(PMRvec)
    R = PMRvec(r1);

for j1 = 1:length(ICshiftvec)-1
    
    ICshift = mod(0.5+ICshiftvec(j1),1);
    
    % Run modeling code
    [t, tot, totSeq, totNonSeq, p] = RunRBC(mu,museq,cycle_len,ICshift,ICshape,R);
    
    % Determine observed PMR by day
    PMRobs = NaN(length(totNonSeq)-480,1);
    for j = 1:length(totNonSeq)-480
        PMRobs(j) = totNonSeq(j+480)/totNonSeq(j);
    end
    
    % Find time points of first 8 days
    td = NaN(1,8);
    for days = 1:8
        td(days) = find(t/24==days);
    end
    
    % Log-linear regression including data from first 7 days
    t1 = t(td(1:7))/24;
    y1 = log10(totNonSeq(td(1:7),:));
    pfit = polyfit(t1,y1,1);
    a = pfit(2);
    b = pfit(1);
    
    % PMR_reg, from log-linear regression
    PMRreg(j1,r1) = 10^(2*b);
    % Max. PMR observed
    maxPMRobs(j1,r1) = max(PMRobs(td(1:5)));
    
end
toc
end

%% Determine max estimated PMRs for each true PMR
maxmaxPMR = NaN(length(PMRvec),1);
maxmaxPMRreg = NaN(length(PMRvec),1);

for r1 = 1:length(PMRvec)
    maxmaxPMR(r1) = max(maxPMRobs(:,r1));
    maxmaxPMRreg(r1) = max(PMRreg(:,r1));
end


%% Figure S4A

f1 = figure(1);
h = heatmap(ICshiftvec(1:6:end-1)*48,PMRvec,maxPMRobs(1:6:end-1,:)');
h.CellLabelFormat = '%.1f';
h.NodeChildren(3).YDir='normal';

set(gca,'fontsize',24)
ylabel('True PMR')
xlabel({'Initial median parasite age'})
caxis([0 5])
h.FontName = 'Arial';
x = h.Colormap;
h.Colormap = [x(:,2) x(:,[1 1])];

annotation('textarrow',[.96,.96],[.97,.97],'string','PMR', ...
    'HeadStyle','none','LineStyle','none',...
    'HorizontalAlignment','center','Fontsize',20);

title('(A) Max. observed PMR')
set(struct(h).Axes.Title,'FontWeight','normal','fontsize',24,'position',[3 5.5])

%% Figure S4B

f2 = figure(2);
h = heatmap(ICshiftvec(1:6:end-1)*48,PMRvec,PMRreg(1:6:end-1,:)');
h.CellLabelFormat = '%.1f';
h.NodeChildren(3).YDir='normal';

set(gca,'fontsize',24)
ylabel('True PMR')
xlabel({'Initial median parasite age'})
caxis([0 5])
h.FontName = 'Arial';


annotation('textarrow',[.96,.96],[.97,.97],'string','PMR', ...
    'HeadStyle','none','LineStyle','none',...
    'HorizontalAlignment','center','Fontsize',20);

title('(B) PMR from regression')
set(struct(h).Axes.Title,'FontWeight','normal','fontsize',24,'position',[3 5.5])
