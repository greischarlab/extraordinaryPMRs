%% Plot  Figure 2C - PMR variation across sampling different initial median parasite age
% Asynchronous infection
clear all; close all; clc;

%% Determine PMR across sampling times

% shift of beta distribution peak
ICshiftvec = linspace(0,1,9*8+1);
% level of initial synchrony (asynchronous)
ICshape = 1;

% Sequestered mortality rate
museq = 0;
% Non-sequestered mortality rate
mu = 0;
% Cycle length (hrs)
cycle_len = 48;

maxPMRobs = NaN(length(ICshiftvec),1);
PMRreg = NaN(length(ICshiftvec),1);

% PMR
R = 6;


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
    PMRreg(j1,1) = 10^(2*b);
    % Max. PMR observed
    maxPMRobs(j1,1) = max(PMRobs(td(1:5)));
    
end

%% Plot Fig 2B - PMR_obs and PMR_reg across sampling times

grayColor = [.7 .7 .7];

f1 = figure();

a = area([0 1],log2([32 32]));
a.FaceColor = grayColor;
a.FaceAlpha = 0.5;
a.EdgeAlpha = 0;

hold on
a1  = area([0 1],log2([16 16]));
a1.FaceColor = grayColor;
a1.EdgeAlpha = 0;

% Maximum PMR_obs
plot(ICshiftvec,log2(maxPMRobs(:,1)),'r-','markersize',20);
% PMR_reg
plot(ICshiftvec,log2(PMRreg(:,1)),'b-');

% Determine and plot sampling times from Fig 2C,D
ICshifttmp = linspace(0,1,9);
for j = 1:length(ICshifttmp)
    id(j) = find(ICshiftvec==ICshifttmp(j));
end

% True PMR
p3 = plot([0 1],log2([p.R p.R]),'k-','linewidth',2);

p2 = plot(ICshiftvec(id),log2(PMRreg(id,1)),'bd','markersize',10,'markerfacecolor','w');
p1 = plot(ICshiftvec(id),log2(maxPMRobs(id,1)),'r.','markersize',20);

ylim([0 log2(512)])
set(gca,'ytick',0:8,'yticklabel',2.^(0:8))
l = legend([p1(1) p2(1) p3],'Max. observed','From regression','True PMR');
set(l,'box','off','location','southeast','fontsize',20)
set(gca,'fontname','Arial','fontsize',24)
ylabel('PMR')
set(gca,'xtick',linspace(0,1,9),'xticklabel',linspace(0,48,9))
xlabel('Initial median parasite age')
set(gca,'box','off')
set(gca,'TickDir','out','box','off','fontname','Arial')
text(0.01,9.5,'(C) Simulated asynchronous infection','fontname','Arial','fontsize',24)

