%% Code for Figure S5 - Loss of synchrony
clear all; close all; clc;

%% Simulate time course

% shift of beta distribution peak (~34 hrs)
ICshiift = 0.34;
ICshift = mod(0.5+ICshiift,1);
% level of initial synchrony (~9 hrs for 99% bursting)
ICshape = 100;

% Sequestered mortality rate
museq = 0;
% Non-sequestered mortality rate
mu = 0;
% Cycle length (hrs)
cycle_len = 48;

% PMR
R = 6;

% Run modeling code
[t, tot, totSeq, totNonSeq, p] = RunRBCLong(mu,museq,cycle_len,ICshift,ICshape,R);

% Determine observed PMR by day
PMRobs = NaN(length(totNonSeq)-480,1);
for j = 1:length(totNonSeq)-480
    PMRobs(j) = totNonSeq(j+480)/totNonSeq(j);
end

% Find time points of first 8 days
td = NaN(1,8);
for days = 1:25
    td(days) = find(t/24==days);
end

% Log-linear regression including data from first 7 days
t1 = t(td(1:7))/24;
y1 = log10(totNonSeq(td(1:7),:));
pfit = polyfit(t1,y1,1);
a = pfit(2);
b = pfit(1);

% PMR_reg, from log-linear regression
PMRreg = 10^(2*b);
% Max. PMR observed
maxPMRobs = max(PMRobs(td(1:end)));

% Display PMR_reg and max PMR_obs
fprintf('\n PMR_reg = %0.2f \n PMR_obs = %0.2f \n',PMRreg,maxPMRobs)

%%

PMRregALL = NaN(25,1);
PMRreg5 = NaN(25,1);

for j = 5:25
    clear t1 y1 pfit a b
    % Log-linear regression 
    t1 = t(td(1:j))/24;
    y1 = log10(totNonSeq(td(1:j),:)); % include all data up to time point
    pfit = polyfit(t1,y1,1);
    a = pfit(2);
    b = pfit(1);
    
    % PMR_reg, from log-linear regression
    PMRregALL(j-4) = 10^(2*b);
end

for j = 7:25
    clear t1 y1 pfit a b
    % Log-linear regression 
    t1 = t(td(j-6:j))/24;
    y1 = log10(totNonSeq(td(j-6:j),:)); % include data from previous five days
    pfit = polyfit(t1,y1,1);
    a = pfit(2);
    b = pfit(1);
    
    % PMR_reg, from log-linear regression
    PMRreg5(j-6) = 10^(2*b);
end

%% Plot PMR measured - 5 day window

maxPMRobstime = NaN(length(td),1);
for j = 1:length(td)-4
    maxPMRobstime(j) = max(PMRobs(td(j):td(j+4)));
end

%% Figure S5A
figure(1)

l3 = semilogy(0:24,R*ones(1,25),'k-','linewidth',2);
hold on
l1 = semilogy(t(td(1:end))/24,maxPMRobstime,'r.-','markersize',24,'linewidth',1);
l2 = semilogy(1:25,PMRreg5,'bd-','markersize',10,'linewidth',1,'markerfacecolor','w');
l4 = semilogy(t(td(1:end))/24,PMRobs(td(1:end)),'x-','color',[54.5/100, 0, 0],'linewidth',1);
semilogy(t(td(1:end))/24,PMRobs(td(1:end)),'x','color',[54.5/100, 0, 0],'linewidth',2);

hold off
xlim([0 21])
ylim([5 100])
set(gca,'fontname','Arial','fontsize',24)
xlabel('Days')
ylabel('Estimated PMR')
l = legend([l1 l4 l2 l3],{'Max. obs. PMR','Obs. PMR','From regression','True PMR'},'location','northeast');
set(l,'box','off','fontsize',20)
set(gca,'TickDir','out','box','off','fontname','Arial')
set(gca,'ytick',[5 10 100])
text(0.3,120,'(A) Estimated PMR through time','fontsize',24)

%% Figure S5B - measurable fraction

figure(2)
y = (tot - totNonSeq)./tot*100;

plot(t/24,y,'linewidth',2)
set(gca,'fontsize',24)
xlim([0 21])
xlabel('Days')
ylabel({'Observed % of true abundance'})
set(gca,'TickDir','out','box','off','fontname','Arial')
text(0.3,105,'(B) Measurable fraction','fontsize',24)
