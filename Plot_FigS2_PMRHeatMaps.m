%% Code for Figure 2CD - Heatmaps of PMR across synchrony and sampling
clear all; close all; clc;

%% Determine PMR across synchrony and sampling
 
% shift of beta distribution peak
ICshiftvec = linspace(0,1,9);
% level of initial synchrony
ICshapevec = fliplr([1 10 50 100 500]);

% Sequestered mortality rate
museq = 0;
% Non-sequestered mortality rate
mu = 0;
% Cycle length (hrs)
cycle_len = 48;

% PMR
R = 6;

maxPMRobs = NaN(length(ICshiftvec)-1,length(ICshapevec));
PMRreg = NaN(length(ICshiftvec)-1,length(ICshapevec));

for j1 = 1:length(ICshiftvec)-1
    
    for j2 = 1:length(ICshapevec)
        
        ICshift = mod(0.5+ICshiftvec(j1),1);
        ICshape = ICshapevec(j2);
        
        % Run modeling code
        [t, tot, totSeq, totNonSeq, p] = RunRBC(mu,museq,cycle_len,ICshift,ICshape,R);
        
        % Determine observed PMR by day
        PMRobs = NaN(length(totNonSeq)-480,1);
        for j = 1:length(totNonSeq)-480
            PMRobs(j) = totNonSeq(j+480)/totNonSeq(j);
        end
        
        % Find time points of first 8 days
        for days = 1:8
            td(days) = find(t/24==days);
        end
        
        % Log-linear regression including data from first 7 days
        t1 = t(td(1:7))/24;
        y1 = log10(totNonSeq(td(1:7),:));
        pfit = polyfit(t1,y1,1);
        a = pfit(2);
        b = pfit(1);
        
        PMRreg(j1,j2) = 10^(2*b);
        maxPMRobs(j1,j2) = max(PMRobs(td(1:5)));
        
    end
    
end

%%  Determine the bursting duration (in hrs)
sp_vec = fliplr([1 10 50 100 500]);
for sp1 = 1:length(sp_vec)
    sp = sp_vec(sp1);
    y = betainv([0.005 .995],sp,sp)*48;
    Burst99(sp1) = y(2) - y(1);
end
B99 = round(Burst99);

%% Plot Fig S2A - max PMR_obs across synchrony and sampling

f1 = figure(1);


h = heatmap(ICshiftvec(1:end-1)*48,B99,maxPMRobs');
h.CellLabelFormat = '%.0f';
h.NodeChildren(3).YDir='normal';

set(gca,'fontsize',24)
ylabel('Bursting duration (hrs)')
xlabel({'Initial median parasite age'})
caxis([0 150])
h.FontName = 'Arial';
x = h.Colormap;
h.Colormap = [x(:,2) x(:,[1 1])];


annotation('textarrow',[.96,.96],[.97,.97],'string','PMR', ...
    'HeadStyle','none','LineStyle','none',...
    'HorizontalAlignment','center','Fontsize',20);


title('(A) Max. observed PMR')
set(struct(h).Axes.Title,'FontWeight','normal','fontsize',24,'position',[3 5.5])


%% Plot Fig S2B - PMR_reg across synchrony and sampling

f2 = figure(2);
%
h = heatmap(ICshiftvec(1:end-1)*48,B99,PMRreg');
h.CellLabelFormat = '%.0f';
h.NodeChildren(3).YDir='normal';


set(gca,'fontsize',24)
ylabel('Bursting duration (hrs)')
xlabel({'Initial median parasite age'})
caxis([0 30])
h.FontName = 'Arial';
annotation('textarrow',[.95,.95],[.97,.97],'string','PMR', ...
    'HeadStyle','none','LineStyle','none',...
    'HorizontalAlignment','center','Fontsize',20)

title('(B) PMR from regression')
set(struct(h).Axes.Title,'FontWeight','normal','fontsize',24,'position',[3 5.5])
