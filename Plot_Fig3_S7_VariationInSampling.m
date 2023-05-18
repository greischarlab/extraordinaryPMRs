%% Plot  Figure 3 (and S7) -  PMR variation across sampling times
% Initial median parasite age = 18 h
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

maxPMRobs = NaN(length(ICshiftvec),1);
maxPMRobs_early = NaN(length(ICshiftvec),1);
maxPMRobs_late = NaN(length(ICshiftvec),1);
PMRreg = NaN(length(ICshiftvec),1);
PMRreg_early = NaN(length(ICshiftvec),1);
PMRreg_late = NaN(length(ICshiftvec),1);

% PMR
R = 6;

f1 = figure(1);
f2 = figure(2);
set(f2,'position',[10 10 560*2 420])


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
    td_early12 = NaN(1,8);
    td_late12 = NaN(1,8);
    td_early = NaN(1,8);
    td_late = NaN(1,8);
    td_early6 = NaN(1,8);
    td_late6 = NaN(1,8);
    for days = 1:8
        td(days) = find(t/24==days);
        td_early12(days) = find(t>(days*24-.1),1,'first');
        td_late12(days) = find(t>(days*24+.1),1,'first');
        td_early(days) = find(t==(days*24-.5));
        td_late(days) = find(t==(days*24+.5));
        td_early6(days) = find(t==(days*24-3));
        td_late6(days) = find(t==(days*24+3));
    end
    
    % Log-linear regression including data from first 7 days
    t1 = t(td(1:7))/24;
    y1 = log10(totNonSeq(td(1:7),:));
    pfit = polyfit(t1,y1,1);
    b = pfit(1);
    a = pfit(2);
    
    
    % Standard time - 6 minutes
    t1_early12 = t(td_early12(1:7))/24;
    y1_early12 = log10(totNonSeq(td_early12(1:7),:));
    pfit_early12 = polyfit(t1_early12,y1_early12,1);
    b_early12 = pfit_early12(1);
    a_early12 = pfit_early12(2);
    
    % Standard time + 6 minutes
    t1_late12 = t(td_late12(1:7))/24;
    y1_late12 = log10(totNonSeq(td_late12(1:7),:));
    pfit_late12 = polyfit(t1_late12,y1_late12,1);
    b_late12 = pfit_late12(1);
    a_late12 = pfit_late12(2);
    
    
    % Standard time - 30 minutes
    t1_early = t(td_early(1:7))/24;
    y1_early = log10(totNonSeq(td_early(1:7),:));
    pfit_early = polyfit(t1_early,y1_early,1);
    b_early = pfit_early(1);
    a_early = pfit_early(2);
    
    % Standard time + 30 minutes
    t1_late = t(td_late(1:7))/24;
    y1_late = log10(totNonSeq(td_late(1:7),:));
    pfit_late = polyfit(t1_late,y1_late,1);
    b_late = pfit_late(1);
    a_late = pfit_late(2);
    
    
    % Standard time - 3 hours
    t1_early6 = t(td_early6(1:7))/24;
    y1_early6 = log10(totNonSeq(td_early6(1:7),:));
    pfit_early6 = polyfit(t1_early6,y1_early6,1);
    b_early6 = pfit_early6(1);
    a_early6 = pfit_early6(2);
    
    % Standard time + 3 hours
    t1_late6 = t(td_late6(1:7))/24;
    y1_late6 = log10(totNonSeq(td_late6(1:7),:));
    pfit_late6 = polyfit(t1_late6,y1_late6,1);
    b_late6 = pfit_late6(1);
    a_late6 = pfit_late6(2);
    
    % PMR_reg, from log-linear regression
    PMRreg(j1,1) = 10^(2*b);
    PMRreg_early(j1,1) = 10^(2*b_early);
    PMRreg_late(j1,1) = 10^(2*b_late);
    % Max. PMR observed
    maxPMRobs(j1,1) = max(PMRobs(td(1:5)));
    maxPMRobs_early(j1,1) = max(PMRobs(td_early(1:5)));
    maxPMRobs_late(j1,1) = max(PMRobs(td_late(1:5)));
    
    if ICshiftvec(j1)*48==18

        % Figure S7
        figure(2)
        subplot(1,2,1)
        grayColor = [.7 .7 .7];
        
        PMRearly = (2*b_early6);
        PMRlate = (2*b_late6);
        
        px(1) = plot([0 7],PMRearly*ones(2,1),'-','linewidth',3,'color','b');%,[0, 0.4470, 0.7410]);
        hold on
        px(2) = plot([0 7],PMRlate*ones(2,1),'--','linewidth',3,'color','b');%,[0, 0.4470, 0.7410]);
        px(5) = plot([0 7],log10([R R]),'k-','linewidth',2);
        for c = 0:6
        a1 = area([1-3/24 1-0.1/24]+c,[3 3]);
        a2 = area([1+3/24 1+0.1/24]+c,[3 3]);
        a1.FaceColor = grayColor;
        a1.FaceAlpha = 0.5;
        a1.EdgeAlpha = 0;
        a2.FaceColor = grayColor;
        a2.FaceAlpha = 0.5;
        a2.EdgeAlpha = 0;
        end
        
        plot(t(10:end-480)/24,log10(PMRobs(10:end)),'linewidth',1,'color','r');
        px(3) = plot(t(211:271)/24,log10(PMRobs(211:271)),'linewidth',3,'color','r');
        plot(t(451:511)/24,log10(PMRobs(451:511)),'linewidth',3,'color','r');
        plot(t(691:751)/24,log10(PMRobs(691:751)),'linewidth',3,'color','r');
        plot(t(931:991)/24,log10(PMRobs(931:991)),'linewidth',3,'color','r');
        hold off
        ylim([.6 2.6])
        xlim([0.5 3.5])
        set(gca,'ytick',0:3,'yticklabel',{'1','10','10^2','10^3'})
        set(gca,'xtick',0:4)
        set(gca,'fontsize',20)
        xlabel('Day of blood-stage infection')
        ylabel('PMR')      
        set(gca,'TickDir','out','box','off','fontname','Arial')
        text(0.25,2.7,'(A) Inference with 6 hr sampling window ','fontsize',24)
        
        subplot(1,2,2)
        grayColor = [.7 .7 .7];
        
        PMRearly12 = (2*b_early12);
        PMRlate12 = (2*b_late12);
        
        px(1) = plot([0 7],PMRearly12*ones(2,1),'-','linewidth',3,'color','b');
        hold on
        px(2) = plot([0 7],PMRlate12*ones(2,1),'--','linewidth',3,'color','b');
        px(5) = plot([0 7],log10([R R]),'k-','linewidth',2);
        for c = 0:6
        a = area([1-0.1/24 1+0.1/24]+c,[3 3]);
        a.FaceColor = grayColor;
        a.FaceAlpha = 0.5;
        a.EdgeAlpha = 0;
        end
        
        plot(t(10:end-480)/24,log10(PMRobs(10:end)),'linewidth',1,'color','r');
        px(3) = plot(t(240:242)/24,log10(PMRobs(240:242)),'linewidth',3,'color','r');
        plot(t(480:482)/24,log10(PMRobs(480:482)),'linewidth',3,'color','r');
        plot(t(720:722)/24,log10(PMRobs(720:722)),'linewidth',3,'color','r');
        hold off
        ylim([.6 2.6])
        xlim([0.5 3.5])
        set(gca,'ytick',0:3,'yticklabel',{'1','10','10^2','10^3'})
        set(gca,'xtick',0:4)
        set(gca,'fontsize',20)
        xlabel('Day of blood-stage infection')
        ylabel('PMR')      
        set(gca,'TickDir','out','box','off','fontname','Arial')
        text(0.25,2.7,'(B) Inference with 12 m sampling window ','fontsize',24)
        l = legend([px(1) px(2) px(3) px(5)],{'From reg. early','From reg. late','Observed PMR','True PMR'});
        set(l,'fontsize',16)
        
        y = log10(sum(totNonSeq,2));
        %
        % Figure 3A
        figure(1)
        
        for cc = 1:4
        a = area([cc+3/24 cc-3/24],[3 3]);
        a.FaceColor = grayColor;
        a.FaceAlpha = 0.5;
        a.EdgeAlpha = 0;   
        a.ShowBaseLine='off';
        hold on
        a = area([cc+3/24 cc-3/24],[-3 -3]);
        a.FaceColor = grayColor;
        a.FaceAlpha = 0.5;
        a.EdgeAlpha = 0;   
        a.ShowBaseLine='off';
        end
        
        plot(t/24,log10(sum(totNonSeq,2)),'k','linewidth',3)
        
        % Log-linear regression
        peR = plot([0; t1_early6],b_early6*[0; t1_early6]+a_early6,'-','linewidth',2,'color','b'); 
        plR = plot([0; t1_late6],b_late6*[0; t1_late6]+a_late6,'--','linewidth',2,'color','b');         
        
        
        for j = 1
            pe = plot([t(td_early6(j)) t(td_early6(j+2))]/24,[y(td_early6(j)) y(td_early6(j+2))],'r-','markersize',8,'linewidth',1.5,'markerfacecolor','w'); 
            pl = plot([t(td_late6(j)) t(td_late6(j+2))]/24,[y(td_late6(j)) y(td_late6(j+2))],'r--','markersize',30,'linewidth',1.5); 
        end
        for j = 1:2
            plot([t(td_early6(j)) t(td_early6(j+2))]/24,[y(td_early6(j)) y(td_early6(j+2))],'ko','markersize',8,'linewidth',1.5,'markerfacecolor','w'); 
            plot([t(td_late6(j)) t(td_late6(j+2))]/24,[y(td_late6(j)) y(td_late6(j+2))],'k.','markersize',30,'linewidth',1.5); 
        end 
        
        hold off
        set(gca,'fontsize',20)
        xlim([0.5 3.5])
        xlabel('Day of blood-stage infection')
        ylabel('Circulating iRBC abundance')
        ylim([-3 3])
        set(gca,'yticklabel',{'10^0','10^1','10^2','10^3','10^4','10^5','10^6'})
        set(gca,'xtick',0:4,'box','off')
        l = legend([pl pe plR peR],{'Max. obs. late','Max. obs. early','From reg. late','From reg. early'});
        set(l,'location','northwest','fontsize',16)
        
    end
    
end

