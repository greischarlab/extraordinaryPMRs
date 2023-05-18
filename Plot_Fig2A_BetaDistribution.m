%% Code for Figure 2A - Initial synchrony distributions
clear; close all; clc

%% Determine the bursting duration (in hrs)

sp_vec = fliplr([1 100]);
for sp1 = 1:length(sp_vec)
    sp = sp_vec(sp1);
    % Beta distribution from 0.5% to 99.5%
    y = betainv([0.005 .995],sp,sp)*48;
    Burst99(sp1) = y(2) - y(1);
end


%% Plot the beta distribution 

cyclelength = 48;
split = 6;
Cycle = round(cyclelength*split);
x = linspace(0,1,Cycle+1);
IC = 1;
ICshapevec = fliplr([1 100 ]);
ICshift = 0.75;
for j2 = 1:length(ICshapevec)
    ICshape = ICshapevec(j2);
    
    ICdist = betapdf(x(2:end),ICshape,ICshape)/Cycle;
    n0 = zeros(Cycle,1);
    id = round(Cycle*ICshift);
    n0(id+1:Cycle) = ICdist(1:Cycle-id)*IC;
    n0(1:id) = ICdist(Cycle-id+1:Cycle)*IC;
    
    f1 = figure(1);
    grayColor = [.7 .7 .7];
    if ICshape==500
        plot((1:Cycle)/split,n0,'--','color',grayColor,'linewidth',3)
    elseif ICshape==100
        plot((1:Cycle)/split,n0,'k-','linewidth',3)
    elseif ICshape==1
        plot((1:Cycle)/split,n0,':','color',0.5*grayColor,'linewidth',3)
    else
        plot((1:Cycle)/split,n0,'linewidth',2)
    end
    
    hold on
    set(gca,'fontsize',24)
    xlim([0 cyclelength])
    f{j2} = sprintf('%d hrs',round(Burst99(j2)));
end
hold off
l = legend(f,'fontsize',20,'location','northeast');
set(l,'box','off')
xlim([0 48])
ylim([0 .05])
set(gca,'xtick',0:12:48)
ylabel('Fraction of population')
xlabel('Age (hrs)')
set(gca,'TickDir','out','box','off','fontname','Arial')
text(30,0.006,'asynchronous','fontsize',20,'color',0.5*grayColor)
text(13,0.038,'synchronous','fontsize',20,'color','k')
text(30,0.05,'Bursting duration','fontsize',20)
text(1,0.053,'(A) Initial level of synchrony','fontsize',24)

