%% Runs ODE system 
function [t, tot, totSeq, totNonSeq, p] = RunRBC(mu,museq,cyclelength,ICshift,ICshape,R)
% Runs ODE system found in 'RBC' (at end of file)
% Input:    mu: mortality of non-sequestered iRBCs
%           museq: mortality of sequestered iRBCs
%           cyclelength: cycle length in hours
%           ICshift: shift in initial parasite age distribution
%           ICshape: level of synchrony
% Output:   t: time points
%           tot: total iRBCs
%           tot: total sequestered iRBCs
%           tot: total non-sequestered iRBCs
%           p: parameter vector


split = 6; % intervals of 10 minutes
mu = mu*split;
museq = museq*split;

% True PMR

p.R = R;

% LR partly filled ring up to 50% of life cycle (WT 12-24hrs)
% ET small fully filled ring up to 62% of life cycle (WT 24-30hrs)
% LT larger ring and up to 2 nuclei up to 75% of life cycle (WT 30-36hrs)
% S 3+ schizonts up to 100% of life cycle (WT 36-48hrs)

Cycle = round(cyclelength*split);
R_end = round(1/4*Cycle);
LR_end = round(1/2*Cycle);
ET_end = round(5/8*Cycle);
T_end = round(3/4*Cycle);
p.Cycle = Cycle;

p.lambda = ones(Cycle,1);
p.lambdaS = ones(Cycle,1);
p.R = R;
p.mu = zeros(Cycle,1);
p.muS = zeros(Cycle,1);

% From fitting of Kriek at al (2003); See Archer et al (2018)
maxp = 467.6209;
p1 = 11.3869/maxp;
p2 = 467.6209;
p2 = p2/maxp;
p3 = 18.5802;
p4 = 0.2242;
x=[1/split:1/split:Cycle/split];
y = 1 - (p1 + (p2-p1) ./ (1 + 10.^((p3-x)*p4)));
p.q = [0 1-y(2:end)./y(1:end-1)];

% Rings 0 - 16
p.mu(1:R_end) = mu;
p.muS(1:R_end) = museq;
p.lambda(1:R_end) = split;
p.lambdaS(1:R_end) = split;
% Late rings 16 - 24
p.mu(R_end+1:LR_end) = mu;
p.muS(R_end+1:LR_end) = museq;
p.lambda(R_end+1:LR_end) = split;
p.lambdaS(R_end+1:LR_end) = split;
% Early Trop 24 - 32
p.mu(LR_end+1:ET_end) = mu;
p.muS(LR_end+1:ET_end) = museq;
p.lambda(LR_end+1:ET_end) = split;
p.lambdaS(LR_end+1:ET_end) = split;
% Trop 32 - 40
p.mu(ET_end+1:T_end) = mu;
p.muS(ET_end+1:T_end) = museq;
p.lambda(ET_end+1:T_end) = split;
p.lambdaS(ET_end+1:T_end) = split;
% Schizont 40-48
p.mu(T_end+1:Cycle) = mu;
p.muS(T_end+1:Cycle) = museq;
p.lambda(T_end+1:Cycle) = split;
p.lambdaS(T_end+1:Cycle) = split;

Maxtime = 200; % hours


% Initial distribution
n0 = zeros(Cycle*2,1);
x = linspace(0,1,Cycle+1);
IC = 1;
ICdist = betapdf(x(2:end),ICshape,ICshape)/Cycle;
id = round(Cycle*ICshift);
n0(id+1:Cycle) = ICdist(1:Cycle-id)*IC;
n0(1:id) = ICdist(Cycle-id+1:Cycle)*IC;

n0(Cycle+LR_end+1:2*Cycle) = n0(LR_end+1:Cycle);
n0(LR_end+1:Cycle) = 0;

options = odeset('NonNegative',1:Cycle*2);
[t,y] = ode45(@RBC,0:.1:Maxtime,n0,options,p);

NonSeq = y(:,1:Cycle); % Non-sequestered
Seq = y(:,Cycle+1:end); % Sequestered

tot = sum(y,2); % total parasitemia
totNonSeq = sum(NonSeq,2); % total non-sequestered parasitemia
totSeq = sum(Seq,2); % total sequestered parasitemia


end

%% ODE system 
function dy = RBC(t,y,p)
% ODE system
% Input:    t: time
%           y: state variables
%           p: parameter struture
% Output:   dy: time derivative of state variables

lam = p.lambda; % transition rate between non-sequestered states
lamS = p.lambdaS; % transition rate between sequestered states
R = p.R; % growth rate
mu = p.mu; % mortalty rate of non-sequestered states
muS = p.muS; % mortalty rate of sequestered states
q = p.q'; % proportin of transition of non-sequestered states
Cycle = p.Cycle; % total number of compartments before cycle is complete

n = y(1:Cycle);
nS = y(Cycle+1:2*Cycle);

dn(1,1) = R*(lam(Cycle)*n(Cycle)+lamS(Cycle)*nS(Cycle)) - (lam(1)+mu(1))*n(1);
ds(1,1) = 0;

dn(2:Cycle,1) = (1-q(1:end-1)).*lam(1:end-1).*n(1:end-1,1) - (lam(2:end)+mu(2:end)).*n(2:end,1);
ds(2:Cycle,1) = q(1:end-1).*lam(1:end-1).*n(1:end-1,1) + lamS(1:end-1).*nS(1:end-1,1) - (lamS(2:end)+muS(2:end)).*nS(2:end,1);

dy = [dn; ds];

end