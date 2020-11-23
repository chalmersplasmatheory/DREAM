% Calculate the runaway rate in various cases with CODE
clc; clear;

if ~exist('CODE_timeDependent', 'file')
    error('CODE has not been loaded. Please load CODE_timeDependent before running this script.');
end

% Constants
c = 299792458;
e = 1.60217662e-19;
eps0 = 8.85418782e-12;
me = 9.10938e-31;

% Plasma parameters
n        = 2e19; %m^{-3}
%B        = 2.5;
Z        = 1;
T        = 5e3;
E        = 0.04;

Barr = [1.5, 1.75, 2, 2.25, 2.5, 3];

Ny  = 2000;
Nxi = 80;
yMax = 400;
Nt = 20;
nSaveSteps = 50;

timeUnit = 'normalized';
tT = 0;
tn = 0;
tE = 0;
tZ = 0;

% Physical parameters
vth = sqrt(2*e*T/me);
lnLambda = 14.9-0.5*log(n/1e20)+log(T/1e3);
Ec       = n*lnLambda*e^3 / (4*pi*eps0^2*me*c^2);
nu0      = n*e^4*lnLambda / (4*pi*eps0^2*me^2*c^3);

bumpP  = zeros(size(Barr));

for i=1:length(Barr)
%for i=3
%for i=1
    % Run CODE
    o = CODESettings();
    o.SetPhysicalParameters(T,n,Z,E,Barr(i),tT,tn,tZ,tE);
    o.timeUnitOfParams = timeUnit;

    %tMax = (yMax*vth/c) / (E/Ec);
    tMax = 400e3;
    %tMaxN = tMax / nu0 * 1000;
    %dt = (tMax/(Nt-1)) / nu0 * 1000;
    dt = tMax/(Nt-1);
    o.SetResolutionParameters(Nxi,Ny,yMax,dt,tMax,timeUnit);

    o.sourceMode = 0;
    o.stepSkip = 10;
    o.nStepsToReturn = nSaveSteps;
    o.collisionOperator = 0;
    o.showDistEvolutionPlot = 2;
    %o.showDistEvolutionPlot = 0;
    o.yGridMode = 6;

    %o.showDistEvolutionPlot = 0;

    s = CODE_timeDependent(o);
    
    % Evaluate f(xi=1)
    distMatrix = reshape(s.f(:,end),s.Ny,s.Nxi);
    plusMin = ones(1,s.Nxi);
    plusMin(2:2:end) = -1;
    dist = sum(distMatrix*diag(plusMin),2);
    
    % Locate bump
    df = [0; diff(dist)];
    %[~,maxIdx] = max(df);
    
    maxIdx = length(df)-1;
    while df(maxIdx) > df(maxIdx+1) || df(maxIdx)*df(maxIdx-1) < 0 || df(maxIdx) < 0
        maxIdx = maxIdx-1;
    end
    
    bumpLocIdx = find(df(maxIdx:end) < 0,1);
    bumpP(i) = s.y(maxIdx+bumpLocIdx) * s.delta;
    
end

B = Barr;
save('CODE-bumps.mat', '-v7.3', 'bumpP', 'B');

%%

bumps = zeros(length(s.times),1);

for it=1:length(s.times)
    distMatrix = reshape(s.f(:,it),s.Ny,s.Nxi);
    plusMin = ones(1,s.Nxi);
    plusMin(2:2:end) = -1;
    dist = sum(distMatrix*diag(plusMin),2);
    
    % Locate bump
    df = [0; diff(dist)];
    [~,maxIdx] = max(df);
    bumpLocIdx = find(df(maxIdx:end) < 0,1);
    bumps(it) = s.y(maxIdx+bumpLocIdx) * s.delta; 
end

plot(s.times, bumps);
