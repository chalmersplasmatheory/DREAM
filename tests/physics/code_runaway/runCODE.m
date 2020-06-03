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
n        = 5e19; %m^{-3}
B        = 0;
T0       = 1e3;
Z        = 1;

Ny  = 2000;
Nxi = 80;
yMax = 250;
Nt = 50;
nSaveSteps = 50;

timeUnit = 'ms';
tT = 0;
tn = 0;
tE = 0;
tZ = 0;

%Tarr = 1e3;
%Zarr = 1;
nE   = 4;
Tarr = [1e2, 5e2, 1e3, 5e3, 1e4];

E = zeros(length(Tarr),nE);
fullRunawayRate = zeros(length(Tarr), nE, min(nSaveSteps, Nt));
runawayRate = zeros(length(Tarr), nE);
structs = cell(length(Tarr), nE);

for i=1:length(Tarr)
    vth = sqrt(2*e*Tarr(i)/me);

    lnLambda = 14.9-0.5*log(n/1e20)+log(Tarr(i)/1e3);
    Ec       = n*lnLambda*e^3 / (4*pi*eps0^2*me*c^2);
    ED       = Ec * (me*c^2) / (e*Tarr(i));
    nu0      = n*e^4*lnLambda / (4*pi*eps0^2*me^2*c^3);

    E0       = 2*Ec;
    EMax     = 0.04*ED;
    
    if EMax < E0
        error('EMax < E0: Perhaps the temperature is too high?');
    end
    
    Earr = linspace(E0, EMax, nE);
    E(i,:) = Earr;
    
    for j=1:nE
        pMax = yMax * (vth/c);
        %tMax = (pMax * (Tarr(i)/T0)^(3/2);
        tMax = 0.7 * pMax / (Earr(j)/Ec);
        
        o = CODESettings();
        o.SetPhysicalParameters(Tarr(i),n,Z,Earr(j),B,tT,tn,tZ,tE);
        o.timeUnitOfParams = timeUnit;

        tMaxN = tMax / nu0 * 1000;
        dt = (tMax/(Nt-1)) / nu0 * 1000;
        o.SetResolutionParameters(Nxi,Ny,yMax,dt,tMaxN,timeUnit);

        o.sourceMode = 0;
        o.stepSkip = 10;
        o.nStepsToReturn = nSaveSteps;
        o.collisionOperator = 0;
        %o.showDistEvolutionPlot = 2;
        o.showDistEvolutionPlot = 0;
        
        s = CODE_timeDependent(o);
        
        structs{i,j} = s;
        
        rr = s.growthRatePerSecond .* s.fracRE .* s.density;
        rr(isnan(rr)) = 0;
        
        % Check that runaway-rate reached steady state
        variation = abs(rr(end)/interp1(s.times, rr, s.times(end)*4/5) - 1);
        if variation > 1e-2
            warning('Runaway rate does not seem to have converged.');
        end
        
        fullRunawayRate(i,j,:) = rr;
        runawayRate(i,j) = rr(end);
    end
end

T = repmat(Tarr', [1,nE]);

save('CODE-rates.mat', '-v7.3', 'runawayRate', 'T', 'E');

