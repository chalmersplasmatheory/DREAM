% Calculate the conductivity in various cases with CODE
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
EEc      = 1e-2;
%EEc      = 50;
T0       = 1e3;

Ny  = 2000;
Nxi = 5;
yMax = 20;
Nt = 20;
nSaveSteps = 50;

timeUnit = 'ms';
tT = 0;
tn = 0;
tE = 0;
tZ = 0;

%Tarr = 1e3;
%Zarr = 1;
Tarr = [1e2, 1e3, 1e4, 45e3];
Zarr = [1, 2, 4, 8, 50];

sigma = zeros(length(Tarr), length(Zarr));
structs = cell(length(Tarr), length(Zarr));

for i=1:length(Tarr)
    for j=1:length(Zarr)
        vth = sqrt(2*e*Tarr(i)/me);
        
        lnLambda = 14.9-0.5*log(n/1e20)+log(Tarr(i)/1e3);
        Ec       = n*lnLambda*e^3 / (4*pi*eps0^2*me*c^2);
        nu0      = n*e^4*lnLambda / (4*pi*eps0^2*me^2*c^3);
        
        E = EEc * Ec;
        tMax = 0.1 * (Tarr(i)/T0)^(3/2);
        
        o = CODESettings();
        o.SetPhysicalParameters(Tarr(i),n,Zarr(j),E,B,tT,tn,tZ,tE);
        o.timeUnitOfParams = timeUnit;

        tMaxN = tMax / nu0 * 1000;
        dt = (tMax/(Nt-1)) / nu0 * 1000;
        o.SetResolutionParameters(Nxi,Ny,yMax,dt,tMaxN,timeUnit);

        o.sourceMode = 0;
        o.stepSkip = 10;
        o.nStepsToReturn = nSaveSteps;
        o.collisionOperator = 0;
        o.showDistEvolutionPlot = 2;
        
        s = CODE_timeDependent(o);
        
        structs{i,j} = s;
        sigma(i,j) = s.currentDensity(end) / (s.EOverEc * Ec);
    end
end

[Z,T] = meshgrid(Zarr, Tarr);

save('CODE-conductivities.mat', '-v7.3', 'sigma', 'T', 'Z');

