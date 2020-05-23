clc; clear;

% Path to PETSc
if exist('/opt/petsc/share/petsc/matlab/', 'dir')
    % Mathias' home computer
    addpath('/opt/petsc/share/petsc/matlab/');
elseif exist('~/petsc/share/petsc/matlab/', 'dir')
    % Ola's laptop
    addpath('~/petsc/share/petsc/matlab/');
elseif exist('/mnt/HDD/lib/petsc-3.11.1/share/petsc/matlab/', 'dir')
    % Mathias' work computer
    addpath('/mnt/HDD/lib/petsc-3.11.1/share/petsc/matlab/');
end

A = PetscBinaryRead('petsc_matrix');

disp(['Matrix density: ',num2str(nnz(A)/numel(A) * 100),'%']);

%% E and PitchScatterTerm

A = PetscBinaryRead('petsc_matrix');
E = PetscBinaryRead('petsc_matrix_E');
P = PetscBinaryRead('petsc_matrix_P');

A = A(231:end,231:end);
E = E(231:end,231:end);
P = P(231:end,231:end);

dA = diag(A);
dE = diag(E);
dP = diag(P);

offs = 100;
figure(1), semilogy(abs(dE)), title('E-term'), xlim([1, 410]+offs);
figure(2), semilogy(abs(dP)), title('Pitch-term'), xlim([1, 100]+offs);
figure(3), semilogy(abs(dA)), title('Combined'), xlim([1, 100]+offs);

%%
semilogy(abs(dA-dP)), title('Combined'), xlim([401, 800]);

%% Solve
load grid

g = sqrt(1+xp1.^2);
Theta = 0.00195695;
fM = repmat(exp(-g/Theta), numel(xp2)*10, 1);
dt = 1e-3;
Vp = reshape(xVprime, [10000,1]);

T = speye(size(P)) / dt;

M = P+E;
%xf = (P+E+T) \ (fM/dt);
xf = (T+M) \ (fM/dt);

figure, hold on;
semilogy(xp1, fM(1:100), 'k:');
semilogy(xp1,xf(1:100));
set(gca, 'YScale', 'log');

dxf = sum((xf-fM).*Vp);

disp(['Conservation:          ',num2str(dxf/max(xf))]);
disp(['Matrix conservativity: ',num2str(sum(M*Vp)/max(max(M)))]);

%% 
load grid
A = PetscBinaryRead('petsc_matrix');
%spy(A), axis([0,200,0,200]);
%plot(diag(A(51:150,51:150)))

M = A(231:end, 231:end);
Vp = reshape(xVprime, [numel(xVprime),1]);

I = sum(M * Vp);
I0 = sum(Vp);

disp(['Conservativity: ',num2str(I/I0*100),'%']);

