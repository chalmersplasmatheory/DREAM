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

J = PetscBinaryRead('petsc_jacobian');
Jn = PetscBinaryRead('petsc_jacobian_num');

disp(['J  density: ',num2str(nnz(J)/numel(J) * 100),'%']);
disp(['Jn density: ',num2str(nnz(Jn)/numel(Jn) * 100),'%']);

%% Plot
MARKERSIZE = 10;
figure(1), spy(J, MARKERSIZE), title('J'), axis([1,6,1,6]);
figure(2), spy(Jn, MARKERSIZE), title('Jn'), axis([1,6,1,6]);

%% Check differences
skipRows = 0;
tol = 1e-2;
MARKERSIZE = 10;

DD = J-Jn;
[ii,jj,ss] = find(DD((skipRows+1):end,:));
ii = ii + skipRows;
for k=1:length(ii)
    Delta = abs(DD(ii(k), jj(k)) / J(ii(k),jj(k)));
    
    if Delta < tol
        DD(ii(k),jj(k)) = 0;
    end
end

spy(DD, MARKERSIZE);
%axis([1,6,1,6]);
