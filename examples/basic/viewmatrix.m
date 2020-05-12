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