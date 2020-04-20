clc; clear;

% Path to PETSc
addpath('/opt/petsc/linux-c-opt/share/petsc/matlab/');
%addpath('/mnt/HDD/lib/petsc-3.11.1/share/petsc/matlab/');

A = PetscBinaryRead('petsc_matrix');

disp(['Matrix density: ',num2str(nnz(A)/numel(A) * 100),'%']);