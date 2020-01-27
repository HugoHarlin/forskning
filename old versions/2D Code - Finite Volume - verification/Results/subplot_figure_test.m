
clear all;
close all;
clc;
az =0;
el = 90;
figure(1)
ylabel('Water depth');
xlabel('Distance from lake center');

A  = [ 1 2 3; 4 5 6; 7 8 9];
B  = 2.*[ 1 2 3; 4 5 6; 7 8 9];
C  = 3.*[ 1 2 3; 4 5 6; 7 8 9];
D  = 4.*[ 1 2 3; 4 5 6; 7 8 9];


subaxis(2, 2,1, 'sh', 0.000, 'sv', 0.000, ...
    'PL', 0.100,  'PT', 0.10, 'PB', 0.100, 'PR', 0.100, ...
    'MT', 0.100,'MB', 0.100, 'ML', 0.1, 'MR', 0.1);
imagesc(A)

subaxis(2, 2,2, 'sh', 0.000, 'sv', 0.000, ...
    'PL', 0.10,  'PT', 0.10, 'PB', 0.100,'PR', 0.100, ...
    'MT', 0.10,'MB', 0.10, 'ML', 0.1, 'MR', 0.10);
imagesc(B)

subaxis(2, 2,3, 'sh', 0.000, 'sv', 0.000, ...
    'PL', 0.100,  'PT', 0.10, 'PB', 0.100,'PR', 0.100, ...
    'MT', 0.100,'MB', 0.100, 'ML', 0.1, 'MR', 0.10);
imagesc(C)

subaxis(2, 2,4, 'sh', 0.000, 'sv', 0.000, ...
    'PL', 0.100,  'PT', 0.10, 'PB', 0.100, 'PR', 0.100, ...
    'MT', 0.100,'MB', 0.100, 'ML', 0.1, 'MR', 0.10);
imagesc(D)