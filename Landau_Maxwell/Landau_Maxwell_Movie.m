clear;
clc;

addpath('..')

% params = [L, A, VT, V0]
params = [4*pi, 0, 1, 0];

DT = 0.025;
NT = 100/DT;
NG = 512;
N = 1e7;

PIC.movie_run(DT, NT, NG, N, 'Landau_Maxwell', params, inf, 1);

rmpath('..')