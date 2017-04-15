clear;
clc;

addpath('..')

% params = [L, A, VT, V0]
params = [4*pi, 0.5, 0, 0.2];

DT = 0.01;
NT = 100/DT;
NG = 128;
N = 1e4;

PIC.movie_run(DT, NT, NG, N, 'Two_Stream_Basic', params, inf, 1);

rmpath('..')