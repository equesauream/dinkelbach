%% run this script to test the dinkelbach algorithm

rng(42)
n = 10;   %  size Sn
rc = 0.5;   % reciprocal cond numb 
kind = 1;   % for pos def generation type
density = .6;

C = 10 * sprandsym(n,density,rc,kind);
lambda0 = 1;
tol = 1e-6;

Xopt = dinkelbach(C, lambda0, tol);
