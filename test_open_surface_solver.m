run ~/git/fmm3dbie/matlab/startup.m
ifac = 3;
S = geometries.disk([1;1], 0.3,  [4;10;10], 8, 11);


zk = 2.1;

src = [0.3;0.3;1.1];
srcinfo = [];
srcinfo.r = src;
f = helm3d.kern(zk, srcinfo, S, 's');

load('fcoefs.mat');

r = sqrt(S.r(1,:).^2 + S.r(2,:).^2);
rr = 2*r -1;
rr = rr(:);
ncoef = length(fcoefs);
p = lege.pols(rr, ncoef-1);
f2 = p.'*fcoefs;
ff = f2(1:24);

errs = surf_fun_error(S, f2);
figure(1)
clf
plot(S, log10(errs)); colorbar();


%%

eps = 1e-5;
ifscal = 0;

Q = helm3d.opensurface.get_quadrature_correction(S, eps, ifscal, zk);
A = reshape(Q.wnear, [S.npts, S.npts]); A = A.';
clear Q

%%
[sig, flag, relres, iter] = gmres(A, f, [], 1e-9, 400);
figure(2)
clf
plot(S, real(sig));

errs = surf_fun_error(S, sig);
figure(3)
clf
plot(S, log10(abs(errs))); colorbar();


%%

ifscal = 1;
Qscal = helm3d.opensurface.get_quadrature_correction(S, eps, ifscal, zk);
Ascal = reshape(Qscal.wnear, [S.npts, S.npts]); Ascal = Ascal.';
[sig_scal, flag_scal, relres, iter_scal] = gmres(Ascal, f2, [], 1e-9, 400);

%%
opts = [];
opts.quadtype = 'ggq';
eps = 1e-6;
Qscal_ggq = helm3d.opensurface.get_quadrature_correction(S, eps, ifscal, zk, opts);
Ascal_ggq = reshape(Qscal_ggq.wnear, [S.npts, S.npts]); Ascal_ggq = Ascal_ggq.';


%%

[sig_scal_ggq, flag_scal_ggq, relres_scal_ggq, iter_scal_scal_ggq] = gmres(Ascal_ggq, f2, [], 1e-9, 400);

%%
figure(4)
clf
plot(S, log10(abs(sig_scal-1))); colorbar();

f3 = Ascal*ones(S.npts,1);
figure(5)
clf
plot(S, log10(abs(f2-f3))); colorbar();

%%
f4 = Ascal_ggq*ones(S.npts,1);
figure(6)
clf
plot(S, log10(abs(f2-f4))); colorbar();

figure(7)
clf
plot(S, abs(f2-f4)); colorbar();

figure(8)
clf
plot(S, abs(f2)); colorbar();

figure(9)
clf
plot(S, abs(f4)); colorbar();



%%

rho = 0.01;
A2 = A + rho*eye(S.npts);
sig2 = A2 \ f;
figure(6)
clf
plot(S, log10(abs(sig2)));

errs2 = surf_fun_error(S, sig2);
figure(7)
clf
plot(S, log10(abs(errs2))); colorbar();




