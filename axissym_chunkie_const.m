zk = 2.1;
Sk = kernel('axissymh', 's', zk);

nleg = 20;
x = lege.exps(nleg);
x = (x+1)./2;
fvals = zeros(nleg,1);
for i = 1:nleg
   f = @(t) get_axissym_kernel(t, x(i), Sk);
   f1 = integral(f, 0, x(i), 'AbsTol', 1e-11, 'RelTol', 1e-8);
   f2 = integral(f, x(i), 1, 'AbsTol', 1e-11, 'RelTol', 1e-8);
   fvals(i) = f1 + f2;
end

%%
[~,w,u,v] = lege.exps(nleg);
fcoefs = u*fvals;
save('fcoefs.mat', 'fcoefs');
fileId = fopen('fcoefs.bin', 'w');
fcoefsr = zeros(2,20);
fcoefsr(1,:) = real(fcoefs.');
fcoefsr(2,:) = imag(fcoefs.');
fwrite(fileId, fcoefsr, 'double');
fclose(fileId);

A = load('~/git/fmm3dbie/examples/helmholtz/fort.34');
[npts, ~] = size(A);
fex = zeros(npts,1);
fcomp = zeros(npts,1);
for i=1:npts
   f = @(t) get_axissym_kernel(t, A(i,3), Sk);
   f1 = integral(f, 0, A(i,3), 'AbsTol', 1e-12, 'RelTol', 1e-8);
   f2 = integral(f, A(i,3), 1, 'AbsTol', 1e-12, 'RelTol', 1e-8);
   fex(i) = f1 + f2;
   fcomp(i) = A(i,1) + 1j*A(i,2);
end

function f = get_axissym_kernel(t, rho, Sk)
    srcinfo = [];
    t = t(:).';
    srcinfo.r = [t;zeros(1,length(t))];
    targinfo = [];
    targinfo.r = [rho; 0];
    f = Sk.eval(srcinfo, targinfo);
    f = f./sqrt(1-t.^2); 
end