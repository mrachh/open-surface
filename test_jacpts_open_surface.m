run ~/git/fmm3dbie/matlab/startup.m
ifac = 3;
npars = [3;3;5];
S = geometries.disk([1;1], 0.3,  npars, 8, 11);

ipatch = npars(1)*npars(1) + npars(2)*7;
figure(1);
clf
plot(S, rand(S.npatches,1), 'FaceAlpha', 0.3);
hold on;
iind = S.ixyzs(ipatch):S.ixyzs(ipatch+1)-1;
xp = S.r(1,iind); xp = xp(:);
yp = S.r(2,iind); yp = yp(:);
zp = S.r(3,iind); zp = zp(:);
scatter3(S.r(1,iind), S.r(2,iind), S.r(3,iind));

f = @(u,v) adap_int(u, v, S.srccoefs{ipatch}, S.norders(ipatch));
% now build adaptive integrator
f1 = integral2(f, -1, 1, -1, 1, 'AbsTol', 1e-12, 'RelTol', 1e-9);

%% now test jacpts
[xj, wj] = jacpts(3, -0.5, 0);
[xl, wl] = legpts(3);

[ujl, vjl] = meshgrid(xl, xj);
[wjlx, wjly] = meshgrid(wl, wj);

w = wjlx(:).*wjly(:);


[~, ff] = adap_int(ujl, vjl, S.srccoefs{ipatch}, S.norders(ipatch));
fcomp = sum(ff(:).*w(:));
fprintf('error in computation=%d\n', norm(fcomp-f1))

function [f, ff] = adap_int(u, v, srccoefs, norder)
    uu = u(:).';
    vv = v(:).';
    uv = [uu; vv];
    srcvals = srccoefs*polytens.lege_pols(norder, uv);
    rr = 1- (srcvals(1,:).^2 + srcvals(2,:).^2);
    rr = rr(:);
    dj = cross(srcvals(4:6,:), srcvals(7:9,:));
    djnorm = vecnorm(dj, 2, 1);
    f = djnorm(:)./sqrt(rr);
    f = reshape(f, size(u));
    if nargout > 1
        ff = f.*sqrt(1-v);
    end
end