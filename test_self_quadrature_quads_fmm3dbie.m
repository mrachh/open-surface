% This code tests the computation of self-quadrature at a point
% on a patch using ggq to adaptive quadrature

norder = 8;
S = geometries.sphere(1, 2, [], norder, 11);
[srcvals, ~, ~, ~, ~, ~] = extract_arrays(S);
npols = S.ixyzs(2) - S.ixyzs(1);
Suse = surfer(1, S.norders(1), srcvals(:,1:npols), S.iptype(1));

ipt = 5;
targinfo = [];
targinfo.r = Suse.r(:,ipt);
targinfo.du = Suse.du(:,ipt);
targinfo.dv = Suse.dv(:,ipt);
targinfo.patch_id = Suse.patch_id(ipt);
targinfo.uvs_targ = Suse.uvs_targ(:,ipt);


zk = 1.1;
rep_pars = [1;0];

eps = 1e-12;
Q = helm3d.dirichlet.get_quadrature_correction(Suse, eps, zk, rep_pars, ...
                                                targinfo);


%% Now split into 4 different patches
rnodes = polytens.lege_nodes(Suse.norders(1));
pols = polytens.lege_pols(Suse.norders(1), rnodes);
umat = polytens.lege_vals2coefs(norder, rnodes);

rnodes_tmp = (rnodes + [1;1])/2;

v0 = [-1;-1];
v1 = [targinfo.uvs_targ(1); -1];
v2 = [-1; targinfo.uvs_targ(2)];

rnodes_use1 = v0 + rnodes_tmp(1,:).*(v1-v0) + rnodes_tmp(2,:).*(v2-v0);

pols1 = polytens.lege_pols(Suse.norders(1), rnodes_use1);
srcvals_tmp = Suse.srccoefs{1}*pols1;
dr = srcvals_tmp(1:3,:);
du = srcvals_tmp(4:6,:)*(v1(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v1(2)-v0(2))/2;
dv = srcvals_tmp(4:6,:)*(v2(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v2(2)-v0(2))/2;
dn = cross(du, dv);
dj = vecnorm(dn, 2, 1);
dn = dn./dj;
srcvals_use = [dr; du; dv; dn];
S1 = surfer(1, Suse.norders(1), srcvals_use, Suse.iptype(1));

%%


v0 = [targinfo.uvs_targ(1); -1];
v1 = [targinfo.uvs_targ(1); targinfo.uvs_targ(2)];
v2 = [1; -1];


rnodes_use2 = v0 + rnodes_tmp(1,:).*(v1-v0) + rnodes_tmp(2,:).*(v2-v0);

pols2 = polytens.lege_pols(Suse.norders(1), rnodes_use2);
srcvals_tmp = Suse.srccoefs{1}*pols2;
dr = srcvals_tmp(1:3,:);
du = srcvals_tmp(4:6,:)*(v1(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v1(2)-v0(2))/2;
dv = srcvals_tmp(4:6,:)*(v2(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v2(2)-v0(2))/2;
dn = cross(du, dv);
dj = vecnorm(dn, 2, 1);
dn = dn./dj;
srcvals_use = [dr; du; dv; dn];

S2 = surfer(1, Suse.norders(1), srcvals_use, Suse.iptype(1));

%%
v0 = [-1; targinfo.uvs_targ(2)];
v1 = [targinfo.uvs_targ(1); targinfo.uvs_targ(2)];
v2 = [-1; 1];


rnodes_use3 = v0 + rnodes_tmp(1,:).*(v1-v0) + rnodes_tmp(2,:).*(v2-v0);

pols3 = polytens.lege_pols(Suse.norders(1), rnodes_use3);
srcvals_tmp = Suse.srccoefs{1}*pols3;
dr = srcvals_tmp(1:3,:);
du = srcvals_tmp(4:6,:)*(v1(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v1(2)-v0(2))/2;
dv = srcvals_tmp(4:6,:)*(v2(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v2(2)-v0(2))/2;
dn = cross(du, dv);
dj = vecnorm(dn, 2, 1);
dn = dn./dj;
srcvals_use = [dr; du; dv; dn];
S3 = surfer(1, Suse.norders(1), srcvals_use, Suse.iptype(1));

%%

v0 = [targinfo.uvs_targ(1); targinfo.uvs_targ(2)];
v1 = [targinfo.uvs_targ(1); 1];
v2 = [1; targinfo.uvs_targ(2)];


rnodes_use4 = v0 + rnodes_tmp(1,:).*(v1-v0) + rnodes_tmp(2,:).*(v2-v0);

pols4 = polytens.lege_pols(Suse.norders(1), rnodes_use4);
srcvals_tmp = Suse.srccoefs{1}*pols4;
dr = srcvals_tmp(1:3,:);
du = srcvals_tmp(4:6,:)*(v1(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v1(2)-v0(2))/2;
dv = srcvals_tmp(4:6,:)*(v2(1) - v0(1))/2 + srcvals_tmp(7:9,:)*(v2(2)-v0(2))/2;
dn = cross(du, dv);
dj = vecnorm(dn, 2, 1);
dn = dn./dj;
srcvals_use = [dr; du; dv; dn];
S4 = surfer(1, Suse.norders(1), srcvals_use, Suse.iptype(1));

%%

targinfo_off = [];
targinfo_off.r = Suse.r(:,ipt);
targinfo_off.du = Suse.du(:,ipt);
targinfo_of.dv = Suse.dv(:,ipt);
Q1 = helm3d.dirichlet.get_quadrature_correction(S1, eps, zk, rep_pars, ...
                                                targinfo_off);
Q2 = helm3d.dirichlet.get_quadrature_correction(S2, eps, zk, rep_pars, ...
                                                targinfo_off);
Q3 = helm3d.dirichlet.get_quadrature_correction(S3, eps, zk, rep_pars, ...
                                                targinfo_off);
Q4 = helm3d.dirichlet.get_quadrature_correction(S4, eps, zk, rep_pars, ...
                                                targinfo_off);


ufinal = umat.'*(pols1*Q1.wnear + pols2*Q2.wnear + pols3*Q3.wnear + pols4*Q4.wnear);


err1 = norm(ufinal - Q.wnear);
fprintf('error = %d\n',err1);


aa2 = log10(abs(pols*(Q.wnear - ufinal)));

[ii, jj] = meshgrid(0:8);
itot = ii(:) + jj(:);
aa2 = aa2(itot<=8);