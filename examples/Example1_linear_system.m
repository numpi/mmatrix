%  Solves a linear system using different factorizations
%
% Methods Compared:
%   - GTH: Pure MATLAB GTH algorithm (Right triplet)
%   - Standard LU: MATLAB's built-in LAPACK-based LU (dgetrf)
%   - Right Triplet MEX: dmmtf2 (unblocked), dmmtrf2 (recursive), 
%                        dmmtrf (blocked), dmmtrff (blocked recursive)
%   - Left Triplet MEX:  dmltf2 (unblocked), dmltrf2 (recursive), 
%                        dmltrf (blocked), dmltrff (blocked recursive)
%

clc;
fprintf('The condition number varies with σ = 10^{-k} for k = 1 to 12.\n\n');

% Problem parameters
n = 128;                   % Matrix size (n x n)
sigmas = 10.^(-(1:.5:12));   % Range of sigma values (affects conditioning)
nt = length(sigmas);
z = zeros(1,nt);
c = 1:nt; conds = z;
% Normwise errors
nrm_gth    = z; nrm_lu      = z;
nrm_dmmtf2 = z; nrm_dmmtrf2 = z; nrm_dmmtrf = z; nrm_dmmtrff = z;
nrm_dmltf2 = z; nrm_dmltrf2 = z; nrm_dmltrf = z; nrm_dmltrff = z;
% Componentwise errors
comp_gth    = z; comp_lu      = z;
comp_dmmtf2 = z; comp_dmmtrf2 = z; comp_dmmtrf = z; comp_dmmtrff = z;
comp_dmltf2 = z; comp_dmltrf2 = z; comp_dmltrf = z; comp_dmltrff = z;

for k = 1:nt

    sigma = sigmas(k);

    % Generate the FEM Matrix and right-hand side vector
    [objM,out] = gettriplet('fem', size=n, left=0, sigma=sigma, gamma=1, h=1/(n+1));
    b = out.b;
    A = objM.M;
    u = objM.u;
    v = objM.v;

    % Compute true solution using High-Precision VPA
    x_exact = double(vpa(objM) \ b);

    conds(k) = cond(A);

    % GTH with right triplet
    [L,U] = gthlu(A,u,v); x = U\(L\b);
    nrm_gth(k)  = norm(x-x_exact)/norm(x_exact);
    comp_gth(k) = max(abs(x-x_exact)./x_exact);

    % LU
    [L,U] = lu(A); x = U\(L\b);
    nrm_lu(k)  = norm(x-x_exact)/norm(x_exact);
    comp_lu(k) = max(abs(x-x_exact)./x_exact);

    % RIGHT TRIPLET MEX----------------------------------------------------
    x = mmex.dmmtf2_solve(A,u,v,b);
    nrm_dmmtf2(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmmtf2(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmmtrf2_solve(A,u,v,b);
    nrm_dmmtrf2(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmmtrf2(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmmtrf_solve(A,u,v,b); 
    nrm_dmmtrf(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmmtrf(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmmtrff_solve(A,u,v,b);% this equals to [L,U]=lu(objM); x=U\(L\b);
    nrm_dmmtrff(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmmtrff(k) = max(abs(x-x_exact)./x_exact);

    %LEFT TRIPLET MEX------------------------------------------------------
    x = mmex.dmltf2_solve(A,u,v,b); 
    nrm_dmltf2(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmltf2(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmltrf2_solve(A,u,v,b); 
    nrm_dmltrf2(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmltrf2(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmltrf_solve(A,u,v,b); 
    nrm_dmltrf(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmltrf(k) = max(abs(x-x_exact)./x_exact);

    x = mmex.dmltrff_solve(A,u,v,b); % this equals to [L,U]=lu(objM'); x=U\(L\b); (note that objM.M=objM.M')
    nrm_dmltrff(k)  = norm(x-x_exact)/norm(x_exact);
    comp_dmltrff(k) = max(abs(x-x_exact)./x_exact);

    fprintf('Sigma index %d/%d: cond(A) = %.2e\n', k, nt, conds(k));
end
figure
loglog(c,conds*eps,'c-','LineWidth', 2.5,'DisplayName','cond(A)\epsilon')
hold on

mk = 1:1:length(c);
loglog(c, nrm_lu,      '-',  'Color',[0.8 0 0], 'LineWidth',1,'MarkerIndices',mk,'DisplayName','LU norm')
loglog(c, comp_lu,     '--', 'Color',[0.8 0 0], 'LineWidth',1.2,'Marker','o','MarkerIndices',mk,'DisplayName','LU comp')
loglog(c, nrm_gth,     '-',  'Color',[1 0.4 0], 'LineWidth',2,'Marker','s','MarkerIndices',mk,'DisplayName','GTH norm')
loglog(c, comp_gth,    '--', 'Color',[1 0.4 0], 'LineWidth',1.2,'Marker','s','MarkerIndices',mk,'DisplayName','GTH comp')
loglog(c, nrm_dmmtf2,  '-',  'Color',[0 0.3 1], 'LineWidth',2,'Marker','o','MarkerIndices',mk,'DisplayName','dmmtf2 norm')
loglog(c, comp_dmmtf2, '--', 'Color',[0 0.3 1], 'LineWidth',1.2,'Marker','o','MarkerIndices',mk,'DisplayName','dmmtf2 comp')
loglog(c, nrm_dmmtrf2, '-',  'Color',[0 0.5 0.8], 'LineWidth',2,'Marker','s','MarkerIndices',mk,'DisplayName','dmmtrf2 norm')
loglog(c, comp_dmmtrf2,'--', 'Color',[0 0.5 0.8], 'LineWidth',1.2,'Marker','s','MarkerIndices',mk,'DisplayName','dmmtrf2 comp')
loglog(c, nrm_dmmtrf,  '-',  'Color',[0 0.7 0.6], 'LineWidth',2,'Marker','^','MarkerIndices',mk,'DisplayName','dmmtrf norm')
loglog(c, comp_dmmtrf, '--', 'Color',[0 0.7 0.6], 'LineWidth',1.2,'Marker','^','MarkerIndices',mk,'DisplayName','dmmtrf comp')
loglog(c, nrm_dmmtrff, '-',  'Color',[0 0.9 0.4], 'LineWidth',2,'Marker','d','MarkerIndices',mk,'DisplayName','dmmtrff norm')
loglog(c, comp_dmmtrff,'--', 'Color',[0 0.9 0.4], 'LineWidth',1.2,'Marker','d','MarkerIndices',mk,'DisplayName','dmmtrff comp')
loglog(c, nrm_dmltf2,  '-',  'Color',[0.6 0 0.8], 'LineWidth',2,'Marker','o','MarkerIndices',mk,'DisplayName','dmltf2 norm')
loglog(c, comp_dmltf2, '--', 'Color',[0.6 0 0.8], 'LineWidth',1.2,'Marker','o','MarkerIndices',mk,'DisplayName','dmltf2 comp')
loglog(c, nrm_dmltrf2, '-',  'Color',[0.8 0 0.8], 'LineWidth',2,'Marker','s','MarkerIndices',mk,'DisplayName','dmltrf2 norm')
loglog(c, comp_dmltrf2,'--', 'Color',[0.8 0 0.8], 'LineWidth',1.2,'Marker','s','MarkerIndices',mk,'DisplayName','dmltrf2 comp')
loglog(c, nrm_dmltrf,  '-',  'Color',[1 0 0.8], 'LineWidth',2,'Marker','^','MarkerIndices',mk,'DisplayName','dmltrf norm')
loglog(c, comp_dmltrf, '--', 'Color',[1 0 0.8], 'LineWidth',1.2,'Marker','^','MarkerIndices',mk,'DisplayName','dmltrf comp')
loglog(c, nrm_dmltrff, '-',  'Color',[1 0.2 1], 'LineWidth',2,'Marker','d','MarkerIndices',mk,'DisplayName','dmltrff norm')
loglog(c, comp_dmltrff,'--', 'Color',[1 0.2 1], 'LineWidth',1.2,'Marker','d','MarkerIndices',mk,'DisplayName','dmltrff comp')
legend('Location','bestoutside')
xlabel('Sigma index')
ylabel('Relative error')
title('Normwise vs componentwise errors for the solution of Mx=b')
hold off

