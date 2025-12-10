function sfwm_maps_signedY_scan() 
% Multi-geometry SFWM scan on (λp, signed-λ) with one map per geometry,
% plus interactive first-order JSA/JSI (user-chosen λp,λs) across all L.
% The code computes Schmidt number per L, keeps only the L with K closest
% to 1, and shows that case. It also reports the Eq.(8) predictor.

%% ===================== [0] USER PARAMETERS =====================
tsv_file        = 'datos_modos_L-3.tsv';
polarization    = 'TE';
TE_frac_min     = 0.99;

save_png        = true;
outdir          = fullfile(pwd,'maps_signedY_scan');

lambda_p_nm     = 720:0.1:890;
y_signed_nm     = -1570:0.1:1570;

L_list_mm       = 0.1:0.01:60;

poly_deg        = 30;
omega_grid_pts  = 4001;

use_brewermap   = false;

% Extra plots para la presentación
plot_fit_debug  = true;   % n_eff(λ) + β1(λ) + β2(λ) + K(L)

%% ===================== [1] UTILS & CONSTANTS ===================
c  = 299792458; um = 1e-6; nm = 1e-9;
to_um = @(m) m*1e6;  to_nm = @(m) m*1e9;
omega_from_lambda = @(lam_m) 2*pi*c./lam_m;
lambda_from_omega = @(w)     2*pi*c./w;
sinc_std = @(x) (sin(x)./x).*(x~=0) + 1*(x==0);

Gamma_eq8 = 0.193;  % from the paper’s Eq.(8)

if save_png && ~exist(outdir,'dir'), mkdir(outdir); end

%% ===================== [2] LOAD TSV & ENUMERATE GEOMS ==========
DATA = readmatrix(tsv_file,'FileType','text','Delimiter','\t');
if isempty(DATA), error('Could not read TSV: %s', tsv_file); end

Alt      = DATA(:,1);
Anc      = DATA(:,2);
lam_um   = DATA(:,3);
te_frac  = nan(size(lam_um)); if size(DATA,2) >= 7, te_frac = DATA(:,7); end

pairs_all = unique([Alt Anc],'rows');
geom_list = [];
for k = 1:size(pairs_all,1)
    m = Alt==pairs_all(k,1) & Anc==pairs_all(k,2);
    if TE_frac_min<=0 || mean(te_frac(m),'omitnan') >= TE_frac_min
        geom_list = [geom_list; pairs_all(k,:)]; %#ok<AGROW>
    end
end
if isempty(geom_list)
    warning('No geometries pass TE_frac_min=%.3f; running first geometry.', TE_frac_min);
    geom_list = pairs_all(1,:);
end

%% ===================== [3] LOOP GEOMETRIES ======================
for ig = 1:size(geom_list,1)

    Alt_um0 = geom_list(ig,1);
    Anc_um0 = geom_list(ig,2);
    mask_g  = Alt==Alt_um0 & Anc==Anc_um0;

    % ---- neff data for this geometry ----
    switch lower(polarization)
        case 'te', neff_all = DATA(:,4);
        case 'tm', assert(size(DATA,2)>=9,'TM neff expected in column 9.'); neff_all = DATA(:,9);
        otherwise, error('polarization must be TE or TM');
    end
    lam_um_g = lam_um(mask_g);
    neff_g   = neff_all(mask_g);

    good = isfinite(lam_um_g) & isfinite(neff_g);
    lam_um_g = lam_um_g(good);  neff_g = neff_g(good);
    if numel(lam_um_g) < 6
        fprintf('Geom Alt=%.3f Anc=%.3f: too few λ points, skipping.\n', Alt_um0, Anc_um0);
        continue;
    end

    [lam_um_u,~,ic] = unique(lam_um_g);
    neff_u = accumarray(ic, neff_g, [], @mean);
    [lam_um_u,Iu] = sort(lam_um_u); neff_u = neff_u(Iu);
    lam_min_um = min(lam_um_u); lam_max_um = max(lam_um_u);

    % ---- fit k(ω) with normalization ----
    lam_m     = lam_um_u*um;
    omega_smp = omega_from_lambda(lam_m);
    k_smp     = (neff_u .* omega_smp) / c;

    deg = poly_deg; fitted=false;
    while deg >= 4 && ~fitted
        try
            [p,~,mu] = polyfit(omega_smp, k_smp, deg);
            fitted = true;
        catch
            deg = deg - 2;
        end
    end
    if ~fitted
        fprintf('Geom Alt=%.3f Anc=%.3f: k(ω) fit failed, skipping.\n', Alt_um0, Anc_um0);
        continue;
    end
    p_der  = polyder(p);
    p_der2 = polyder(p_der);  % segunda derivada del polinomio

    % evaluators (no interp cropping)
    k_of_omega      = @(w)        polyval(p,     (w - mu(1))/mu(2));
    beta1_of_omega  = @(w) (1/mu(2))*polyval(p_der, (w - mu(1))/mu(2));       % s/m
    beta2_of_omega  = @(w) (1/mu(2)^2)*polyval(p_der2,(w - mu(1))/mu(2));     % s^2/m

    %% ---------- Optional debug plots: n_eff(λ), β1(λ), β2(λ) ----------
    if plot_fit_debug
        % Fine grid in ω over the data range
        omega_fine = linspace(min(omega_smp), max(omega_smp), 1000);
        
        % Polynomial k(ω) on fine grid
        k_fit = k_of_omega(omega_fine);         % [1/m]
        
        % Convert back to λ and n_eff
        lam_fit_m  = lambda_from_omega(omega_fine);
        lam_fit_um = lam_fit_m / um;
        neff_fit   = (k_fit .* c) ./ omega_fine;   % since k = n_eff * ω / c
        
        % β1(ω) and β2(ω) from polynomial derivatives
        beta1_fit = beta1_of_omega(omega_fine);      % [s/m]
        beta2_fit = beta2_of_omega(omega_fine);      % [s^2/m]
        % Units conversion:
        % β1: s/m -> ps/mm
        beta1_ps_per_mm = beta1_fit * 1e12 * 1e-3;
        % β2: s^2/m -> fs^2/mm
        beta2_fs2_per_mm = beta2_fit * 1e30 * 1e-3;
        
        % Also compute β1 on the original sample points (for comparison)
        beta1_data = beta1_of_omega(omega_smp);
        beta1_data_ps_per_mm = beta1_data * 1e12 * 1e-3;
        
        % ===== Figure 1: n_eff(λ) with polynomial fit + residuals =====
        figure('Color','w','Position',[100 100 900 380]);
        subplot(1,2,1);
        plot(lam_um_u*1e3, neff_u, 'o', 'MarkerSize', 5, 'DisplayName','data'); hold on;
        plot(lam_fit_um*1e3, neff_fit, '-', 'LineWidth', 1.8, ...
             'DisplayName',sprintf('poly fit (deg %d)',deg));
        grid on;
        xlabel('\lambda (nm)');
        ylabel('n_{\mathrm{eff}}');
        title(sprintf('n_{eff} vs. \\lambda  —  Alt=%.3f \\mum, Anc=%.3f \\mum', Alt_um0, Anc_um0));
        legend('Location','best');
        
        % Residuales (neff_data - neff_fit at sample points)
        k_fit_on_data      = k_of_omega(omega_smp);
        neff_fit_on_data   = (k_fit_on_data .* c) ./ omega_smp;
        res_neff           = neff_u - neff_fit_on_data;
        
        subplot(1,2,2);
        stem(lam_um_u*1e3, res_neff, 'filled');
        grid on;
        xlabel('\lambda (nm)');
        ylabel('\Delta n_{\mathrm{eff}} (data - fit)');
        title('Fit residuals');
        
        % ===== Figure 2: β1(λ) (group delay per length) =====
        figure('Color','w','Position',[120 120 720 400]);
        hold on;
        plot(lam_fit_um*1e3, beta1_ps_per_mm, '-', 'LineWidth', 1.8, ...
             'DisplayName','\beta_1 from polynomial');
        plot(lam_um_u*1e3, beta1_data_ps_per_mm, 'o', 'MarkerSize', 5, ...
             'DisplayName','\beta_1 on data \omega');
        grid on;
        xlabel('\lambda (nm)');
        ylabel('\beta_1 (ps/mm)');
        title(sprintf('\\beta_1 vs. \\lambda  —  Alt=%.3f \\mum, Anc=%.3f \\mum', Alt_um0, Anc_um0));
        legend('Location','best');
        
        % ===== Figure 3: β2(λ) (GVD) =====
        figure('Color','w','Position',[140 140 720 400]);
        plot(lam_fit_um*1e3, beta2_fs2_per_mm, 'LineWidth', 1.8);
        grid on;
        xlabel('\lambda (nm)');
        ylabel('\beta_2 (fs^2/mm)');
        title(sprintf('\\beta_2 vs. \\lambda  —  Alt=%.3f \\mum, Anc=%.3f \\mum', Alt_um0, Anc_um0));
    end

    % ---- phasematching map once (independent of L) ----
    [Xp_nm, Ysigned_nm] = meshgrid(lambda_p_nm, y_signed_nm);
    lambda_p_m = Xp_nm*nm;
    omega_p    = omega_from_lambda(lambda_p_m);

    omega_s = nan(size(Xp_nm)); omega_i = nan(size(Xp_nm));
    posY = Ysigned_nm > 0;   negY = Ysigned_nm < 0;

    lambda_s_pos_m = (Ysigned_nm(posY))*nm;
    omega_s(posY)  = omega_from_lambda(lambda_s_pos_m);
    omega_i(posY)  = 2*omega_p(posY) - omega_s(posY);

    lambda_i_neg_m = abs(Ysigned_nm(negY))*nm;
    omega_i(negY)  = omega_from_lambda(lambda_i_neg_m);
    omega_s(negY)  = 2*omega_p(negY) - omega_i(negY);

    lambda_s_m = lambda_from_omega(omega_s);
    lambda_i_m = lambda_from_omega(omega_i);

    bad = ~isfinite(omega_s) | ~isfinite(omega_i) | (omega_s<=0) | (omega_i<=0);
    lambda_s_m(bad) = NaN; lambda_i_m(bad) = NaN;

    in_support = (to_um(lambda_s_m) >= lam_min_um) & (to_um(lambda_s_m) <= lam_max_um) & ...
                 (to_um(lambda_i_m) >= lam_min_um) & (to_um(lambda_i_m) <= lam_max_um) & ...
                 ~isnan(lambda_s_m) & ~isnan(lambda_i_m);

    DeltaK = nan(size(Xp_nm)); Tau_s  = nan(size(Xp_nm));
    Tau_i  = nan(size(Xp_nm)); Tau_si = nan(size(Xp_nm));

    idx = find(in_support);
    kp  = k_of_omega(omega_p(idx));
    ks  = k_of_omega(omega_s(idx));
    ki  = k_of_omega(omega_i(idx));
    DeltaK(idx) = ks + ki - 2*kp;

    b1p = beta1_of_omega(omega_p(idx));
    b1s = beta1_of_omega(omega_s(idx));
    b1i = beta1_of_omega(omega_i(idx));
    Tau_s(idx)  = (b1s - b1p);
    Tau_i(idx)  = (b1i - b1p);
    Tau_si(idx) = (b1s + b1i - 2*b1p);

    absDeltaK = abs(DeltaK); absDeltaK(~in_support) = NaN;

    vals = absDeltaK(in_support); vals = vals(isfinite(vals));
    if isempty(vals)
        clim=[0 1];
    else
        pr = prctile(vals,[5 95]);
        if any(~isfinite(pr)) || pr(1)>=pr(2)
            pr=[min(vals) max(vals)];
        end
        clim = pr;
    end

    figure('Name',sprintf('Map Alt=%.3f Anc=%.3f',Alt_um0,Anc_um0),...
           'Color','w','Position',[80 80 1120 620]);
    imagesc(lambda_p_nm, y_signed_nm, absDeltaK); set(gca,'YDir','normal'); hold on;
    cb = colorbar; ylabel(cb,'|Δk| (1/m)'); caxis(clim);
    xlabel('\lambda_p (nm)'); ylabel('signed \lambda (nm)');
    title(sprintf('Alt = %.3f \\mum, Anc = %.3f \\mum   |   map shown once', Alt_um0, Anc_um0));
    if use_brewermap
        try, colormap(flipud(brewermap(256,'RdBu'))); catch, end
    end
    contour(lambda_p_nm, y_signed_nm, DeltaK, [0 0], 'k', 'LineWidth', 1.6);
    contour(lambda_p_nm, y_signed_nm, Tau_s,  [0 0], 'g', 'LineWidth', 1.1);
    contour(lambda_p_nm, y_signed_nm, Tau_i,  [0 0], 'm', 'LineWidth', 1.1);
    contour(lambda_p_nm, y_signed_nm, Tau_si, [0 0], 'c', 'LineWidth', 1.3);
    legend({'|Δk|','Δk=0','τ_s=0','τ_i=0','τ_s+τ_i=0'}, ...
           'Location','southoutside','Orientation','horizontal');

    if save_png
        fn = sprintf('map_ONCE_Alt%.3f_Anc%.3f.png', Alt_um0, Anc_um0);
        exportgraphics(gcf, fullfile(outdir, fn), 'Resolution', 300);
    end

    %% --------- INTERACTIVE: pick λp, λs and build JSA ----------
    do_jsa_user = input('\nCompute JSA/JSI+Schmidt for THIS geometry? [y/N]: ','s');
    do_jsa_user = ~isempty(do_jsa_user) && lower(do_jsa_user(1))=='y';
    if ~do_jsa_user, continue; end

    lp_user = input('  Enter pump wavelength λp (nm)  : ');
    ls_user = input('  Enter signal wavelength λs (nm): ');
    if isempty(lp_user) || isempty(ls_user) || ~isfinite(lp_user) || ~isfinite(ls_user)
        warning('  Skipping: invalid λp/λs.');
        continue;
    end

    op = 2*pi*c / (lp_user*nm);
    os = 2*pi*c / (ls_user*nm);
    oi = 2*op - os;
    li_user = (2*pi*c/oi)/nm;

    pulse_tau_fs = input('  Pump FWHM (fs) [100]: ');       if isempty(pulse_tau_fs), pulse_tau_fs = 100; end
    N_jsa        = input('  Grid points per axis N [121]: '); if isempty(N_jsa), N_jsa = 121; end
    sigma_window = input('  Half-window in σ_p units [3.0]: '); if isempty(sigma_window), sigma_window = 3.0; end
    save_jsa_png = input('  Save figures? [y/N]: ','s');
    save_jsa_png = ~isempty(save_jsa_png) && lower(save_jsa_png(1))=='y';
    if save_jsa_png
        jsa_outdir = fullfile(outdir,'jsa_bestL');
        if ~exist(jsa_outdir,'dir'), mkdir(jsa_outdir); end
    end

    % Pump spectral sigma (transform-limited Gaussian)
    tau = pulse_tau_fs*1e-15;
    delta_nu_FWHM = 0.441 / tau;
    sigma_nu    = delta_nu_FWHM / (2*sqrt(2*log(2)));
    sigma_omega = 2*pi * sigma_nu;

    % Grid around (os,oi)
    w_span = sigma_window * sigma_omega;
    ws_vec = linspace(os - w_span, os + w_span, N_jsa);
    wi_vec = linspace(oi - w_span, oi + w_span, N_jsa);
    [WS, WI] = meshgrid(ws_vec, wi_vec);
    SUM   = WS + WI;
    nu_s  = WS - os;  nu_i = WI - oi;

    % α(ωs+ωi)
    Alpha = exp( - ((SUM - 2*op).^2) / (2*sigma_omega^2) );

    % group-delay slopes (per length)
    Ts0 = beta1_of_omega(os) - beta1_of_omega(op);   % s/m
    Ti0 = beta1_of_omega(oi) - beta1_of_omega(op);   % s/m
    theta_deg = -atan2d(Ts0, Ti0);  % purely for reporting

    % Eq.(8) predictor: 2 Γ σ^2 |τs τi| = 1, with τ• = L * (β1• - β1p)
    % => L_opt = 1 / ( σ * sqrt( 2 Γ |Ts0*Ti0| ) )
    Lopt_m  = 1 / ( sigma_omega * sqrt( 2*Gamma_eq8 * max(eps,abs(Ts0*Ti0)) ) );
    Lopt_mm = 1e3 * Lopt_m;

    % ridge line Δk_lin=0
    slope = -Ts0/Ti0;
    ws_line = linspace(min(ws_vec), max(ws_vec), 400);
    wi_line = oi + slope*(ws_line - os);
    ws_line_THz = ws_line/1e12; wi_line_THz = wi_line/1e12;

    % Sweep all L, compute K, keep best (closest to 1)
    dws = mean(diff(ws_vec)); dwi = mean(diff(wi_vec));
    best = struct('Lmm',NaN,'K',Inf,'pur',NaN,'Phi',[],'Alpha',[],'JSA',[],'JSI',[], ...
                  'U',[],'V',[],'lambda',[],'ws_THz',ws_vec/1e12,'wi_THz',wi_vec/1e12);

    K_list   = nan(size(L_list_mm));
    pur_list = nan(size(L_list_mm));
    R_list   = nan(size(L_list_mm));

    fprintf('\n--- L sweep (Schmidt) ---\n');
    for iL = 1:numel(L_list_mm)
        Lmm = L_list_mm(iL);  Lm = Lmm*1e-3;

        DeltaK_lin = Ts0.*nu_s + Ti0.*nu_i;      % 1/m
        Arg = (Lm/2) .* DeltaK_lin;
        Phi = sinc_std(Arg) .* exp(1i*Arg);

        JSA = Alpha .* Phi;
        JSA(~isfinite(JSA)) = 0;

        % Weighted SVD: M_ij = JSA_ij * sqrt(Δωs Δωi)
        M = JSA * sqrt(dws) * sqrt(dwi);  % equivalent to .*sqrt(dws*dwi)
        svals = svd(M,'econ');
        lambdas = (svals.^2);
        if ~any(isfinite(lambdas))
            K = NaN; purity = NaN;
        else
            lambdas = lambdas / sum(lambdas,'omitnan');
            K = 1/sum(lambdas.^2,'omitnan');
            purity = 1/K;
        end

        % Eq.(8) residual for this L (uses τs=L Ts0, τi=L Ti0)
        R_eq8 = abs( 2*Gamma_eq8 * sigma_omega^2 * (Lm^2) * (Ts0*Ti0) ) - 1;

        fprintf('  L=%.3f mm  ->  K=%.3f (pur=%.3f),  Eq(8) residual R=%.3g\n', ...
                Lmm, K, purity, R_eq8);

        K_list(iL)   = K;
        pur_list(iL) = purity;
        R_list(iL)   = R_eq8;

        % keep best |K-1|
        if abs(K-1) < abs(best.K-1)
            best.Lmm   = Lmm;
            best.K     = K;
            best.pur   = purity;
            best.Phi   = Phi;
            best.Alpha = Alpha;
            best.JSA   = JSA;
            best.JSI   = abs(JSA).^2;
            % store a few modes
            [U,S,V] = svd(M,'econ');
            lam = (diag(S).^2); lam = lam/sum(lam,'omitnan');
            best.lambda = lam;
            best.U = U; best.V = V;
        end
    end

    % --------- PLOT K(L) vs L with Eq.(8) predictor ----------
    if plot_fit_debug
        figure('Color','w','Position',[90 90 780 420]);
        hold on;
        hK    = plot(L_list_mm, K_list, '-','LineWidth',1.6,'DisplayName','K(L)');
        hBest = plot(best.Lmm, best.K,'ro','MarkerFaceColor','r',...
                     'DisplayName','Best (min |K-1|)');
        hK1   = yline(1,'k--','DisplayName','K=1');
        hLopt = xline(Lopt_mm,'r--','LineWidth',1.4,'DisplayName','L_{opt} Eq.(8)');
        grid on;
        xlabel('L (mm)');
        ylabel('Schmidt number K');
        title(sprintf('K vs L — Alt=%.3f \\mum, Anc=%.3f \\mum',Alt_um0,Anc_um0));
        legend([hK hBest hK1 hLopt],'Location','best');
    end

    % --------- SHOW ONLY THE BEST L ----------
    Lmm = best.Lmm; Lm = Lmm*1e-3;
    ws_THz = best.ws_THz; wi_THz = best.wi_THz;

    % 1) |Φ| panel (helps judge the diagonal width/angle)
    figure('Color','w','Position',[70 70 1180 620]);
    imagesc(ws_THz, wi_THz, abs(best.Phi)); set(gca,'YDir','normal'); colorbar;
    xlabel('\omega_s (THz)'); ylabel('\omega_i (THz)'); title('| \Phi | (first-order)');
    hold on;
    xl = xlim; yl = ylim;
    msk = ws_line_THz>=xl(1) & ws_line_THz<=xl(2) & ...
          wi_line_THz>=yl(1) & wi_line_THz<=yl(2);
    plot(ws_line_THz(msk), wi_line_THz(msk), 'w--', 'LineWidth', 1.3);
    sgtitle(sprintf(['Alt=%.3f, Anc=%.3f  |  L=%.3f mm  |  λ_p=%.1f nm, λ_s=%.1f nm, λ_i=%.1f nm\n' ...
                     'θ≈%.1f°;  Eq(8): L_{opt}≈%.3f mm,  K=%.3f (pur=%.3f)'], ...
        Alt_um0, Anc_um0, Lmm, lp_user, ls_user, li_user, ...
        theta_deg, Lopt_mm, best.K, best.pur), 'FontWeight','bold');

    if save_jsa_png
        tag = sprintf('BEST_Phi_Alt%.3f_Anc%.3f_L%.3fmm_lp%.1f_ls%.1f.png', ...
                      Alt_um0, Anc_um0, Lmm, lp_user, ls_user);
        exportgraphics(gcf, fullfile(jsa_outdir, tag), 'Resolution', 300);
    end

    % 2) JSA/JSI + Schmidt spectrum + first two modes
    figure('Color','w','Position',[80 80 1300 760]);
    t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    nexttile; imagesc(ws_THz, wi_THz, abs(best.JSA)); set(gca,'YDir','normal'); colorbar;
    xlabel('\omega_s (THz)'); ylabel('\omega_i (THz)'); title('|JSA| = |\alpha\Phi|');

    nexttile; imagesc(ws_THz, wi_THz, best.JSI); set(gca,'YDir','normal'); colorbar;
    xlabel('\omega_s (THz)'); ylabel('\omega_i (THz)'); title('JSI = |JSA|^2');

    % Schmidt spectrum
    nexttile; stem(1:numel(best.lambda), best.lambda, 'filled'); grid on;
    xlabel('Mode index n'); ylabel('\lambda_n');
    title(sprintf('Schmidt spectrum — K=%.3f, Pur=%.3f', best.K, best.pur));
    xlim([0 max(10,nnz(best.lambda>1e-4))+1]);

    % First two idler/signal modes (magnitudes)
    u1 = best.U(:,1); v1 = best.V(:,1); 
    u2 = best.U(:,2); v2 = best.V(:,2);
    nexttile; plot(wi_THz, abs(u1)/max(abs(u1)), 'LineWidth', 1.4); grid on;
    xlabel('\nu_i (THz)'); ylabel('|u_1|'); title('Idler u_1');
    nexttile; plot(ws_THz, abs(v1)/max(abs(v1)), 'LineWidth', 1.4); grid on;
    xlabel('\nu_s (THz)'); ylabel('|v_1|'); title('Signal v_1');
    nexttile; hold on;
    plot(wi_THz, abs(u2)/max(abs(u2)), 'LineWidth', 1.4);
    plot(ws_THz, abs(v2)/max(abs(v2)), '--', 'LineWidth', 1.4);
    grid on; xlabel('THz'); ylabel('|u_2| / |v_2|');
    legend('Idler u_2','Signal v_2'); title('Second modes');

    sgtitle(t, sprintf(['BEST L = %.3f mm  (λ_p=%.1f nm, λ_s=%.1f nm, ' ...
                        'λ_i=%.1f nm, θ≈%.1f°;  Eq(8) L_{opt}≈%.3f mm)'], ...
        Lmm, lp_user, ls_user, li_user, theta_deg, Lopt_mm), ...
        'FontWeight','bold');

    if save_jsa_png
        tag = sprintf('BEST_JSA_Schmidt_Alt%.3f_Anc%.3f_L%.3fmm_lp%.1f_ls%.1f.png', ...
                      Alt_um0, Anc_um0, Lmm, lp_user, ls_user);
        exportgraphics(gcf, fullfile(jsa_outdir, tag), 'Resolution', 300);
    end

end % geometry loop
end
