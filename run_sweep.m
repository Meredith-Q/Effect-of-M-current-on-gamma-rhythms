%run_sweep.m
clear; clc; close all;
global tau_w_max

g_m_list = [0 0.1 0.25 0.5 0.75];
tau_list = [50 100 200 400 800];
seed_list = [5489];          % baseline only. Use [5489 1 2] for 3 network realizations
nG = length(g_m_list); nT = length(tau_list); nS = length(seed_list);

% results arrays: (g_M, tau, seed)
MPC = nan(nG,nT,nS); GB = nan(nG,nT,nS); CC = nan(nG,nT,nS);
NF  = nan(nG,nT,nS);   % network (E population) frequency, Hz
EF  = nan(nG,nT,nS);   % mean E cell frequency, Hz
IFR = nan(nG,nT,nS);   % mean I cell frequency, Hz

for isd = 1:nS
for ia = 1:nG
for ib = 1:nT
    g_m = g_m_list(ia);
    tau_w_max = tau_list(ib);
    try
        gamma_simulator_M_Current_PING;
        PINGSynchronyMeasures;
        MPC(ia,ib,isd) = mean_mpc;
        GB(ia,ib,isd)  = B;
        CC(ia,ib,isd)  = mean_crcorr;
        NF(ia,ib,isd)  = avnfreq;
        EF(ia,ib,isd)  = mean(avecfreq);
        IFR(ia,ib,isd) = mean(avicfreq);
        if isd == 1   % save rasters for the first seed only
            saveas(1, sprintf('raster_gm%g_tau%g.png', g_m, tau_w_max));
        end
    catch ME
        fprintf('Run failed (g_m=%g, tau=%g, seed %d): %s (line %d of %s)\n', ...
            g_m, tau_w_max, isd, ME.message, ME.stack(1).line, ME.stack(1).name);
    end
    close all;   % avoid piling up figures
end
end
end

save sweep_results.mat g_m_list tau_list seed_list MPC GB CC NF EF IFR

% ---------- sanity check: g_M = 0 row should not depend on tau ----------
row0 = find(g_m_list == 0);
if ~isempty(row0)
    fprintf('\ng_M = 0 check (should be ~0 if tau is irrelevant):\n');
    spread = @(x) max(x) - min(x);   % not range(): the simulator defines a variable named range
    fprintf('  max deviation across tau: MPC %.3g, B %.3g, CC %.3g\n', ...
        spread(MPC(row0,:,1)), spread(GB(row0,:,1)), spread(CC(row0,:,1)));
end

% ---------- summary table ----------
[GG,TT] = meshgrid(g_m_list, tau_list);
M = @(A) reshape(mean(A,3,'omitnan')', [], 1);   % mean over seeds, ordered like meshgrid
T = table(GG(:), TT(:), M(MPC), M(GB), M(CC), M(NF), M(EF), M(IFR), ...
    'VariableNames', {'g_M','tau_w_max','MPC','B','CC','NetFreq_Hz','Efreq_Hz','Ifreq_Hz'});
writetable(T, 'sweep_summary.csv');

% ---------- color plots ----------
measures = {MPC, GB, CC, NF, EF, IFR};
titles   = {'Mean phase coherence','Golomb bursting measure B','Mean cross-correlation', ...
            'Network frequency (Hz)','Mean E cell frequency (Hz)','Mean I cell frequency (Hz)'};
for q = 1:numel(measures)
    figure;
    imagesc(mean(measures{q},3,'omitnan')); 
    cb = colorbar; ylabel(cb, titles{q}); colormap(jet);
    set(gca,'YDir','normal','XTick',1:nT,'XTickLabel',tau_list, ...
            'YTick',1:nG,'YTickLabel',g_m_list);
    xlabel('\tau_{w,max} (ms)'); ylabel('g_M');
    title(titles{q});
    saveas(gcf, sprintf('colorplot_%d.png', q));
end