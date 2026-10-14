% compute network synchrony measures

% compute mean phase coherence
[mean_mpc, mpc_cellpairs] = mpc_PINGnetwork(num_e,spiketimes_e2);
mean_mpc

% construct spiketime traces with spkgauss.m
[timevec,traces, traces_all] = spiketraces(num_e,spiketimes_e2);

% compute Golomb bursting measure
sigma = zeros(num_e,1);
for i=1:num_e
    sigma(i) = var(traces(i,:));
end
sigma_all = var(traces_all);
B = sigma_all/mean(sigma)

% compute pairwise cross correlations with zero time lag
[mean_crcorr, crcorr_cellpairs] = crcorr_network(num_e,traces);
mean_crcorr
