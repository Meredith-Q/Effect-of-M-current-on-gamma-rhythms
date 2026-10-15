clear;

params0_MCurrent;

dt05=dt/2; m_steps=floor(t_final/dt+0.1); r_e=exp(-dt05/tau_d_stoch_e); r_i=exp(-dt05/tau_d_stoch_i);

t_discard = 50; % ms

g_m_vals = (0:0.05:1)';
avecfreq_vals = zeros(size(g_m_vals));
avnfreq_vals  = zeros(size(g_m_vals));
avicfreq_vals = zeros(size(g_m_vals));
mpc_vals      = zeros(size(g_m_vals));
B_vals        = zeros(size(g_m_vals));
crcorr_vals   = zeros(size(g_m_vals));

for gi = 1:length(g_m_vals)
    g_m = g_m_vals(gi);

    rand('twister',5489); randn('state',0);

    g_ee=g_hat_ee*ones(num_e,num_e).*(sign(p_ee-rand(num_e,num_e))+1)/2/(p_ee*num_e);
    g_ei=g_hat_ei*ones(num_e,num_i).*(sign(p_ei-rand(num_e,num_i))+1)/2/(p_ei*num_e);
    g_ie=g_hat_ie*ones(num_i,num_e).*(sign(p_ie-rand(num_i,num_e))+1)/2/(p_ie*num_i);
    g_ii=g_hat_ii*ones(num_i,num_i).*(sign(p_ii-rand(num_i,num_i))+1)/2/(p_ii*num_i);

    v_e=-80+rand(num_e,1)*30; n_e=0*ones(num_e,1); m_e=m_e_inf(v_e); h_e=0*ones(num_e,1); s_e=zeros(num_e,1);
    v_i=-80+rand(num_i,1)*30; n_i=0*ones(num_i,1); m_i=m_i_inf(v_i); h_i=0*ones(num_i,1); s_i=zeros(num_i,1);
    s_stoch_e=zeros(num_e,1); s_stoch_i=zeros(num_i,1);
    num_spikes_e=0; num_spikes_i=0;
    w_e = zeros(num_e,1);

    t_e_spikes=[]; i_e_spikes=[];
    t_i_spikes=[]; i_i_spikes=[];

    for k=1:m_steps
        t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;

        v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +g_m*w_e.*(-100-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e) ...
               +I_e(t_old)+g_stoch_e*s_stoch_e.*(v_rev_e-v_e);
        n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
        h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
        w_e_inc=(w_inf(v_e)-w_e)./tau_w(v_e);
        s_e_inc=0.5*(1+tanh(v_e/4)).*(1-s_e)./tau_r_e-s_e./tau_d_e;

        v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ei'*s_e).*(v_rev_e-v_i)+(g_ii'*s_i).*(v_rev_i-v_i) ...
               +I_i(t_old)+g_stoch_i*s_stoch_i.*(v_rev_e-v_i);
        n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
        h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
        s_i_inc=0.5*(1+tanh(v_i/4)).*(1-s_i)./tau_r_i-s_i./tau_d_i;

        v_e_tmp=v_e+dt05*v_e_inc;
        n_e_tmp=n_e+dt05*n_e_inc;
        m_e_tmp=m_e_inf(v_e_tmp);
        h_e_tmp=h_e+dt05*h_e_inc;
        w_e_tmp=w_e+dt05*w_e_inc;
        s_e_tmp=s_e+dt05*s_e_inc;
        v_i_tmp=v_i+dt05*v_i_inc;
        n_i_tmp=n_i+dt05*n_i_inc;
        m_i_tmp=m_i_inf(v_i_tmp);
        h_i_tmp=h_i+dt05*h_i_inc;
        s_i_tmp=s_i+dt05*s_i_inc;
        s_stoch_e=s_stoch_e*r_e; s_stoch_i=s_stoch_i*r_i;

        v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +g_m*w_e_tmp.*(-100-v_e_tmp) ...
               +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp)+(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +I_e(t_mid)+g_stoch_e*s_stoch_e.*(v_rev_e-v_e_tmp);
        n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
        h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
        w_e_inc=(w_inf(v_e_tmp)-w_e_tmp)./tau_w(v_e_tmp);
        s_e_inc=0.5*(1+tanh(v_e_tmp/4)).*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e;
        v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp)+35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ei'*s_e_tmp).*(v_rev_e-v_i_tmp)+(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +I_i(t_mid)+g_stoch_i*s_stoch_i.*(v_rev_e-v_i_tmp);
        n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
        h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
        s_i_inc=0.5*(1+tanh(v_i_tmp/4)).*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;

        v_e_old=v_e; v_i_old=v_i;

        v_e=v_e+dt*v_e_inc; n_e=n_e+dt*n_e_inc; m_e=m_e_inf(v_e);
        h_e=h_e+dt*h_e_inc; w_e=w_e+dt*w_e_inc; s_e=s_e+dt*s_e_inc;
        v_i=v_i+dt*v_i_inc; n_i=n_i+dt*n_i_inc; m_i=m_i_inf(v_i);
        h_i=h_i+dt*h_i_inc; s_i=s_i+dt*s_i_inc;

        s_stoch_e=s_stoch_e*r_e; s_stoch_i=s_stoch_i*r_i;
        u_e=rand(num_e,1); u_i=rand(num_i,1);
        s_stoch_e=s_stoch_e+max(sign(f_stoch_e*dt/1000-u_e),0).*(1-s_stoch_e);
        s_stoch_i=s_stoch_i+max(sign(f_stoch_i*dt/1000-u_i),0).*(1-s_stoch_i);

        which_e=find(v_e_old<0 & v_e>=0); which_i=find(v_i_old<0 & v_i>=0);
        l_e=length(which_e); l_i=length(which_i);

        if l_e>0
            range=[num_spikes_e+1:num_spikes_e+l_e];
            i_e_spikes(range)=which_e;
            t_e_spikes(range)=(v_e(which_e)*(k-1)*dt-v_e_old(which_e)*k*dt)./(v_e(which_e)-v_e_old(which_e));
            num_spikes_e=num_spikes_e+l_e;
        end
        if l_i>0
            range=[num_spikes_i+1:num_spikes_i+l_i];
            i_i_spikes(range)=which_i;
            t_i_spikes(range)=(v_i(which_i)*(k-1)*dt-v_i_old(which_i)*k*dt)./(v_i(which_i)-v_i_old(which_i));
            num_spikes_i=num_spikes_i+l_i;
        end
    end

    spiketimes_e=[t_e_spikes', i_e_spikes'];
    spiketimes_i=[t_i_spikes', i_i_spikes'];

    spiketimes_e2=spiketimes_e(spiketimes_e(:,1)>t_discard,:);
    spiketimes_i2=spiketimes_i(spiketimes_i(:,1)>t_discard,:);

    % ---- 放电频率 ----
    nisis=diff(spiketimes_e2(:,1)); nisis=nisis(nisis>0.1);
    if ~isempty(nisis), avnfreq_vals(gi)=1/(mean(nisis)/1000); end

    avecfreq=zeros(num_e,1);
    for j=1:num_e
        jspikes=spiketimes_e2(spiketimes_e2(:,2)==j,:);
        if ~isempty(jspikes)
            jisis=diff(jspikes(:,1));
            avecfreq(j)=1/(mean(jisis)/1000);
        end
    end
    avecfreq_vals(gi)=mean(avecfreq);

    avicfreq=zeros(num_i,1);
    for j=1:num_i
        jspikes=spiketimes_i2(spiketimes_i2(:,2)==j,:);
        if ~isempty(jspikes)
            jisis=diff(jspikes(:,1));
            avicfreq(j)=1/(mean(jisis)/1000);
        end
    end
    avicfreq_vals(gi)=mean(avicfreq);

    % measure
    if size(spiketimes_e2,1) > num_e
        [mean_mpc, ~] = mpc_PINGnetwork(num_e, spiketimes_e2);
        [~, traces, traces_all] = spiketraces_quiet(num_e, spiketimes_e2);

        sigma = zeros(num_e,1);
        for i=1:num_e
            sigma(i) = var(traces(i,:));
        end
        sigma_all = var(traces_all);
        B = sigma_all/mean(sigma);

        [mean_crcorr, ~] = crcorr_network(num_e, traces);

        mpc_vals(gi)    = mean_mpc;
        B_vals(gi)      = B;
        crcorr_vals(gi) = mean_crcorr;
    else
        mpc_vals(gi) = NaN; B_vals(gi) = NaN; crcorr_vals(gi) = NaN;
    end

    fprintf('g_m=%.3f | avnfreq=%.2f Hz | MPC=%.3f | B=%.3f | crcorr=%.3f\n', ...
        g_m, avnfreq_vals(gi), mpc_vals(gi), B_vals(gi), crcorr_vals(gi));
end

% plot
figure;
subplot(2,2,1); hold on;
plot(g_m_vals, avecfreq_vals, '-o','DisplayName','Ave_{EF}');
plot(g_m_vals, avicfreq_vals, '-s','DisplayName','Ave_{IF}');
plot(g_m_vals, avnfreq_vals,  '-^','DisplayName','Overall');
xlabel('g_M'); ylabel('F [Hz]'); title('Firing Frequency vs g_M'); legend show; grid on;

subplot(2,2,2);
plot(g_m_vals, mpc_vals, '-o','LineWidth',2,'Color',[0.1 0.4 0.8]);
xlabel('g_M'); ylabel('Mean Phase Coherence'); title('MPC vs g_M'); ylim([0 1.05]); grid on;

subplot(2,2,3);
plot(g_m_vals, B_vals, '-o','LineWidth',2,'Color',[0.85 0.1 0.1]);
xlabel('g_M'); ylabel('B'); title('Golomb Synchrony Measure B vs g_M'); grid on;

subplot(2,2,4);
plot(g_m_vals, crcorr_vals, '-o','LineWidth',2,'Color',[0.1 0.6 0.3]);
xlabel('g_M'); ylabel('Mean cross-correlation'); title('Pairwise Cross-Correlation vs g_M'); grid on;

%spiketraces
function [timevec,traces,traces_all] = spiketraces_quiet(n,spiketimes)
    srate=10; min_timevec = floor(min(spiketimes(:,1))); max_timevec = ceil(max(spiketimes(:,1)));
    sigma=2; peak=1;
    traces = zeros(n,(max_timevec-min_timevec)*srate+1);
    for i = 1:n
        stimes=spiketimes(spiketimes(:,2)==i,1);
        [trace,timevec,~]=spikegauss_quiet(stimes,srate,min_timevec,max_timevec,sigma,peak);
        traces(i,:) = trace';
    end
    traces_all = mean(traces);
end

function [spkvec,timevec,updatedpeak]=spikegauss_quiet(timestamps,srate,min_timevec,max_timevec,sigma,peak)
    timevec = min_timevec : 1/srate : max_timevec;
    spike_count = histc( timestamps, timevec );
    gk_mu=0; gk_sigma=sigma;
    gk_x = -10*sigma+1/srate : 1/srate : +10*sigma;
    gk = 1/(sqrt(2*pi)*gk_sigma) * exp( -(gk_x-gk_mu).^2 / (2*gk_sigma^2) );
    if peak==0, gk=gk/sum(gk); else, gk=peak*gk/max(gk); end
    updatedpeak=max(gk);
    spkvec=conv(spike_count, gk, 'same');
end