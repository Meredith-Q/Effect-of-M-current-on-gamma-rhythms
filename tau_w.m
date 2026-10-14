%function tau_w=tau_w(v);
%tau_w_max = 400; %NEED TO VARY %M-current gating variable time constant
%tau_w=tau_w_max./(3.3*exp((v+35)/20)+exp(-(v+35)/20));

function tau_w=tau_w(v);
global tau_w_max
tau_w=tau_w_max./(3.3*exp((v+35)/20)+exp(-(v+35)/20));