function [b1_label,b1_control,g_label,g_control]=pCASL_pulses(g_max, dt, gamma, pw_rfperf, rfperf_sep, a_gzperf, gpav, b1pav, tlbl, tag_dist)

% INPUTS

% g_max: maximum gradient amplitude (G/cm)
% dt: simulation time-step (us)
% gamma: proton gyromagnetic ratio (Hz/G)
% pw_rfperf: duration of the Hanning RF pulse (us)
% rfperf_sep: time between RF pulses (us)
% a_gzperf: amplitude of the gradient (G/cm)
% gpav: time average of the gradient (G/cm)
% b1pav: time average of the RF (G)
% tlbl: label duration (us)
% tag_dist: distance from isocentre to labelling plane (cm)

% OUTPUTS

% b1_label: RF for label
% b1_control: RF for control
% g_label: gradient for label
% g_control: gradient for control

% Written by Stephen Wastling October 2022

% area of the gradient lobe that accompanies the RF pulse 
pw_gzperf = pw_rfperf;
area_pgrad = pw_rfperf * a_gzperf; 

% timings of a refocussing gradient with the -area_pgrad
pw_gzperfr = area_pgrad / g_max;
a_gzperfr = -g_max;

% amplitude of the refocussing gradient during the labelling such that the time average gradient is gpav
a_gzperfr1 = (1 - (gpav * rfperf_sep) / area_pgrad) * a_gzperfr;

% amplitude of the Hanning RF pulse (G)
a_rfperf = 2 * b1pav * (rfperf_sep / pw_rfperf); 

% number of pulses that will fit in the PCASL pulse train
nbl1 = floor(tlbl ./ rfperf_sep); 

%% generate the gradient waveforms

% gradient simultaneous with the RF pulse
gz_perf = ones(1, pw_gzperf ./ dt) .* a_gzperf;

% refocussing pulse in the control state
gz_perfr = ones(1, pw_gzperfr ./ dt) .* a_gzperfr;

% refocussing pulse in the labelling state
gz_perfr1 = ones(1,pw_gzperfr ./ dt) .* a_gzperfr1;

tpad_gz = rfperf_sep - (pw_rfperf + pw_gzperfr);
gz_pad = zeros(1,tpad_gz ./ dt);

% build the gradients by concatenating the different sections together
g_c=[gz_perf gz_perfr(1 : end) gz_pad]; 
g_control=repmat(g_c,1,nbl1);

g_l=[gz_perf gz_perfr1(1 : end) gz_pad]; 
g_label=repmat(g_l,1,nbl1);

%% generate the RF pulse

N = pw_rfperf / dt;
n=0 : N - 1;
b1_single_pulse = a_rfperf .* 0.5 .* (1 - cos(2 .* pi .* n ./ (N - 1)));

tpad_rf = rfperf_sep - pw_rfperf;
rf_pad = zeros(1, tpad_rf ./ dt);
b1_single_pulse = [b1_single_pulse rf_pad];

b1_label=repmat(b1_single_pulse, 1, nbl1);
b1_control_temp=[b1_single_pulse -b1_single_pulse];
b1_control = repmat(b1_control_temp, 1, nbl1 / 2);

% shift in excitation k-space needed required to move the tagging plane 
% away from the isocentre
% using eqns. 5.25 and 5.37 in Handboook of MRI Pulse Sequences by Bernstein
kz_l = gamma .* fliplr(cumsum(fliplr(g_label))) .* dt .* 1E-6; 
kz_c = gamma .* fliplr(cumsum(fliplr(g_control))) .* dt .* 1E-6; 

b1_label = b1_label .* exp(-2 * pi * 1i * tag_dist * kz_l); 
b1_control = b1_control .* exp(-2 * pi * 1i * tag_dist * kz_c); 
