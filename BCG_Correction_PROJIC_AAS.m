% BCG_Correction_PROJIC_AAS() - perform BCG artefact correction by correction the 
%                               BCG-related ICs selected by PROJIC for the artefact occurrences using
%                               AAS [Allen et al., NeuroImage 1998], before back-projecting all the ICs in the reconstruction of the EEG signal
% Usage:
%    >> [ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC(dataset, ...
%       activations, mix_matrix, limits_ecg, filters, k_clusters, art_harm, harm_thr, ...
%       TR, win_hz, plt, nb_plts)
% 
% Inputs:
%   dataset     - input data structure with the following mandatory fields:
%                  1) data -> EEG signal (each row represent an EEG channel)
%                  2) ecg -> ECG signal (row vector)
%                  3) srate [double] -> sampling rate
%                  4) Kp -> R peak annotations
%                  5) Cne [cell] -> name of each EEG channel
%                  6) TR -> annotations defining the begining of each volume acquisition
%   activations - IC time-courses returned by a given ICA decomposition
%                 algorithm
%   mix_matrix  - mixing matrix returned by a given ICA decomposition algorithm
%                (typically, activations = mix_matrix * dataset.data)
%   limits_ecg  - upper and lower limits, in milliseconds, for epoching the BCG artefact 
%                 occurrences using the R peak annotations as triggers {[ lower_lim upper_lim ]} 
%   filters     - low and high cutoff values of the band-pass filtering of EEG data {[ low_cutoff high_cutoff ]}
%   n_win       - number of windows (can be an array, and BCG correction 
%                 will be performed for all number of windows specified, using AAS [Allen et al., NI 1998])
%   k_clusters  - number of clusters (can be an array, and BCG correction
%                 will be performed for all number of clusters specified, upon IC
%                 selection using PROJIC.m) -> if n_win > 1 and k_clusters > 1, 
%                 all possible combinations of parameters will be used for BCG artefact correction
%   art_harm    - number of cardiac frequency harmonics for BCG artefact correction evaluation {recommended = 5}
%   harm_thr    - threshold for inclusion of harmonics {recommended = 0.25; for higher variability across subjects and channels, the threshold should be lower}
%   TR          - repetition time of the fMRI acquisition [s]
%   win_hz      - window length for which the BCG artefact correction will be assessed (recommended = 0.065} [Hz]
%   delay_QRS   - time delay between the R peak annotations and the BCG
%                 artefact occurrences {recommended: 0.21, according to Allen et. al, NeuroImage 1998} [s]
%   plt         - 1/0 = do/don't display the results
%   nb_plts     - number of subplots within each plot for displaying the selected IC time-courses
%     
% Outputs:
%   PowerFreq [cell]     - BCG artefact correction quantification in terms of
%                          relative reduction of average spectral power within artefact
%                          representative windows. It has (n_win, k_clusters + 1) elements, the first
%                          representing this quantification before correction
%   PowerFreq_bkg [cell] - the same for PowerFreq, but storing the BCG 
%                          artefact correction quantification in terms of
%                          relative reduction of average spectral power within 
%                          physiological background representative windows
% 
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC_AAS(dataset, ...
    activations, mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, ...
    harm_thr, TR, win_hz, delay_QRS, plt, nb_plts)

warning off

delay_QRS = round(delay_QRS * dataset.srate); % convert sec to number of samples

if ~isempty(filters) % filter the EEG data
    eeg = eegfilt(dataset.data, dataset.srate, filters(1), 0);
    eeg = eegfilt(eeg, dataset.srate, 0, filters(2));
end

nw = length(n_win);
nk = length(k_clusters);
PowerFreq = cell(3, nw, nk + 1); % 1->power, 2->param, 3->ICs
PowerFreq_bkg = cell(3, nw, nk + 1); % 1->power, 2->param, 3->ICs

R = dataset.Kp; % retrive the R peak annotations 
N = size(activations, 1);

% Frequency domain "parameters"
Fs = dataset.srate; L = length(eeg);  
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = Fs/2 * linspace(0, 1, NFFT/2 + 1);

tr_f = (1:floor(f(end) / (1 / TR)) - 1) .* (1 / TR); % remove from the BCG artefact correction quantification TR-related frequencies
tr_fm = tr_f - 0.02; tr_fM = tr_f + 0.02;

f_deep = ones(1, length(f));
for i = 1:length(tr_f)
    f1 = find(f < tr_fm(i)); f1_s = f1(end);
    f2 = find(f < tr_fM(i)); f2_s = f2(end);
    f_deep(f1_s:f2_s) = 0;
end

cardiac_f = 1 / (L / Fs / length(R)); % cardiac frequency
w_aux = find(f < win_hz); w = w_aux(end); % window length
[ max_conv, harmonics ] = harmonics_Lorentz(dataset, eeg, cardiac_f, Fs, art_harm, 0.1, 0.3, 0.01, win_hz, 0); % fine-tune the position of the spectral peaks
rel_conv = (max_conv ./ max(max(max_conv))) > harm_thr; % criterion for the inclusion of harmonics in the quantification

% Apply ICA-based BCG correction
for n = 1:nw
    for k = 1:(nk + 1)
        act_copy = activations;
        
        if k < 2 % quantification before correction
            eeg_bcg = eeg; ICs = []; param = [];
        else
            unmix_matrix = mix_matrix^(-1); 

            [ ~, ~, ~, ICs ] = PROJIC(dataset, eeg, activations, mix_matrix, ...
                limits_ecg, k_clusters(k - 1), plt, nb_plts); % call PROJIC for selecting the BCG-related ICs
            
            act_corrected = fmrib_pas_EEGActivations(activations, sort(ICs), R, 'mean', [], n_win(n), Fs, delay_QRS);
            act_copy(sort(ICs), :) = act_corrected; % correct the selected ICs for artefact occurrences using AAS
            
            eeg_bcg = unmix_matrix * act_copy; param = [ k_clusters(k - 1) n_win(n) ]; % back-reconstruct the EEG signal using all the ICs, including the AAS-corrected ones
        end
        
        [ power, power_bkg ] = BCG_Correction_Assessment_PowerSpectrum... % quantify the average power within artefact and physiological background representative windows
            (N, art_harm, eeg_bcg, NFFT, L, harmonics, w, f_deep, rel_conv);
        
        PowerFreq{1, n, k} = power; PowerFreq{2, n, k} = param; PowerFreq{3, n, k} = ICs;
        PowerFreq_bkg{1, n, k} = power_bkg; PowerFreq_bkg{2, n, k} = param; PowerFreq_bkg{3, n, k} = ICs;
    end
end