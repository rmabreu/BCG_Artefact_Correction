% BCG_Correction_PROJIC() - perform BCG artefact correction by removing the 
%                           contribution of the BCG-related ICs selected by 
%                           PROJIC in the back-reconstruction step of the EEG signal
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
%   k_clusters  - number of clusters (can be an array, and BCG correction
%                 will be performed for all number of clusters specified, upon IC
%                 selection using PROJIC.m)
%   art_harm    - number of cardiac frequency harmonics for BCG artefact correction evaluation {recommended = 5}
%   harm_thr    - threshold for inclusion of harmonics {recommended = 0.25; for higher variability across subjects and channels, the threshold should be lower}
%   TR          - repetition time of the fMRI acquisition [s]
%   win_hz      - window length for which the BCG artefact correction will be assessed (recommended = 0.065} [Hz]
%   plt         - 1/0 = do/don't display the results
%   nb_plts     - number of subplots within each plot for displaying the selected IC time-courses
%     
% Outputs:
%   PowerFreq [cell]     - BCG artefact correction quantification in terms of
%                          relative reduction of average spectral power within artefact
%                          representative windows. It has (k_clusters + 1) elements, the first
%                          representing this quantification before correction
%   PowerFreq_bkg [cell] - the same for PowerFreq, but storing the BCG 
%                          artefact correction quantification in terms of
%                          relative reduction of average spectral power within 
%                          physiological background representative windows
% 
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC(dataset, ...
    activations, mix_matrix, limits_ecg, filters, k_clusters, art_harm, harm_thr, ...
    TR, win_hz, plt, nb_plts)

warning off

if ~isempty(filters) % filter the EEG data
    eeg = eegfilt(dataset.data, dataset.srate, filters(1), 0);
    eeg = eegfilt(eeg, dataset.srate, 0, filters(2));
end

nk = length(k_clusters);
PowerFreq = cell(3, nk + 1); % 1->power, 2->param, 3->ICs
PowerFreq_bkg = cell(3, nk + 1); % 1->power, 2->param, 3->ICs

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
for k = 1:(nk + 1)
    if k < 2 % quantification before correction
        eeg_bcg = eeg; ICs = []; param = [];
    else % quantification after correction, for different number of clusters
        unmix_matrix = mix_matrix^(-1); 

        [ ~, ~, ~, ICs ] = PROJIC(dataset, eeg, activations, mix_matrix, ...
            limits_ecg, k_clusters(k - 1), plt, nb_plts); % call PROJIC for selecting the BCG-related ICs
        
        I = 1:size(activations, 1); ICs_remove = ~ismember(I, ICs);
        Z = diag(ICs_remove);
        
        eeg_bcg = unmix_matrix * Z * activations; param = k_clusters(k - 1); % remove the ICs from the back-reconstruction step of the EEG signal
    end
    
    [ power, power_bkg ] = BCG_Correction_Assessment_PowerSpectrum... % quantify the average power within artefact and physiological background representative windows
        (N, art_harm, eeg_bcg, NFFT, L, harmonics, w, f_deep, rel_conv);
    
    % store the results
    PowerFreq{1, k} = power; PowerFreq{2, k} = param; PowerFreq{3, k} = ICs;
    PowerFreq_bkg{1, k} = power_bkg; PowerFreq_bkg{2, k} = param; PowerFreq_bkg{3, k} = ICs;
end

