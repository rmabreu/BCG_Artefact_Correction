% MAIN FUNCTION - EEG BCG Correction
%
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

% change the following parameters accordingly
art_harm = 5; % number of harmonics to include in the BCG artefact correction quantification
win_hz = 0.065; % window length for which the BCG artefact correction will be assessed
filters = [ 0.5 45 ]; % low and high cutoff values of the band-pass filtering of EEG data
TR = 2.5; % repetition time of the fMRI acquisition [s]
harm_thr = 0.25; % threshold for inclusion of harmonics {recommended = 0.25; for higher variability across subjects and channels, the threshold should be lower}
plt = 0; % 1/0 = do/don't display the results
nb_plts = 4; % number of subplots within each plot for displaying the selected IC time-courses

% define the path in which data is stored
data_path = '/home/rmabreu/Data/MATLAB/Scripts/BCG_Artefact_Correction/NI_Code/';

% load all necessary matrices (data, activations and mixing matrix from ICA)
dt = load([ data_path, 'data.mat' ]); dataset = dt.dataset; clear dt;
sph = load([ data_path, 'sphere' ], '-ascii');
wghts = load([ data_path, 'weights' ], '-ascii');
activations = load([ data_path, 'activations' ], '-ascii');

mix_matrix = wghts * sph;

% retrieve R peak annotations and compute appropriate limits for epoching
% EEG data using the R peaks as triggers (dependent on the subjects' heart
% rate)
R = dataset.Kp; m_R = (min(diff(R)) / dataset.srate) * 1000;
lim_inf_ms = -100;
if m_R < 4 * abs(lim_inf_ms)
    lim_sup_ms = lim_inf_ms + (4 * abs(lim_inf_ms));
else
    lim_sup_ms = lim_inf_ms + m_R;
end
limits_ecg = [ lim_inf_ms lim_sup_ms ];

% define the standard time-delay (210 ms) between R peak and BCG artefact
% occurrences (according to Allen et al., NeuroImage 1998)
delay_QRS = 0.21;

% PROJIC
k_clusters = 2:5; 
tic; fprintf('\n--->>>--->BCG Artefact Correction using PROJIC...\n\n');
[ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC(dataset, activations, ...
    mix_matrix, limits_ecg, filters, k_clusters, art_harm, harm_thr, TR, win_hz, ...
    plt, nb_plts);

mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
if mttoc < 60
    fprintf('--->>>--->BCG Artefact Correction using PROJIC finished in %d min %d sec.\n', mttoc, sttoc);
else
    httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
    fprintf('--->>>--->BCG Artefact Correction using PROJIC finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
end

% PROJIC-OBS
k_clusters = 2:5; npc = 3:6; 
tic; fprintf('\n--->>>--->BCG Artefact Correction using PROJIC-OBS...\n\n');
[ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC_OBS(dataset, activations, ...
    mix_matrix, limits_ecg, filters, npc, k_clusters, art_harm, harm_thr, ...
    TR, win_hz, delay_QRS, plt, nb_plts);

mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
if mttoc < 60
    fprintf('--->>>--->BCG Artefact Correction using PROJIC-OBS finished in %d min %d sec.\n', mttoc, sttoc);
else
    httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
    fprintf('--->>>--->BCG Artefact Correction using PROJIC-OBS finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
end

% PROJIC-AAS
k_clusters = 2:5; n_win = 10:10:40; 
tic; fprintf('\n--->>>--->BCG Artefact Correction using PROJIC-AAS...\n\n');
[ PowerFreq, PowerFreq_bkg ] = BCG_Correction_PROJIC_AAS(dataset, activations, ...
    mix_matrix, limits_ecg, filters, n_win, k_clusters, art_harm, harm_thr, ...
    TR, win_hz, delay_QRS, plt, nb_plts);

mttoc = floor(toc/60); sttoc = round(toc - mttoc*60);
if mttoc < 60
    fprintf('--->>>--->BCG Artefact Correction using PROJIC-AAS finished in %d min %d sec.\n', mttoc, sttoc);
else
    httoc=floor(mttoc / 60); mttoc=round(mttoc - httoc*60);
    fprintf('--->>>--->BCG Artefact Correction using PROJIC-AAS finished in %d hrs %d min %d sec.\n', httoc, mttoc, sttoc);
end