% PROJIC() -  perform BCG-related IC selection by projecting BCG
%             artefact-triggered EEG data onto the IC space and clustering 
%             them as having the highest-powered projections. 
% Usage:
%    >> [ P, power, ERP, ICs ] = PROJIC(dataset, activations, wghts, sph, ...
%                                limits_ecg, filters, k_cluster, plt, nb_plts)
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
%                occurrences using the R peak annotations as triggers {[ lower_lim upper_lim ]} 
%   filters     - low and high cutoff values of the band-pass filtering of EEG data 
%                 {[ low_cutoff high_cutoff ]}
%   k_clusters  - number of clusters
%   plt         - 1/0 = do/don't display the results
%   nb_plts     - number of subplots within each plot for displaying the selected IC time-courses
%     
% Outputs:
%   P     - the projections of the BCG artefact-triggered EEG data by each IC
%   power - the power of each projection as the squared sum of the projection instances
%   ERP   - the BCG artefact-triggered EEG data
%   ICs   - the selected BCG-related ICs
%
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ P, power, ERP, ICs ] = PROJIC(dataset, eeg, activations, ...
    mix_matrix, limits_ecg, k_cluster, plt, nb_plts)

warning off

R = dataset.Kp; % retrive the R peak annotations

cardiac_R = epoch(eeg, R, limits_ecg .* (dataset.srate / 1000), 'verbose', 'off'); % epoch EEG data
ERP = MeanEpochs(cardiac_R); N = size(ERP, 1); % average the epoched EEG across events

P = mix_matrix * ERP; % compute the projections
power = sum(P .^ 2, 2); % compute the power of each projection
power_perc = int8((power ./ sum(power)) .* 100);

% perform clustering of the ICs
[ cl, cn ] = kmeans(power, k_cluster, 'Replicates', 5000); % replicate 5000 times to ensure stable solution
[ ~, cn_ind ] = sort(cn, 'descend');

ICs = 1:N; ICs = ICs(ismember(cl, cn_ind(1:k_cluster - 1)));
[ ~, power_id ] = sort(power, 'descend');
ICs = power_id(1:length(ICs)); % find the selected BCG-related ICs

% information for plot purposes
[ ~, m ] = max(sum(ERP .^ 2, 2));
minS = limits_ecg(1) / 1000;
time = (1:length(ERP)) ./ dataset.srate;
time_onset = (time + minS) * 1000;

if plt % DISPLAY RESULTS
    
    % Average ERP for R peaks as time-locking events
    figure('Name', 'ICs');
    plot(time_onset, ERP');
    hold on
    a = plot(time_onset, ERP(m, :));
    set(a, 'LineWidth', 3);
    axis('tight');
    xlabel('Time [s]');
    ylabel('Amplitude [\mu V]');
    hl = legend(a, [ 'Channel #', num2str(m) ]);
    hold off
    
    % Final Result: Projections and percentage of explained variance
    b = cell(1, length(ICs));
    if length(ICs) > 1
        for i = 1:length(ICs)
            b{i} = ['IC #', num2str(ICs(i))];
        end
    else
        b{1} = ['IC #', num2str(ICs(1))];
    end
    
    figure('Name', 'ICs');
    
    subplot(1,2,1)
    plot(time_onset, P');
    hold on
    a = plot(time_onset, P(ICs, :)');
    set(a, 'LineWidth', 3);
    axis('tight');
    xlabel('Time [ms]', 'FontSize', 12);
    ylabel('Amplitude [a.u.]', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    hl = legend(a, b);
    hold off
    
    subplot(1,2,2)
    plot(power_perc, '-o');
    xlabel('ICs', 'FontSize', 12);
    ylabel('Relative Power [%]', 'FontSize', 12);
    axis('tight');
    hold on
    plot(ICs, power_perc(ICs), 'ro');
    set(gca, 'FontSize', 12);
    % set(gca, 'YTick', 1:2:19);
    hold off
    
    % Time-course of the BCG-related components
    
    l = length(ICs); T = (1:length(activations)) ./ dataset.srate;
    count = 1;
    figure('Name', 'ICs', 'Color', 'white'), hold on
    for i = 1:l
        if count <= nb_plts
            subplot(nb_plts, 1, count), plot(T, activations(ICs(i), :));
            ylabel([ 'IC#', num2str(ICs(i)) ]);
            xlabel('Time [s]');
            line([ T(R); T(R) ], [ min(activations(ICs(i), :)) max(activations(ICs(i), :)) ], 'Color', 'red');
            axis('tight');
            
            count = count + 1;
        else
            count = 1;
            
            hold off
            figure('Name', 'ICs', 'Color', 'white'), hold on
            subplot(nb_plts, 1, count), plot(T, activations(ICs(i), :));
            ylabel([ 'IC#', num2str(ICs(i)) ]);
            xlabel('Time [s]');
            line([ T(R); T(R) ], [ min(activations(ICs(i), :)) max(activations(ICs(i), :)) ], 'Color', 'red');
            axis('tight');
            
            count = count + 1;
        end
    end
end