% harmonics_Lorentz() - determines the exact value of the fundamental
% frequency and first 'art_harm' harmonics of the ECG signal, f_n, by convolving
% a Lorentzian function with power spectrum windows centered in f_n, and extracting 
% the convolution peaks (Rothlübbers et al., Brain Topography, 2014), for
% each EEG channel and harmonic
% 
% Usage:
%    >> [ max_conv, real_harmonics ] = harmonics_Lorentz(eeg, f0, fs, art_harm, ...
%     PeakRadiusf, LorentzRf, LorentzSigma, win_hz, plt)
% 
% Inputs:
%   dataset      - input data structure with the following mandatory fields:
%                  1) data -> EEG signal (each row represent an EEG channel)
%                  2) ecg -> ECG signal (row vector)
%                  3) srate [double] -> sampling rate
%                  4) Kp -> R peak annotations
%                  5) Cne [cell] -> name of each EEG channel
%                  6) TR -> annotations defining the begining of each volume acquisition
%   f0           - approximation of the cardiac fundamental frequency (based on the
%       heart rate, for example); if empty, f0 is estimated by determining
%       the second peak of the autocorrelation function of the ECG signal
%   fs           - sampling rate
%   art_harm     - number of cardiac frequency harmonics for BCG artefact correction evaluation {recommended = 5}
%   PeakRadiusf  - length of the power spectrum windows {recommended for BCG = 0.1} [Hz]
%   LorentzRf    - radius of the Lorentzian function {recommended for BCG = 0.3} [Hz]
%   LorentzSigma - width of the Lorenzian function {recommended for BCG = 0.01} [Hz]
%   win_hz       - window length for which the BCG artefact correction will be assessed (recommended = 0.065} [Hz]
%   plt          - 1/0 = do/don't display the results
%     
% Outputs:
%   max_conv     - maximum value of the convolution peak, across channels
%   and harmonics
%   preal_harmonics - the exact values of f_n for each EEG channel and
%   harmonic
% 
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ max_conv, real_harmonics ] = harmonics_Lorentz(dataset, eeg, f0, fs, ...
    art_harm, PeakRadiusf, LorentzRf, LorentzSigma, win_hz, plt)

ecg = dataset.ecg; % retrieve the ECG signal
[ N, L ] = size(eeg); real_harmonics = zeros(N, art_harm);

if isempty(f0) % estimate f0 from the auto-correlation function, after band-pass filtering the ECG between 1 and 45 Hz
    f0 = freq_from_autocorr(ecg, [ 1 45 ], fs);
end

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = fs/2 * linspace(0, 1, NFFT/2 + 1);

f0_index = find(f <= f0); f0_index = f0_index(end);
harmonics_index = (1:art_harm) .* f0_index;

harmonics = (1:art_harm) .* f0;
max_conv = zeros(N, length(harmonics));

for n = 1:N
    channel = eeg(n, :);
    Y = fft(channel, NFFT) / L; 
    FFTData = 2 * abs(Y(1:NFFT/2 + 1)); % compute the FFT for each EEG channel
    
    df = diff(f); df = df(1);
    PeakSearchR = ceil(PeakRadiusf / df);
    
    lorentzR = ceil(LorentzRf / df); %% create a lorentzian peak
    lorentzX = (-lorentzR:lorentzR) * df;
    lorentz  = LorentzSigma^2 ./ (LorentzSigma^2 + (lorentzX).^2); % build the Lorentzian function according to inputs
    
    if plt
        f_10Hz = find(f < 10); f_10Hz = f_10Hz(end);
        figure('Name', 'PSD'), plot(f(1:f_10Hz), FFTData(1:f_10Hz)), hold on
    end
    
    for p = 1:length(harmonics) %% find optimal convolution of lorentz with given peaks
        
        artRemPeak = FFTData(harmonics_index(p) + (-PeakSearchR:PeakSearchR)); % local spctrum around peak
        c = conv(artRemPeak,lorentz, 'same'); % convolve with lorentzian
        optimalShift = find(c == max(c)); max_conv(n, p) = max(c); % find optimal shift
        maximumPeakIndex = harmonics_index(p) - PeakSearchR + optimalShift - 1;
        real_harmonics(n, p) = maximumPeakIndex; % determine the exact values of f_n
        
        if plt 
            w_aux = find(f < win_hz); w = w_aux(end); % window length
            
            f_samples = real_harmonics(n, p);
            f_hz = horzcat(f(f_samples - w), f(f_samples + w));
            f_hz_bkg = horzcat(f(f_samples + ceil(5 * w)), f(f_samples + ceil(11 * w)));
            
            plot(f(real_harmonics(n, p)), FFTData(real_harmonics(n, p)), 'ro', 'MarkerFaceColor', 'r');
            line([ f_hz; f_hz ], [ min(FFTData) max(FFTData) ], 'Color', 'red', 'LineWidth', 1);
            line([ f_hz_bkg; f_hz_bkg ], [ min(FFTData) max(FFTData) ], 'Color', 'black', 'LineWidth', 1);
        end
    end
end



