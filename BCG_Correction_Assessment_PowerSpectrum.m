% BCG_Correction_Assessment_PowerSpectrum() - perform BCG artefact correction
% quantification as the average spectral power within artefact and
% physiological representative windows
% 
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function [ power, power_bkg ] = BCG_Correction_Assessment_PowerSpectrum...
    (N, art_harm, eeg_bcg, NFFT, L, harmonics, w, f_deep, rel_conv)

% FREQUENCY DOMAIN
power = zeros(N, art_harm);
power_bkg = zeros(N, art_harm);

for n = 1:N
    y = eeg_bcg(n, :);
    Y = fft(y, NFFT) / L;
    PowerSpectrum = 2 * abs(Y(1:NFFT/2 + 1));
    
    for h = 1:art_harm
        f_samples = harmonics(n, h);
        win = (f_samples - w:f_samples + w);
        win_deep = f_deep(win); true_win = win(win_deep == 1);
        
        W = PowerSpectrum(true_win);
        % Normalize the power according to the window size
        power(n, h) = sum(W) / length(W);
        
        win = (f_samples + ceil(5 * w)):(f_samples + ceil(11 * w));
        win_deep = f_deep(win); true_win = win(win_deep == 1);
        
        W_bkg = PowerSpectrum(true_win);
        % Normalize the power according to the window size
        power_bkg(n, h) = sum(W_bkg) / length(W_bkg);
    end
end

power = power .* rel_conv; power_bkg = power_bkg .* rel_conv;