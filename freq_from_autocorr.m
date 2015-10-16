% freq_from_autocorr() - computes the auto-correlation function of a signal
% and extracts the associated second peak, as an estimative of its
% fundamental frequency. A quadratic interpolation for estimating the true
% position of an inter-sample maximum when nearby samples are known was
% used
% 
% Usage:
%    >> fundamental_freq = freq_from_autocorr(signal, filters, fs)
% 
% Inputs:
%   dataset - input data matrix (row or column vector)
%   filters - low and high cutoff values of the band-pass filtering of EEG data {[ low_cutoff high_cutoff ]}
%   fs      - sampling rate
%     
% Outputs:
%  fundamental_freq - estimative of the signal fundamental frequency
% 
% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function fundamental_freq = freq_from_autocorr(signal, filters, fs)

if ~isempty(filters)
    signal = eegfilt(signal, fs, filters(1), 0); % freq_max = 45 Hz
    signal = eegfilt(signal, fs, 0, filters(2));
end

corr = conv(signal, fliplr(signal)); %convolution between the signal and the associated time-inverted
[ ~, l ] = max(corr);
corr = corr(l:end);

d = diff(corr);
start = find(d > 0); start = start(1);

[ ~, peak ] = max(corr(start:end));
peak = peak + start;

[ px, ~ ] = parabolic(corr, peak); % perform quadric interpolation for more accurate estimation

fundamental_freq = fs / px;

end
