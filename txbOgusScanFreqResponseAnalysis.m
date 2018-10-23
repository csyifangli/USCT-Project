%% TxB Frequency Response Further Analysis

load('AmpSpectLowEnergy.mat');
load('AmpSpectHighEnergy.mat');
load('AmpSpectXyScan.mat');

L = 0.001; % [m]
c = 4010;  % [m/s]

freq_low       = AmpSpectLowEnergy.frequency;
amp_spect_low  = AmpSpectLowEnergy.amplitudeSpectrum;
freq_high      = AmpSpectHighEnergy.frequency;
amp_spect_high = AmpSpectHighEnergy.amplitudeSpectrum;
freq_scan      = AmpSpectXyScan.frequency;
amp_spect_scan = AmpSpectXyScan.amplitudeSpectrum;

% Normalize

k = max([amp_spect_low, amp_spect_high, amp_spect_scan]);
k_amp_spect_scan = 1 / max(amp_spect_scan);
k_amp_spect_low = 1 / max(amp_spect_low);
k_amp_spect_high = 1 / max(amp_spect_high);

figure(1);
plot(freq_low / 1e6,  amp_spect_low  * k_amp_spect_low,  'k',   ...
     freq_high / 1e6, amp_spect_high * k_amp_spect_high, 'k--', ...
     freq_scan / 1e6, amp_spect_scan * k_amp_spect_scan, 'b-.');
 xlabel('Frequency [MHz]');
 xlim([0, 4]);
 ylabel('Normalized Amplitude Spectrum');
 legend('Frequency Response to Low Energy OGUS Pulse', ...
        'Frequency Response to High Energy OGUS Pulse', ...
        ['Amplitude Spectrum averaged over X-Y measurement plane for '...
                                                     'pulsed excitation']); 