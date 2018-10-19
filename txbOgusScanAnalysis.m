%% TxB OGUS Scan Analysis (Low and High Energy Pulse Data)
% Author: Morgan Roberts
% Date: 19/10/2018
%
% DESCRIPTION
%     A script that extracts data from a raw oscilliscope data file,
%     creates a time series starting at zero, bandpass filters the voltage
%     signal with given parameters, pads the signal to centre the impulse
%     response in the time series, applies a Gaussian window function with
%     given parameters, and plots the amplitude spectrum of the impulse
%     response. There are two available data sets, which contain the
%     impulse response of 'TxB' to a low energy, and high energy, optically
%     generated ultrasound source (OGUS) pulse. The script also calculates
%     both the full-width-half-maximum bandwidth as a percentage and the
%     centre frequency, and plots these for reference.
%
% REFERENCES
%    [1] applyFilter is from the k-wave Toolbox, and was written by Ben
%    Cox and Bradley Treeby. 
%    [2] getWin, spect, and fwhm are from the k-wave Toolbox, and were
%    written by Bradley Treeby.
%
%% Section 1: Defining parameters and reading files.

close all
TEMPERATURE = 23.4; % [degC]

% Set to 1 for high energy pulse data, and 0 for low energy pulse data.
energy = 1;

% Choosing which data to load.
if energy == 1
    scope_data = 'TxB_OGUS_140msflp_19_10_18.xlsx';
elseif energy == 0
    scope_data = 'TxB_OGUS_263msflp_19_10_18.xlsx';
end

% Extracting data from excel file.
t          = xlsread(scope_data, 'A3:A4');
trig_volt  = xlsread(scope_data, 'B3:B2002');
volt       = xlsread(scope_data, 'C3:C2002');

% Calculating sample period and frequency.
dt = t(2) - t(1);
fs = 1 / dt;

% Creating time array, starting at zero.
time = 0:dt:dt * (length(volt) - 1);

%% Section 2: Applying filters.

% Setting passband cutoff frequencies (60dB attentuation in stop band).
band_pass = [0.5e6, 5e6];

% Filtering voltage data.
filt_volt = applyFilter(volt, fs, band_pass, 'BandPass');

%% Section 3: Padding to align impulse with centre of window function.

% Finding maximum of impulse response.
index = find(abs(filt_volt) == max(abs(filt_volt)));

% Forward-Padding to align maximum with the centre.
pad_volt = [zeros(1, length(filt_volt) - (2 * index) - 1), filt_volt];
pad_time = 0:dt:dt * ((2 * (length(filt_volt) - 1)) - (2 * index));

% Applying Gaussian window function with std dev parameter.
filt_param = 0.05;
win        = getWin(length(pad_time), 'Gaussian', 'Param', filt_param);
win_volt   = win' .* pad_volt;

%% Section 4: Amplitude Sprectrum analysis.

[freq, amp_spect] = spect(win_volt, 1/dt, 'FFTLength', length(win_volt) * ...
                                                    8, 'PowerTwo', 'true');
                                                   
%% Section 5: Finding Bandwidth. (Adapted from fwhm.m, see References)

% Find the maximum value of the amplitude spectrum.
[f_max, i_max] = max(amp_spect);

% Setup the indexing variables for finding the leading edge.
k1 = i_max;

% Looping until the index at half maximum is found.
while amp_spect(k1) > (0.5 * f_max)
    k1 = k1 - 1;
    if k1 < 1
        error('Left half maximum not found.');
    end
end

% Linear interpolation between the previous values to find the frequency
% at the leading edge at half maximum.
m = (amp_spect(k1+1) - amp_spect(k1)) / (freq(k1+1) - freq(k1));
c = amp_spect(k1) - freq(k1) * m;
i_leading = (0.5 * f_max  - c) / m;

% Setup the indexing variables for finding the trailing edge.
k2 = i_max;

% Looping until the index at half maximum is found.
while amp_spect(k2) > (0.5 * f_max)
    k2 = k2 + 1;
    if k2 > length(amp_spect)
        error('Right half maximum not found.');
    end
end

% Linear Interpolation between the previous values to find the frequency
% at the trailing edge at half maximum.
m = (amp_spect(k2 - 1) - amp_spect(k2)) / (freq(k2 - 1) - freq(k2));
c = amp_spect(k2) - freq(k2) * m;
i_trailing = (0.5 * f_max  - c) / m;
    
% Computing the FWHM.
fwhm_val = abs(i_trailing - i_leading);

% Finding the centre frequency of the FWHM.
centre_freq = mean([i_trailing, i_leading]);

% Calculating and displaying the bandwidth as a percentage.
bandwidth = (fwhm_val / centre_freq) * 100;
fprintf('The centre frequency is %.2fMHz.\nThe bandwidth is %.f%%',...
                                            centre_freq / 1e6, bandwidth);

%% Section 6: Plotting.

% Setting axis colors.
fig1 = figure(1);
left_color  = [0, 0, 0];
right_color = [0, 0, 0];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);

% Plotting Original Signal.
subplot(2,2,1);
plot(time, volt, 'k');
ylabel('Voltage [V]');
xlabel('Time [s]');
%xlim(xlim_param);
title('Original signal');

% Plotting BandPass Filtered Signal.
subplot(2,2,2);
plot(time, filt_volt, 'k');
ylabel('Voltage [V]');
xlabel('Time [s]');
%xlim(xlim_param);
title(sprintf('Bandpass Filtered signal %.2fMHz - %.0fMHz', ...
                                  band_pass(1) / 1e6, band_pass(2) / 1e6));

% Plotting Padded and Windowed Signal.
subplot(2,2,3);
yyaxis left
plot(pad_time, win_volt, 'k');
ylabel('Voltage [V]');
yyaxis right
plot(pad_time, win, 'b');
xlabel('Time [s]');
title(sprintf('Front-Padded and Gaussian Windowed (std dev = %.2f)', ...
                                                              filt_param));

% Plotting Amplitude Spectrum.
subplot(2,2,4);
hold on
plot(freq / 1e6, amp_spect, 'k');
h = plot([i_leading, i_trailing] / 1e6, [0.5 * f_max, 0.5 * f_max], 'k--');
g = plot(centre_freq / 1e6, amp_spect(k1 + ((k2 - k1) / 2)), 'kx ');
xlabel('Frequency [MHz]');
ylabel('FFT Amplitude');
title('Amplitude Spectrum');
xlim([0, 5]);
legend([h, g], {sprintf('FWHM Bandwidth = %.f%%', bandwidth), ...
               sprintf('Centre Frequency = %.2fMHz%', centre_freq / 1e6)});
hold off