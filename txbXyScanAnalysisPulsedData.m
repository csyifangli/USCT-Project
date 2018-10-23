%% TxB X-Y Scan Analysis (Pulse Data)
% Author: Morgan Roberts
% Date: 18/10/2018
%
% DESCRIPTION
%     A script to extract and save measured data, apply hydrophone 
%     sensitivity, extract the magnitude and phase of the pressure at 
%     the transducer driving frequency then project the field back to the 
%     transducer face. The pressure magnitude and phase at the driving
%     frequency are plotted in the measurement plane, and at the transducer
%     face.
%
% REFERENCES
%     [1] extract_time_series, and getHydrophoneSensitivity are from the 
%     BUG Measurement Toolbox, and are written by Elly Martin.
%     [2] applyCalibration is from the BUG Measurement Toolbox and the 
%     k-wave Toolbox, and is written by Elly Martin and Bradley Treeby.
%     [3] extractAmpPhase is from the k-wave Toolbox, and is written by 
%     Bradley Treeby and Yan To Ling.
%     [4] speedSoundWater is part of the k-wave Toolbox, and is written by
%     Bradley Treeby. It calculates the speed of sound in distilled water 
%     at a given temperature using the 5th order polynomial given by 
%     Marczak (1997) "Water as a standard in the measurements of speed of 
%     sound in liquids," J. Acoust. Soc. Am., 102, 2776-2779. 
%     [5] angularSprectrumCW is part of the k-wave Toolbox, and is written
%     by Bradley Treeby. The implementation follows the spectral propagator
%     with angular restriction described in: Zeng, X., & McGough, R. J. 
%     (2008). Evaluation of the angular spectrum approach for simulations 
%     of near-field pressures. The Journal of the Acoustical Society of 
%     America, 123(1), 68-76.

%% Section 1: Defining parameters, file names, and file paths.

% Scan and processing parameters.
TRANSDUCER_FREQ = 2e6;                 % [Hz]
FILTER_FREQ     = [0.5e6, 5e6, 0.5e6]; % [Hz] (start, end, taper width)
HYDROPHONE      = 'PA2042';            % 0.2mm Needle Hydrophone #2042

% Setting directories containing hydrophone sensitivity data, scan text
% files and scan '.ssf' file.
SENS_PATH       = ['C:\Users\BUG\Synology Drive\Biomedical Ultrasound '...
                            'Group\Equipment\Hydrophones\PA needle 2042\'];
INPUT_DATA_DIR = ['C:\Users\BUG\Synology Drive\Morgan Synology Drive\'...
                                'Scanning\Scan Data\TxB\xy scan 10th oct'];

% Setting filename to give the extracted data.
FILENAME = 'Extract_TxB_xy_scan_10_10_18';

%% Section 2: Determining whether data needs to be extracted or loaded.

% Set to 1 if data needs to be extracted. Set to zero if data needs to be
% processed and analysed.
EXTRACT = 0;

if EXTRACT == 1
    ScanData = extract_time_series(INPUT_DATA_DIR, INPUT_DATA_DIR, ...
                                             FILENAME, TRANSDUCER_FREQ, 0);
elseif EXTRACT == 0
    load(FILENAME)
end

%% Section 3: Truncating the CW data to a whole number of cycles.

% Extracting sample period, point spacing, and water temperature.
dt          = ScanData.samplePeriod(1);
dx          = ScanData.PointSpacing(1);
temperature = mean(ScanData.Temperature);

%% Section 4: Extracting the pressure data from voltage data.

% Extract the hydrophone sensitivity information.
[sensitivity, ~, ~, ~] = getHydrophoneSensitivity(HYDROPHONE, SENS_PATH);

% Extracting the band pass filtered pressure data, using the frequency 
% response of the measurement device as defined by the sensitivity input.
% Waveforms are first filtered in the frequency domain using the given 
% filter parameters.
measured_p = applyCalibration(ScanData.Voltage, 1/dt, ...
                        sensitivity, 'Dim', 2, 'FilterParam', FILTER_FREQ);

% Reshaping the 1-D pressure and position arrays into 2-D grid.
measured_p = reshape(measured_p,     ScanData.NumberPoints(1), ...
                            ScanData.NumberPoints(1), size(measured_p, 2));
X          = reshape(ScanData.posX,  ScanData.NumberPoints);
Y          = reshape(ScanData.posY,  ScanData.NumberPoints);

avgs = 50
rand_pos_x = ceil(rand(1,avgs) * 26) + 32;
rand_pos_y = ceil(rand(1,avgs) * 26) + 32;
amp_spect = zeros(avgs, 17773);
win_amp_spect = zeros(avgs, 17773);
figure(1);
xlim([0 , 1e7]);
for i = 1:avgs
    pressure     = permute(measured_p(rand_pos_x(i), rand_pos_y(i), :), ...
                                                                [3, 2, 1]);
    index        = find(abs(pressure) == max(abs(pressure)));
    pad_pressure = [zeros(1, length(pressure) - (2 * index) - 1), pressure'];
    time         = 0:dt:dt * ((2 * (size(measured_p, 3) - 1)) - (2 * index));
    win          = getWin(length(pressure), 'Gaussian', 'Param', 0.5);
    win_pressure = win .* pressure;
    [f, func_as, func_ps] = spect(pressure, 1/dt, 'FFTLength', ...
                                     length(pressure) * 8, 'PowerTwo', 'true');
    [f_win, func_as_win] = spect(win_pressure, 1/dt, ...
                        'FFTLength', length(pressure) * 8, 'PowerTwo', 'true');
    fprintf('Iteration %1f, length of fft is %1f\n', i, length(f_win));
    amp_spect(i,:) = func_as;
    win_amp_spect(i,:) = func_as_win;
    plot(f, win_amp_spect, f, amp_spect)
    drawnow
    xlim([0, 1e7]);
    hold on
end


avg_amp_spect = mean(amp_spect, 1);
avg_win_amp_spect = mean(win_amp_spect, 1);

plot(f, avg_amp_spect, 'k', f, avg_win_amp_spect, 'b');
xlim([0, 1e7]);


OUTPUT_DATA_DIR = ['C:\Users\BUG\Synology Drive\Morgan Synology Drive\'...
                               'Scanning\Scan Analysis\TxB\TxB OGUS Data'];
AmpSpectXyScan = struct('frequency', f, 'amplitudeSpectrum', ...
                                                            avg_amp_spect);
save([OUTPUT_DATA_DIR, '/', 'AmpSpectXyScan', '.mat'], ...
                                            'AmpSpectXyScan', '-v7.3');

% 
% figure(1);
% subplot(1,2,1);
% yyaxis left
% plot(time, pressure, 'k', time, win_pressure, 'r');
% ylabel('Pressure [Pa]');
% yyaxis right
% plot(time, win, 'b')
% xlabel('Time [s]');
% ylabel('Window Amplitude [au]');
% 
% subplot(1,2,2);
% plot(f, func_as, 'k', f, func_as_win, 'b');
% xlim([0 1e7]);
% xlabel('Frequency [Hz]');
% ylabel('FFT Amplitude')
% 

% Find the maximum value of the amplitude spectrum.
[f_max, i_max] = max(avg_amp_spect);

% Setup the indexing variables for finding the leading edge.
k1 = i_max;

% Looping until the index at half maximum is found.
while avg_amp_spect(k1) > (0.5 * f_max)
    k1 = k1 - 1;
    if k1 < 1
        error('Left half maximum not found.');
    end
end

% Linear interpolation between the previous values to find the frequency
% at the leading edge at half maximum.
m = (avg_amp_spect(k1+1) - avg_amp_spect(k1)) / (f(k1+1) - f(k1));
c = avg_amp_spect(k1) - f(k1) * m;
i_leading = (0.5 * f_max  - c) / m;

% Setup the indexing variables for finding the trailing edge.
k2 = i_max;

% Looping until the index at half maximum is found.
while avg_amp_spect(k2) > (0.5 * f_max)
    k2 = k2 + 1;
    if k2 > length(avg_amp_spect)
        error('Right half maximum not found.');
    end
end

% Linear Interpolation between the previous values to find the frequency
% at the trailing edge at half maximum.
m = (avg_amp_spect(k2 - 1) - avg_amp_spect(k2)) / (f(k2 - 1) - f(k2));
c = avg_amp_spect(k2) - f(k2) * m;
i_trailing = (0.5 * f_max  - c) / m;
    
% Computing the FWHM.
fwhm_val = abs(i_trailing - i_leading);

% Finding the centre frequency of the FWHM.
centre_freq = mean([i_trailing, i_leading]);

% Calculating and displaying the bandwidth as a percentage.
bandwidth = (fwhm_val / centre_freq) * 100;
fprintf('The centre frequency is %.2fMHz.\nThe bandwidth is %.f%%\n',...
                                            centre_freq / 1e6, bandwidth);