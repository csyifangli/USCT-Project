scope_data = 'TxB_OGUS_140msflp_19_10_18.xlsx';

% Extracting data from excel file.
t          = xlsread(scope_data, 'A3:A4');
trig  = xlsread(scope_data, 'B3:B2002');
volt       = xlsread(scope_data, 'C3:C2002');

% Calculating sample period and frequency.
dt = t(2) - t(1);
fs = 1 / dt;

% Creating time array, starting at zero.
time = 0:dt:dt * (length(volt) - 1);

% Setting passband cutoff frequencies (60dB attentuation in stop band).
band_pass = [0.5e6, 5e6];

% Filtering voltage data.
filt_volt = applyFilter(volt, fs, band_pass, 'BandPass');
filt_trig= applyFilter(trig, fs, band_pass, 'BandPass');

figure(1);
plot(time, volt, 'k', time, trig / 100, 'b');

[r, lags] = xcorr(volt, gradient(trig));
figure(2);
plot(lags, abs(r));

tof_samples = lags(abs(r) == max(abs(r)));
tof         = tof_samples * dt; 

