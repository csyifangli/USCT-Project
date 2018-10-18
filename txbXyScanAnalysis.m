%% TxB X-Y Scan Analysis
% Original Author: Elly Martin
% Updated by: Morgan Roberts
% Date: 15/10/2018
%
% DESCRIPTION
%     Script to extract and save measured data, apply hydrophone sensitivity,
%     extract mag and phase of pressure at fundamental frequency then project
%     the field back to the transducer face.
%     required toolboxes: k-wave-matlab, bug-measurement-toolbox


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
    'Scanning\Scan Data\TxB\xy scan 2 10th oct'];

% Setting filename to give the extracted data.
FILENAME = 'Extract_TxB_xy_scan2_10_10_18.mat';

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

% Extracting sample period, poiunt spacing, and water temperature.
dt          = ScanData.samplePeriod(1);
dx          = ScanData.PointSpacing(1);
temperature = mean(ScanData.Temperature);

% Calculating the number of data points per cycle, and the number of 
% cycles contained in the data.
pts_cycle = 1 / (TRANSDUCER_FREQ * dt);
n_cycles  = round(size(ScanData.Voltage, 2) / pts_cycle);

% Finding the index to cut the data at to get integer cycles.
data_cut = round(n_cycles * pts_cycle);
if data_cut > length(size(ScanData.Voltage, 2))
    n_cycles = n_cycles - 1;
    data_cut = round(n_cycles * pts_cycle);
end

%% Section 4: Extracting the pressure data from voltage data.

% Extract the hydrophone sensitivity information.
[sensitivity, ~, ~, ~] = getHydrophoneSensitivity(HYDROPHONE, SENS_PATH);

% Extracting the band pass filtered pressure data, using the frequency 
% response of the measurement device as defined by the sensitivity input.
% Waveforms are first filtered in the frequency domain using the given 
% filter parameters.
measured_p = applyCalibration(ScanData.Voltage(:,1:data_cut), 1/dt, ...
                        sensitivity, 'Dim', 2, 'FilterParam', FILTER_FREQ);

% Extracting the amplitude and phase of the pressure data at the transducer
% driving frequency from the frequency spectrum, which is calculated using 
% a windowed and zero padded FFT.                    
[mag, phase] = extractAmpPhase(measured_p, 1/dt, TRANSDUCER_FREQ, 'Dim',...
                                                       2, 'FFTPadding', 3);
%%
% make the input pressure into a complex 2D plane: 
input_pressure = mag .* exp(1i.* phase);
% reshape to make 2D

input_pressure = reshape(input_pressure, ScanData.NumberPoints);
X = reshape(ScanData.posX, ScanData.NumberPoints);
Y = reshape(ScanData.posY, ScanData.NumberPoints);

% Plotting magnitude pressure map on 2-D measurement plane.
figure(1);
subplot(1,2,1);
surf(X, Y, abs(input_pressure), 'LineStyle', 'None');
grid off
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title(sprintf('TxB Pressure Magnitude at z = %.1fmm', ScanData.posZ(1)));
axis image;
view(2);
c = colorbar('eastoutside');
c.Label.String = 'Pressure Magnitude [Pa]';


% calculate Nz based on the point spacing dx, and the distance to the
% transducer - you know this already, or you can extract it from the posZ
% field in ScanData.
Nz = floor(ScanData.posZ(1)/dx);

c0 = speedSoundWater(temperature);
% calculate the pressure at (or near) the transducer face:
pressure = angularSpectrumCW(input_pressure, dx/1000, Nz, TRANSDUCER_FREQ, c0, 'Reverse', true);
 
% the output will be a complex 2D plane. Have a look at both the magnitude
% and phase (mag and angle) and plot.

%figure(2);
subplot(1,2,2);
hold on
surf(X, Y, abs(pressure(:,:,1)), 'LineStyle', 'None');
t=0:0.01:2*pi;
% Drawing Reference circles with radii
xPZT = 19*sin(t);
yPZT = 19*cos(t);
xPLA = 21*sin(t);
yPLA = 21*cos(t);
xQWML = 16.65*sin(t);
yQWML = 16.65*cos(t);
h1 = plot3(xPLA, yPLA, max(max(abs(pressure(:,:,1))))*ones(size(xPLA)), 'k');
h2 = plot3(xPZT, yPZT, max(max(abs(pressure(:,:,1))))*ones(size(xPZT)), 'k--');
h3 = plot3(xQWML, yQWML, max(max(abs(pressure(:,:,1))))*ones(size(xQWML)), 'k-.');
hold off
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title('TxB Pressure Magnitude projected back to transducer face');
axis image;
legend([h1, h2, h3], {'Housing Boundary', 'PZT Boundary', 'QWML Boundary'});
c = colorbar('eastoutside');
c.Label.String = 'Pressure Magnitude [Pa]';

