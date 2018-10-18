%% TxB X-Y Scan Analysis
% Original Author: Elly Martin
% Updated by: Morgan Roberts
% Date: 15/10/2018
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
%     with angular restriction described inZeng, X., & McGough, R. J. 
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

% Extracting sample period, point spacing, and water temperature.
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
                                                   
%% Section 5: Plotting complex pressure in the measurement plane.

% Converting the input pressure into a complex 2D plane. 
input_pressure = mag .* exp(1j .* phase);

% Reshaping the 1-D pressure and position arrays into 2-d grid.
input_pressure = reshape(input_pressure, ScanData.NumberPoints);
X              = reshape(ScanData.posX,  ScanData.NumberPoints);
Y              = reshape(ScanData.posY,  ScanData.NumberPoints);

% Plotting pressure magnitude in the measurement plane.
figure(1);
subplot(1,2,1);
surf(X, Y, abs(input_pressure), 'LineStyle', 'None');
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title(sprintf('TxB Pressure Magnitude at z = %.1fmm', ScanData.posZ(1)));
axis image;
view(2);
c = colorbar('eastoutside');
c.Label.String = 'Pressure Magnitude [Pa]';

% Plotting pressure phase in the measurement plane.
subplot(1,2,2);
surf(X, Y, angle(input_pressure), 'LineStyle', 'None');
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title(sprintf('TxB Pressure Phase at z = %.1fmm', ScanData.posZ(1)));
axis image;
view(2);
c = colorbar('eastoutside');
c.Label.String = 'Pressure Phase [rad]';

%% Section 6: Projecting pressure data back to the transducer face.

% Calculating number of backward projections needed using point spacing.
Nz = floor(ScanData.posZ(1) / dx);

% Calculating the sound speed in water at the given temperature
c0 = speedSoundWater(temperature);

% Projecting the 2-D measurement plane onto a 3D volume using the angular 
% spectrum method. 
pressure = angularSpectrumCW(input_pressure, dx/1000, Nz, ...
                                     TRANSDUCER_FREQ, c0, 'Reverse', true);
 
%% Section 7: Plotting complex pressure at the transducer face.

% Plotting pressure magnitude at the transducer face.
figure(2);
subplot(1,2,1);
hold on
surf(X, Y, abs(pressure(:,:,1)), 'LineStyle', 'None');

% Drawing reference circles with radii matched to transducer geometry.
t = 0:0.01:(2 * pi);

% PZT Element Boundary.
x_pzt = 19*sin(t);
y_pzt = 19*cos(t);

% Transducer Housing Boundary.
x_pla = 21*sin(t);
y_pla = 21*cos(t);

% Quarter-Wavelength-Matching-Layer Boundary.
x_qwml = 16.65*sin(t);
y_qwml = 16.65*cos(t);

% Plotting circles.
h1 = plot3(x_pla,  y_pla,  max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_pla)),  'k');
h2 = plot3(x_pzt,  y_pzt,  max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_pzt)),  'k--');
h3 = plot3(x_qwml, y_qwml, max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_qwml)), 'k-.');
hold off

% Adding axes labels, title, legend, and colorbar.
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title('TxB Pressure Magnitude projected back to transducer face');
axis image;
legend([h1, h2, h3], {'Housing Boundary', 'PZT Boundary', ...
                                                         'QWML Boundary'});
c = colorbar('eastoutside');
c.Label.String = 'Pressure Magnitude [Pa]';

% Plotting pressure phase at the transducer face.
subplot(1,2,2);
hold on
surf(X, Y, angle(pressure(:,:,1)), 'LineStyle', 'None');

% Plotting circles.
h1 = plot3(x_pla,  y_pla,  max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_pla)),  'k');
h2 = plot3(x_pzt,  y_pzt,  max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_pzt)),  'k--');
h3 = plot3(x_qwml, y_qwml, max(max(abs(pressure(:,:,1)))) * ...
                                                ones(size(x_qwml)), 'k-.');
hold off

% Adding axes labels, title, legend, and colorbar.
xlabel('X Position [mm]');
ylabel('Y Position [mm]');
title('TxB Pressure Phase projected back to transducer face');
axis image;
legend([h1, h2, h3], {'Housing Boundary', 'PZT Boundary', ...
                                                         'QWML Boundary'});
c = colorbar('eastoutside');
c.Label.String = 'Pressure Phase [rad]';

