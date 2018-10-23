%% Backing and Matching Layer Composition Analysis
% Author: Morgan Roberts
% Date: 19/10/18
%
% DESCRIPTION
%
% REFERENCES
%     Half wavelength mode
%     https://www.americanpiezo.com/images/stories/content_images/pdf/apc_materials_properties.pdf
%     https://universe.ida.dk/meetupfiles/downloadfile/?meetupNumber=311638&fileid=2cf310a7-cbb8-48dd-ba40-a297bd0fdb5d
%     http://usglobalimages.stratasys.com/Main/Files/Material_Spec_Sheets/MSS_PJ_PJMaterialsDataSheet.pdf?v=635785205440671440
%     https://arxiv.org/ftp/arxiv/papers/1005/1005.4304.pdf
%% Section 1: Defining Parameters.

PZT_MATERIALS       = {'840', '841', '850', '854', '855', '880', 'Pz27'};
FREQUENCY_CONSTANTS = [2005, 2005, 2040, 2000, 2079, 2110, 1950];
PZT_DENSITIES       = [7600, 7600, 7600, 7600, 7600, 7600, 7700];
pzt_frequencies     = containers.Map(PZT_MATERIALS,FREQUENCY_CONSTANTS);
pzt_densities       = containers.Map(PZT_MATERIALS,PZT_DENSITIES);

pzt_material = '840';
c_pzt        = pzt_frequencies(pzt_material) * 2; % [m/s]
rho_pzt      = pzt_densities(pzt_material);       % [kg/m3]
Z_pzt        = c_pzt * rho_pzt;                   % [Rayls]

HOUSING_MATERIALS   = {'PLA', 'VeroBlack'};
SOUND_SPEEDS        = {2279, 2495};
HOUSING_DENSITIES   = {1170, 1175};
housing_sound_speed = containers.Map(HOUSING_MATERIALS,SOUND_SPEEDS);
housing_density     = containers.Map(HOUSING_MATERIALS,HOUSING_DENSITIES);

housing_material = 'VeroBlack';
c_housing        = housing_sound_speed(housing_material); % [m/s]
rho_housing      = housing_density(housing_material);     % [kg/m3]
Z_housing        = c_housing * rho_housing;               % [Rayls] 

%% Section 2: Extract time-of-flight data from different composite samples.
%
%% Section 3: Calculate sound speeds for samples and store.
%

temperature      = 23.4;                         % [degC]
c_water          = speedSoundWater(temperature); % [m/s]
sample_thickness = 0.005;                        % [m]
tof_without      = 35.13e-6;                     % [s]
tof_with         = 33e-6;                        % [s]

c_sample = (c_water * sample_thickness) / (sample_thickness - ...
                                     (c_water * (tof_without - tof_with)));


%% Section 4: Choose desired reflection coefficients at boundaries.
%

Z_qwml = sqrt(Z_pzt * Z_housing);
%Z_qwml    = 2.6e6;
Z_backing = 5980*2551;

pzt_qwml_R    = ((Z_qwml    - Z_pzt) / (Z_qwml    + Z_pzt))^2;
pzt_backing_R = ((Z_backing - Z_pzt) / (Z_backing + Z_pzt))^2;

%% Section 5: Calculate required impedance values.
%

Z_qwml = sqrt(Z_housing * Z_pzt);
Z_backing = Z_pzt;

%% Section 6: Interpolate material properties data to decide compositions.
%

% Data to collect for each sample

% TOF without sample  |  TOF with sample  |  Sample thickness  | Sample mass  |  Sample Volume
% ============================================================================================
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |
%                     |                   |                    |              |