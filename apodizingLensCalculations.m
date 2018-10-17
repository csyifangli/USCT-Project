%% Apodizing Lens Calculations
% Author: Morgan Roberts
% Date: 17/10/2018
%
% DESCRIPTION:
%
% REFERENCES:
%
%speedSoundWater calculates the speed of sound in distilled water at a
%     a given temperature using the 5th order polynomial given by Marczak
%     (1997) "Water as a standard in the measurements of speed of sound in
%     liquids," J. Acoust. Soc. Am., 102, 2776-2779.
%
%% Section 1: Defining Parameters.

RING_DIAMETER    = 0.2;  % [m]
TEMPERATURE      = 22;   % [degC]
LENS_SOUND_SPEED = 2495; % [m/s]

% Calculating the speed of sound in water for given temperature.
c_water = speedSoundWater(TEMPERATURE);

% Calculating the focal distance as the centre of the ring array.
focal_length = RING_DIAMETER / 2;

%% Section 2: Determining lens shape and calculating lens radius.

if LENS_SOUND_SPEED > c_water
    lens_radius = focal_length * (1 - (c_water / LENS_SOUND_SPEED));
    fprintf('Lens is convex, and has a lens radius of %.3fmm', ...
                                                       lens_radius * 1000);
elseif LENS_SOUND_SPEED < c_water
    lens_radius = focal_length * ((c_water / LENS_SOUND_SPEED) - 1);
    fprintf('Lens is concave, and has a lens radius of %.3fmm', ...
                                                       lens_radius * 1000);
elseif LENS_SOUND_SPEED == c_water
    disp(['The speed of sound in the lens material is equal to the '...
        'speed of sound in water. Choose a new lens material']);
end
