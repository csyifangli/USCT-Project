%% Impedance Matching Calculations
% Author:   Morgan Roberts
% Date:     12/10/18
% 
% DESCRIPTION:
%     A script that chooses the component values needed to build a step-up
%     conjugate-matched L-C network for a transducer at a single design 
%     frequency under steady state excitation. The Maximum Power Transfer
%     Theorem is used [1], along with matching networks based on filter
%     structures for high frequency ultrasound transducers [2].
% 
% REFERRENCES:
%     [1] Moura, L., & Darwazeh, I. (2005). Sinusoidal AC electrical 
%     analysis. In Introduction to Linear Circuit Analysis and Modelling: 
%     From DC to RF pp. 54–65). Elsevier.
%     [2] Moon, J. Y., Lee, J., & Chang, J. H. (2016). Electrical 
%     impedance matching networks based on filter structures for high 
%     frequency ultrasound transducers. Sensors and Actuators, A: Physical,
%     251, 225–233. https://doi.org/10.1016/j.sna.2016.10.025
%
%% Section 1: Defining Parameters.

close all
TRANSDUCER_IMPEDANCE = 5.9996 - 5.9771j;   % [ohms]
GENERATOR_IMPEDANCE  = 50 + 0j;            % [ohms]
FREQUENCY            = 1e6;                % [Hz]
GENERATOR_AMPLITUDE  = 50;                 % [V]
GENERATOR_PHASE      = 0;                  % [rad]
w                    = 2 * pi * FREQUENCY; % [rad/s] 

%% Section 2: Generating time series and generator voltage waveforms.

t  = 0:1e-8:2e-6;        % Time data
vg = real(GENERATOR_AMPLITUDE * exp(1j * (w * t - GENERATOR_PHASE)));

%% Section 3: Performance for direct connection (no matching network).

% Use series combination of transducer and generator impedances to
% calculate generator current waveform ig(t).
ig = real((GENERATOR_AMPLITUDE / abs(GENERATOR_IMPEDANCE + ...
    TRANSDUCER_IMPEDANCE)) * exp(1j * (w * t - GENERATOR_PHASE - ...
    angle(GENERATOR_IMPEDANCE + TRANSDUCER_IMPEDANCE))));

% Calculate resistive (R) and reactive (X) components of transducer and 
% generator impedances.
Rt = real(TRANSDUCER_IMPEDANCE); 
Xt = imag(TRANSDUCER_IMPEDANCE); 
Rg = real(GENERATOR_IMPEDANCE); 
Xg = imag(GENERATOR_IMPEDANCE); 

% Calculate average power delivered by generator for direct connection
Pt = (GENERATOR_AMPLITUDE^2 / 2) * Rt / ((Rg + Rt)^2 + (Xg + Xt)^2);
fprintf('Unmatched average power delivered = %.3f watts.\n',Pt);

%% Section 4: Designing the conjugate-matched network

% We model the combination of the signal generator and the L-C matching
% network as a 'Voltage Source'. 
% Using source transformations, we determine the Thevenin equivalent 
% Voltage, vs(t), and Impedance, Zs, for this 'Voltage Source'.

% We can derive values of C and L for which the Thevenin equivalent 
% Impedance is the complex conjugate of the Transducer Impedance.
C = sqrt((Rg - Rt) / (Rt * Rg^2 * w^2));
L = (Rg^2 * w * C - Xt * (Rg^2 * w^2 * C^2 + 1)) / (w + Rg^2 * w^3 * C^2);

% We calculate the Impedance for these components.
ZC = 1 / (1j * w* C); 
ZL = 1j * w * L;

% We calculate the Thevenin equivalent Impedance of the source, with its
% resistive (Rs) and reactive (Xs) components.
Zs = ZL + (Rg * ZC) / (Rg + ZC);
Rs = real(Zs); 
Xs = imag(Zs);

% We calculate the Thevenin equivalent Voltage waveform vs(t), and its
% amplitude, Vs.
vs = real(Vg * abs(ZC / (Rg + ZC)) * exp(1j * (w * t - GENERATOR_PHASE ...
                                                + angle(ZC / (Rg + ZC)))));
Vs = max(vs);

%%
% To check that the power delivered to the transducer has been maximised,
% we will determine the phase difference between the voltage and the
% current, and calculate the maximum average power delivered to the load.


is = real((Vg*abs(ZC/(Rg + ZC))/abs(Zs + TRANSDUCER_IMPEDANCE))*exp(j*(w*t - GENERATOR_PHASE + angle(ZC/(Rg + ZC)) - angle(Zs + TRANSDUCER_IMPEDANCE))));

fig1 = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
hold on
plot(t,vg,'k'); ylabel('Volts [V]'); ylim([-max([vs vg])*1.2 max([vs vg])*1.2]);
[pkvg,tvg] = findpeaks(vg);
stem(t(tvg(1)),pkvg(1),'k-')
 
yyaxis right
plot(t,ig,'k--'); ylabel('Current [A]');

xlabel('Time [s]'); ylim([-max([is ig])*1.2 max([is ig])*1.2]); 
title('Signal Generator v_g(t) and i_g(t) for un-matched circuit');
[pkig,tig] = findpeaks(ig);
stem(t(tig(1)),pkig(1),'k--')
legend('Sig. Gen. Voltage v_g(t)','V_g','Sig. Gen. Current i_g(t)','I_g');
hold off

fig2 = figure(2);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig2,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(t,vs,'k'); ylabel('Volts [V]'); ylim([-max([vs vg])*1.2 max([vs vg])*1.2]);
yyaxis right
hold on
plot(t,is,'k--'); ylim([-max([is ig])*1.2 max([is ig])*1.2]);
ylabel('Current [A]'); legend('Source Voltage v_s(t)','Source Current i_s(t)'); 
xlabel('Time [s]'); title('Source v_s(t) and i_s(t) for conjugate-matched Thevenin equivalent circuit');
Ptmatch = (Vs^2/2)*Rt/((Rs + Rt)^2 + (Xs + Xt)^2); % Average power delivered to transducer

fprintf('Unmatched current phase angle =  %.2f radians.\n',-angle(GENERATOR_IMPEDANCE + TRANSDUCER_IMPEDANCE));

fprintf('Conjugate matched average power delivered = %.3f watts.\n',Ptmatch);

