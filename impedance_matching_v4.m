%% Impedance Matching Calculations
% Author:   Morgan Roberts
% Date:     12/10/18
% 
% DESCRIPTION:
%     A script that chooses the component values needed to build a step-up
%     conjugate-matching network for a transducer at a single design 
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
ANGULAR_FREQUENCY    = 2 * pi * FREQUENCY; % [rad/s] 
GENERATOR_AMPLITUDE  = 50;                 % [V]
GENERATOR_PHASE      = 0;                  % [rad]

%% Section 2: Generating time series and sig. gen. voltage waveforms.

t  = 0:1e-8:2e-6;        % Time data
vg = real(GENERATOR_AMPLITUDE * exp(1j * (ANGULAR_FREQUENCY * t - ...
                                                        GENERATOR_PHASE)));

%% Section 3: Performance for direct connection (no matching network)
% First, we will analyse the case when the transducer is directly connected to the
% source, determine whether the load voltages and current are in phae, 
% and calculate the maximum average power delivered to the load.

%vt = real(max(vg)*abs(Zt/(Zg + Zt))*exp(j*(w*t + angle(Zt/(Zg + Zt))))); % Transducer voltage vt(t)
ig = real((max(vg)/abs(GENERATOR_IMPEDANCE + TRANSDUCER_IMPEDANCE))*exp(j*(ANGULAR_FREQUENCY*t - GENERATOR_PHASE - angle(GENERATOR_IMPEDANCE + TRANSDUCER_IMPEDANCE))));         % Transducer current it(t)

Rt = real(TRANSDUCER_IMPEDANCE); Xt = imag(TRANSDUCER_IMPEDANCE); % Resistive and Reactive components of transducer impedance
Rg = real(GENERATOR_IMPEDANCE); Xg = imag(GENERATOR_IMPEDANCE); % Resistive and Reactive components of signal generator impedance
Pt = (Vg^2/2)*Rt/((Rg + Rt)^2 + (Xg + Xt)^2); % Average power delivered to transducer

%%
% Next, we introduce a matching network to maximise the power delivered to
% the transducer by bringing the transducer current in phase with the
% voltage.
% We can model the combination of the signal generator and the matching
% network as a 'voltage source'. Using source transformations, we can
% determine the Thevenin equivalent circuit for this 'Voltage Source'.
% We need to derive values of C and L for which Zs = Zt*.

C = sqrt((Rg - Rt)/(Rt*Rg^2*ANGULAR_FREQUENCY^2));
L = (Rg^2*ANGULAR_FREQUENCY*C - Xt*(Rg^2*ANGULAR_FREQUENCY^2*C^2 +  1))/(ANGULAR_FREQUENCY + Rg^2*ANGULAR_FREQUENCY^3*C^2);

Zc = 1/(j*ANGULAR_FREQUENCY*C); % Matching capacitor impedance
Zl = j*ANGULAR_FREQUENCY*L;     % Matching inductor impedance

Zs = Zl + (Rg*Zc)/(Rg + Zc); % 'Voltage source' Impedance (Thevenin equivalent)
Rs = real(Zs); Xs = imag(Zs); % Resistive and Reactive components of 'Voltage source' impedance

% The total impedance seen by the Source voltage (Thevenin equivalent) is
% now purely real.

%%
% To check that the power delivered to the transducer has been maximised,
% we will determine the phase difference between the voltage and the
% current, and calculate the maximum average power delivered to the load.

vs = real(Vg*abs(Zc/(Rg + Zc))*exp(j*(ANGULAR_FREQUENCY*t - GENERATOR_PHASE + angle(Zc/(Rg + Zc)))));
Vs = max(vs);
is = real((Vg*abs(Zc/(Rg + Zc))/abs(Zs + TRANSDUCER_IMPEDANCE))*exp(j*(ANGULAR_FREQUENCY*t - GENERATOR_PHASE + angle(Zc/(Rg + Zc)) - angle(Zs + TRANSDUCER_IMPEDANCE))));

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
fprintf('Unmatched average power delivered = %.3f watts.\n',Pt);
fprintf('Conjugate matched average power delivered = %.3f watts.\n',Ptmatch);

