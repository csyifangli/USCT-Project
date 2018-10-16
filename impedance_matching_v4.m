close all
Zt = 5.9996 - 5.9771j; % Measured impedance of transducer
Zg = 50 + 0j;           % Impedance of Signal Generator
f = 1e6; 
w = 2*pi*f;    % Frequency 
Vg = 50;                % Signal generator amplitude
theta = 1.6;              % Signal generator phase shift
t = 0:1e-8:2e-6;        % Time data
vg = real(Vg * exp(j * (w * t - theta))); % Signal generator waveform

%%
% First, we will analyse the case when the transducer is directly connected to the
% source, determine whether the load voltages and current are in phae, 
% and calculate the maximum average power delivered to the load.

%vt = real(max(vg)*abs(Zt/(Zg + Zt))*exp(j*(w*t + angle(Zt/(Zg + Zt))))); % Transducer voltage vt(t)
ig = real((max(vg)/abs(Zg + Zt))*exp(j*(w*t - theta - angle(Zg + Zt))));         % Transducer current it(t)

Rt = real(Zt); Xt = imag(Zt); % Resistive and Reactive components of transducer impedance
Rg = real(Zg); Xg = imag(Zg); % Resistive and Reactive components of signal generator impedance
Pt = (Vg^2/2)*Rt/((Rg + Rt)^2 + (Xg + Xt)^2); % Average power delivered to transducer

%%
% Next, we introduce a matching network to maximise the power delivered to
% the transducer by bringing the transducer current in phase with the
% voltage.
% We can model the combination of the signal generator and the matching
% network as a 'voltage source'. Using source transformations, we can
% determine the Thevenin equivalent circuit for this 'Voltage Source'.
% We need to derive values of C and L for which Zs = Zt*.

C = sqrt((Rg - Rt)/(Rt*Rg^2*w^2));
L = (Rg^2*w*C - Xt*(Rg^2*w^2*C^2 +  1))/(w + Rg^2*w^3*C^2);

Zc = 1/(j*w*C); % Matching capacitor impedance
Zl = j*w*L;     % Matching inductor impedance

Zs = Zl + (Rg*Zc)/(Rg + Zc); % 'Voltage source' Impedance (Thevenin equivalent)
Rs = real(Zs); Xs = imag(Zs); % Resistive and Reactive components of 'Voltage source' impedance

% The total impedance seen by the Source voltage (Thevenin equivalent) is
% now purely real.

%%
% To check that the power delivered to the transducer has been maximised,
% we will determine the phase difference between the voltage and the
% current, and calculate the maximum average power delivered to the load.

vs = real(Vg*abs(Zc/(Rg + Zc))*exp(j*(w*t - theta + angle(Zc/(Rg + Zc)))));
Vs = max(vs);
is = real((Vg*abs(Zc/(Rg + Zc))/abs(Zs + Zt))*exp(j*(w*t - theta + angle(Zc/(Rg + Zc)) - angle(Zs + Zt))));

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

fprintf('Unmatched current phase angle =  %.2f radians.\n',-angle(Zg + Zt));
fprintf('Unmatched average power delivered = %.3f watts.\n',Pt);
fprintf('Conjugate matched average power delivered = %.3f watts.\n',Ptmatch);

