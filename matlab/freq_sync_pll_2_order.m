symbol_rate = 97.65625e3/2;
oversample = 1; sig_level = 1;

fs=symbol_rate*oversample; 
Ts=1/fs;

% Loop filter parameters
Bw=135;              % BW of PLL / carrier synchronization
zeta=1/sqrt(2);      % 1/sqrt(2)
kp_c=1;              % Proportional constant
kp_c=kp_c*sig_level; % Proportional constant - account signal level - carrier synchronization
k0_c=1;              % k0

alfa_p1=2*zeta*((2*Bw*Ts)/(zeta+(1/(4*zeta))));
beta_p2=((2*Bw*Ts)/(zeta+(1/(4*zeta))))^2;

alfa=alfa_p1*(1/kp_c)*(1/k0_c);
beta=beta_p2*(1/kp_c)*(1/k0_c);

alfa
beta