function system_analysis(num, den, wp)
    w = linspace(-pi, pi, 1024);
    figure('Name', 'Filter Analysis Dashboard', 'Position', [100, 100, 1200, 800]);
    H = freqz(num, den, w);
    mag_dB = 20*log10(abs(H));
    % 1. Pole-Zero Plot
    subplot(2, 3, 1);
    zplane(num, den);
    grid on;
    title("Pole-Zero Map");
    % 2. Full Magnitude Response (-pi to pi)
    subplot(2, 3, 2);
    plot(w, mag_dB);
    grid on;
    xlim([-pi, pi]);
    xlabel("\omega (radians)");
    ylabel("Magnitude (dB)");
    title('Magnitude Response (-\pi \rightarrow \pi)');
    % 3. Passband Magnitude Response (-wp to wp) to check ripples
    subplot(2, 3, 3);
    plot(w, mag_dB);
    grid on;
    xlim([-wp, wp]);
    xlabel("\omega (radians)");
    ylabel("Magnitude (dB)");
    title('Magnitude Response (-\omega_p \rightarrow \omega_p)');
    % 4. Phase Response
    subplot(2, 3, 4);
    phi = angle(H);
    plot(w, phi);
    grid on;
    xlim([-pi, pi]);
    xlabel("\omega (radians)");
    ylabel("Phase (radians)");
    title('Phase Response');
    % 5. Group Delay
    subplot(2, 3, 5);
    gd = grpdelay(num, den, w);
    plot(w, gd);
    grid on;
    xlim([-pi, pi]);
    xlabel('\omega (radians)');
    ylabel('Group Delay (samples)');
    title('Group Delay');
    % 6. Impulse Response
    subplot(2, 3, 6);
    [h, t] = impz(num, den, 50);
    stem(t, h);
    grid on;
    xlabel('n');
    ylabel('h[n]');
    title('Impulse Response');
end
%%%%%%%%%%%

figure
%% plotting zeros and poles of transfer function
A= -0.6;
a = 1;
b = [1 -A];
subplot(3,2,1)
zplane(a,b);
title('Z-plane with Poles and Zeros');
%% Magnitude response in dB in freq range from -pi to pi
w = -pi:pi/100:pi;
mag = 1./sqrt ((1-A*cos(w)).^2+(A*sin(w)).^2);
mag_dB = 20*log10(mag);
subplot(3,2,2)
grid on;
plot (w,mag_dB);
xlabel ('w');
ylabel ('mag_ dB');
title('Magnitude response in dB');
%% Magnitude response in dB in freq range from wp to pi & evaluating the ripples
w = 2.91:pi/1000:pi;
mag = 1./sqrt ((1-A*cos(w)).^2+(A*sin(w)).^2);
mag_dB = 20*log10(mag);
subplot(3,2,3)
grid on;
plot (w,mag_dB);
xlabel ('w');
ylabel ('mag_ dB');
title('Magnitude response in dB in freq range from wp to pi');
%% Phase response
w = -pi:pi/100:pi;
real_part = (1 - A*cos(w));
imag_part = (A*sin(w));
phase = -atan2(imag_part, real_part);
subplot(3,2,4)
plot(w, phase);
grid on;
xlabel('\omega');
ylabel('Phase (radians)');
title('Phase Response');
% numerical derivative
dphi = diff(phase);
dw = diff(w);
gd = -dphi ./ dw;
% plotting
subplot(3,2,5)
plot(w(1:end-1), gd);
grid on;
xlabel('\omega');
ylabel('Group Delay');
title('Group Delay');
%% calculation of inverse z-transform
syms z n
H_z = 1/(1 - A*z^(-1));
h_n = iztrans(H_z, z, n);
%% Impulse response plot
n_plot = 0:36;
h_plot = double(subs(h_n, n, n_plot));
subplot(3,2,6)
stem(n_plot,h_plot,'filled');
grid on;
xlabel('n');
ylabel('h[n]');
title('Impulse Response');




%%%%%%
pole1=0.6 ;%%

wc = 0.25 * pi;          % Cutoff frequency
width = 0.1 * pi;        % Transition band width
ws = wc + (width / 2);   % Stopband edge
wp = wc - (width / 2);   % Passband edge

%%% zeros
zero1 = 0.3*exp(j*ws); 
zero2 = 0.3*exp(-j*ws);
zeros_poly = poly([zero1 zero2]);
W = linspace(-(wp), wp, 1024);
Ripples = 1;  %%dB
Counter = 0;

for r = 0 : 0.001 : 1
        Counter = Counter + 1;
        pole2 = r*exp(j*wp);
        pole3 = r*exp(-j*wp);
        poles = poly([pole1 pole2 pole3]);
    
        H = freqz(zeros_poly,poles,W);
        Mag_dB = 20*log10(abs(H));
        mag(Counter) = r;
        Storage(Counter) = max(Mag_dB) - min(Mag_dB);
   
end   
[minR,i] = min(Storage);
optimal_r = mag(i); 
fprintf('minimum ripples: %.3f dB\n', minR);
fprintf('optimal r: %.3f\n', optimal_r);

% plotting
num = poly([zero1 zero2]);
den = poly([pole1 optimal_r*exp(1j*wp) optimal_r*exp(-1j*wp)]);
system_analysis(num, den, wp);

%%% Section 6: Fifth Order Lowpass Filter
wp2 = wp/2;
ws2 = (pi + ws) / 2;

pole2 = optimal_r * exp(j * wp);
pole3 = optimal_r * exp(-j * wp);
pole4 = optimal_r * exp(j * wp2);
pole5 = optimal_r * exp(-j * wp2);
poles = poly([pole1 pole2 pole3 pole4 pole5]);

zero3 = exp(j * ws2);
zero4 = exp(-j * ws2);
zeros_poly = poly([zero1 zero2 zero3 zero4]);

W = linspace(-wp2, wp2, 1024);
H = freqz(zeros_poly, poles, W);
Mag_dB = 20*log10(abs(H));
Ripples = range(Mag_dB);
fprintf("Passband ripples before adjusting magnitude of p1: %.3f dB\n", Ripples);

Counter = 0;
for r = 0 : 0.001 : 1
   Counter = Counter + 1;
   pole1 = r;
   poles = poly([pole1 pole2 pole3 pole4 pole5]);
   H = freqz(zeros_poly, poles, W);
   Mag_dB = 20*log10(abs(H));
   mag2(Counter) = r;
   Storage2(Counter) = max(Mag_dB) - min(Mag_dB);
end
   
[minR,i] = min(Storage2);
optimal_r2 = mag2(i); 
fprintf('minimum ripples: %.3f dB\n', minR);
fprintf('optimal r2: %.3f\n', optimal_r2);

% plotting
num = poly([zero1 zero2 zero3 zero4]);
den = poly([optimal_r2 optimal_r*exp(1j*wp) optimal_r*exp(-1j*wp) optimal_r*exp(1j * wp2) optimal_r*exp(-1j * wp2)]);
system_analysis(num, den, wp2);


%HPF
[z_lpf, p_lpf, k_lpf] = tf2zpk(num, den);

z_hpf = z_lpf * exp(1j * pi); %rotating by π
p_hpf = p_lpf * exp(1j * pi);

[num_hpf, den_hpf] = zp2tf(z_hpf, p_hpf, k_lpf);
num_hpf = real(num_hpf);
den_hpf = real(den_hpf);

%HPF Ripples
W_hpf = linspace(pi - wp2, pi+wp2, 1024); 
H_hpf = freqz(num_hpf, den_hpf, W_hpf);
Mag_dB_hpf = 20*log10(abs(H_hpf));
Ripples_hpf = max(Mag_dB_hpf) - min(Mag_dB_hpf);  
fprintf('\n--- HPF Analysis ---\n');
fprintf('Passband ripples for HPF: %.3f dB\n', Ripples_hpf);

system_analysis(num_hpf, den_hpf, pi - wp);


% BPF centered at pi/2
z_bpf = [z_lpf * exp( 1j * pi/2); %we rotate by both +pi/2 and -pi/2 to restore symmetry (real coefficients)
         z_lpf * exp(-1j * pi/2)];

p_bpf = [p_lpf * exp( 1j * pi/2);
         p_lpf * exp(-1j * pi/2)];

k_bpf = k_lpf^2;

[num_bpf, den_bpf] = zp2tf(z_bpf, p_bpf, k_bpf);
num_bpf = real(num_bpf);
den_bpf = real(den_bpf);

%BPF Ripples
W_bpf = linspace(pi/2 - wp2, pi/2 +wp2, 1024); 
H_bpf = freqz(num_bpf, den_bpf, W_bpf);
Mag_dB_bpf = 20*log10(abs(H_bpf));
Ripples_bpf = max(Mag_dB_bpf) - min(Mag_dB_bpf); 
fprintf('Passband ripples for BPF: %.3f dB\n', Ripples_bpf);

system_analysis(num_bpf, den_bpf, wp);


%%% Section 8: Comb Filter%% Section 8: Comb Filter
L = 8; 
num_5th = poly([zero1 zero2 zero3 zero4]);
den_5th = poly([optimal_r2 pole2 pole3 pole4 pole5]);

len_num = (length(num_5th) - 1) * L + 1;
len_den = (length(den_5th) - 1) * L + 1;

num_comb = zeros(1, len_num);
den_comb = zeros(1, len_den);

num_comb(1:L:end) = num_5th;
den_comb(1:L:end) = den_5th;

system_analysis(num_comb, den_comb, wp);

W_comb = linspace(-wp2/L, wp2/L, 1024); % The passband is compressed by L
H_comb = freqz(num_comb, den_comb, W_comb);
Mag_dB_comb = 20*log10(abs(H_comb));
Ripples_comb = range(Mag_dB_comb);

fprintf('Comb Filter passband ripples: %.3f dB\n', Ripples_comb);