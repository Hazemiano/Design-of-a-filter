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
pole1=0.6 ;%% or -0.6

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


