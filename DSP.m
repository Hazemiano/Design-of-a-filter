
function system_analysis(num, den, wp)
w = linspace(-pi, pi,1024);

figure();
zplane(num,den);
grid on;
title("Pole-zero in Z domain");

H=freqz(num,den, w);

mag_dB =20*log10(abs(H));

figure();
plot(w,mag_dB);
grid on;
xlim([-pi,pi]);
xlabel("\omega (radian)");
ylabel("Magnitude");
title('Magnitude Response(-pi->pi)');

figure();
plot(w,mag_dB);
grid on;
xlim([-wp,wp]);
xlabel("\omega (radians)");
ylabel("Magnitude");
title(' Magnitude Response (-wp->wp)');

phi = angle(H);

figure();
plot(w,phi);
grid on;
xlim([-pi, pi]);
xlabel("\omega (radians)");
ylabel("Phase (radians)");
title('Phase Response');

dphi = diff(phi);
dw = diff(w);

gd = -dphi ./ dw;

figure();
plot(w(1:end-1), gd);
grid on;
xlim([-pi,pi]);
xlabel('\omega (radians)');
ylabel('Group Delay');
title('Group Delay');

figure();
[h,t] = impz(num, den,50);
stem(t,h);
grid on;
xlabel('n');
ylabel('h[n]');
title('Impulse Response');

end

%%syms r

pole1=0.6 ;%% or -0.6

wc = 0.25*pi;
width = 0.1*pi;
ws = wc + (width/2);
wp = wc - (width/2);

%%zeros
zero1 =exp(j*ws);
zero2 =exp(-j*ws);
zeros =poly([zero1 zero2]);
W=linspace(0,wp,1024);
Ripples =1;  %%dB
Counter =0;
for r = 0 : 0.001 : 1
    Counter = Counter + 1;
    pole2 = r*exp(j*wp);
    pole3 = r*exp(-j*wp);
    poles = poly([pole1 pole2 pole3]);
    
    H =freqz(zeros,poles,W);
    Mag_dB=20*log10(abs(H));
    mag(Counter)=r;
    Storage(Counter)=range(Mag_dB);
   
end
[minR,i]=min(Storage);
optimal_r=mag(i); 
fprintf('minimum ripples: %.3f dB\n',minR);
fprintf('optimal r: %.3f\n',optimal_r);



%% Point 6
wp2 = wp / 2;
ws2 = (pi + ws) / 2;

pole2 = optimal_r * exp(j * wp);
pole3 = optimal_r * exp(-j * wp);

pole4 = optimal_r * exp(j * wp2);
pole5 = optimal_r * exp(-j * wp2);
poles = poly([pole1 pole2 pole3 pole4 pole5]);

zero3 = exp(j * ws2);
zero4 = exp(-j * ws2);
zeros = poly([zero1 zero2 zero3 zero4]);

W = linspace(0,wp,1024);
H = freqz(zeros,poles,W);
Mag_dB = 20*log10(abs(H));
Ripples = range(Mag_dB);
fprintf("Passband ripples before adjusting magnitude of p1: %.3f dB\n", Ripples);

Counter = 0;
for r = -1 : 0.001 : 1
    Counter = Counter + 1;
    pole1 = r;
    poles = poly([pole1 pole2 pole3 pole4 pole5]);

    H = freqz(zeros,poles,W);
    Mag_dB = 20 * log10(abs(H));
    mag2(Counter) = r;
    Storage2(Counter) = range(Mag_dB);

end
[minR,i]=min(Storage2);
optimal_r=mag2(i); 
fprintf('minimum ripples: %.3f dB\n',minR);
fprintf('optimal r: %.3f\n',optimal_r);
