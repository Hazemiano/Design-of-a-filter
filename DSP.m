
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

