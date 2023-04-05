function [imkx, dphi, alpha_loss] = calculateImgKx(eps1,eps2,eps3)
c = 3e8;

d = 15E-6;

Ef = 0.4;
hbar = 6.582119569E-16; %eV.s

fmin = 0.01;
fmax = 10;

vf = (fmin:0.05:fmax)*1e12;

kmax = real(sqrt(eps2)*2*pi*vf(end)/c);

pq = 0.005;

vkx = (0.001:pq:1.5)*kmax;

[kx f] = meshgrid(vkx,vf);

w = 2*pi*f;

kz1 = sqrt( eps1.*(w/c).^2 - kx.^2);
kz2 = sqrt( eps2.*(w/c).^2 - kx.^2);
kz3 = sqrt( eps3.*(w/c).^2 - kx.^2);

%sigma = condgrafPeresNuno(Ef,w);

sigma = 0;% condgrafDrude(Ef,w);

rp12 = crp(kz1,kz2,eps1,eps2,w,sigma);

%rp12 = crp(kz1,kz2,eps1,eps2,w,0);
rp23 = crp(kz2,kz3,eps2,eps3,w,0);

fase = kz2*d;

num = (rp12 + rp23.*exp(i*2*fase));
den = (1 + rp12.*rp23.*exp(i*2*fase));

rp13 = num./den;

dphi = diff(arg(den),1,2)/pq;

[pico, ic] = max(dphi);

Vwdphi = dphi(ic,:);
valwf = vf(ic);
Vwq = kx(ic,1:end-1);

vsup = 0.98*sqrt(eps2)*2*pi*valwf/c;
vinf = sqrt(eps3)*2*pi*valwf/c;

Iinf = find(Vwq <= vinf);
Isup = find(Vwq >= vsup);

Vwdphi(Iinf) = 0;
Vwdphi(Isup) = 0;

%% recta para eliminar los otros modos para alta frequencia

q0 = 0;
f0 = 2e12;
q1 = 6e5;
f1 = 10e12;
penA = (f1-f0)/(q1-q0);
const = f0;

recta = (valwf - const)/penA;
Izr = find(Vwq < recta);

Vwdphi(Izr) = 0;

hold on;

[dphimax posm] = max(Vwdphi);
[dphimin posmin] = min(Vwdphi);

Inz = find(Vwdphi > 0);
NVwdphi = Vwdphi(Inz);
NVwq = Vwq(Inz);

[yprime1 params1 resnorm1 residual1 jacobain1] = lorentzfit(NVwq,NVwdphi);

%plot(NVwq/(2*pi*valwf/c),NVwdphi,'*k',NVwq/(2*pi*valwf/c),yprime1);xlim([0 4])

P1 = params1(1);
P2 = params1(2);
P3 = params1(3);
C = params1(4);

kxm = P2;

cay = 0.5*(P1/P3 - C);
imkx = sqrt(P1/cay - P3);
alpha_loss = 20*(imkx/log(10));

endfunction
