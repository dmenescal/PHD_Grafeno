%% este programa ayuda na encontrar el valor sobre el cual se van a extraer la parte reale e imaginaria de kx



eps1 = 1;
eps2 = (3.42)^2 + 0i;
eps3 = (1.98)^2;

fmin = 0.01;
fmax = 10;


c = 3e8;

d = 15E-6;


vf = (fmin:0.05:fmax)*1e12;

pq = 0.005;

kmax = real(sqrt(eps2)*2*pi*vf(end)/c);

vq = (0.001:pq:1.5)*kmax;

% % selecionar dependendo da posicao no vetro de frequancias vf onde começa o modo 1. f = 2.4600e+12
valini = 50;

[Q F] = meshgrid(vq,vf);

w = 2*pi*F;


kz1 = sqrt(eps1.*(w/c).^2 - Q.^2);
kz2 = sqrt(eps2.*(w/c).^2 - Q.^2);
kz3 = sqrt(eps3.*(w/c).^2 - Q.^2);

rp12 = crp(kz1,kz2,eps1,eps2,w,0);
rp23 = crp(kz2,kz3,eps2,eps3,w,0);

fase = kz2*d;

num = rp12 + rp23.*exp(1i*2*fase);
den = 1 + rp12.*rp23.*exp(1i*2*fase);

rp13 = num./den;

phi = angle(den);

dphi = diff(arg(den),1,2)/pq;

[pico, ic] = max(max(dphi(1:50, :)));

ic

surf(real(Q(:,1:ic)),F(:,1:ic), dphi(:,1:ic),'FaceAlpha',1,'EdgeColor','none'); colorbar;view(2);

Vwdphi = dphi(ic,:);
valwf = vf(ic);
Vwq = Q(ic,:);

%% se elimina cono de luz guia

vsup = 0.98*sqrt(real(eps2))*2*pi*valwf/c;
vinf = sqrt(real(eps3))*2*pi*valwf/c;

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

imkx
#surf(real(Q(:,1:end-1)),F(:,1:end-1), new_dphi,'FaceAlpha',1,'EdgeColor','none'); colorbar;view(2);
#save -ascii realkxvsf_15um_TM01_Lorentz.dat REkxf
#save -ascii imagkxvsf_15um_TM01_Lorentz.dat Imkxf


