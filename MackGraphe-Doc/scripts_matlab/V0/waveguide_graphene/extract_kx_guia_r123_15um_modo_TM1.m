%% este programa ayuda na encontrar el valor sobre el cual se van a extraer la parte reale e imaginaria de kx

eps1 = 1;
eps2 = (3.42)^2;
eps3 = (1.98)^2;

% valmax = 1e-3; %valor se saca de grafico de dphi

% valmax = 1e-3; %%EF = 0,0.1

valmax = 0.3e-3;


c = 3e8;
i = sqrt(-1);

d = 15E-6;

Ef = 0.4;
hbar = 6.582119569E-16; %eV.s


vRre = load('Mreflec_real_300n_R4um_np500_EF04.dat');
vRim = load('Mreflec_imag_300n_R4um_np500_EF04.dat');
vf = load('vector_frequencias_500_R4um.dat');
vq = load('vector_kx_500_R4um.dat');


valini = 106;

% % vRre = load('Mreflec_imag_300n_R4um.dat');
% % vRim = load('Mreflec_real_300n_R4um.dat');
% % vf = load('vector_frequencias_100_R4um.dat');
% % vq = load('vector_kx_100_R4um.dat');

[Q F] = meshgrid(vq,vf);

%%%%%%% generando matriz refleccion sistema de nanoribbons
[nf nc] = size(vRre);

ic = 1;
while ic <= nf

      ni = vRre(ic,1);
      nj = vRre(ic,2);
      valre = vRre(ic,3);
      valim = vRim(ic,3);

      rp12(ni,nj) = valre + i*valim;

      ic = ic + 1;
end

clear vRre;
clear vRim;
%clear vf;
%clear vq;

%%%%%%%%%%%%%%%%%%%%%%%%

w = 2*pi*F;

kz2 = sqrt(eps2.*(w/c).^2 - Q.^2);
kz3 = sqrt(eps3.*(w/c).^2 - Q.^2);

rp23 = crp(kz2,kz3,eps2,eps3,w,0);

fase = kz2*d;
%se coloca -r12 porque en fortran se calculo r21
num = -rp12 + rp23.*exp(i*2*fase);
den = 1 - rp12.*rp23.*exp(i*2*fase);

rp13 = num./den;

phi = angle(den);

dphi = diff(phi,1,2)./diff(Q,1,2);

[nf nc] = size(dphi);

ic = valini;

ickx = 1;

while (ic <= nf)

Vwdphi =   dphi(ic,:);
valwf = vf(ic);
Vwq = Q(ic,1:end-1);


%% se eleimina cono de luz guia

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


[dphimax posm] = max(Vwdphi);
[dphimin posmin] = min(Vwdphi);

kxm = Vwq(posm);

valdphi = dphimax;

icimag = posm + 1;

while (valdphi >= 0.5*(dphimax))

    valdphi = Vwdphi(icimag);

    icimag = icimag +1;

end

posim = icimag-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (posim == posm +1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imkx =  0;%Vwq(posim) - kxm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (Vwdphi(posim) < 0)
     posim = posim + 1;
 end

    Dkxsq = (Vwq(posim) - kxm)^2;

    % V1 = dphimax;
    % V2 = Vwdphi(posim);
    V1 = dphimax - dphimin;
    V2 = Vwdphi(posim) - dphimin;

    Cay = V1/V2;

    imkx = sqrt(Dkxsq/(Cay - 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REkxf(ickx,1) = kxm;
REkxf(ickx,2) = valwf;
REkxf(ickx,3) = imkx;

Imkxf(ickx,1) = valwf;
Imkxf(ickx,2) = imkx;
Imkxf(ickx,3) = kxm;

ickx = ickx + 1;
ic = ic + 1;
end


surf(real(Q(:,1:end-1)),F(:,1:end-1),dphi,'FaceAlpha',1,'EdgeColor','none');view(2);shading interp;colorbar;
#save -ascii realkxvsf_15um_TM01_EF04.dat REkxf
#save -ascii imagkxvsf_15um_TM01_EF04.dat Imkxf


% %% se cambio de 0.99 para 0.95
%
% Iz = find(Q(:,1:end-1)>=0.95*sqrt% V1 = dphimax;
% V2 = Vwdphi(posim);(eps2)*2*pi*F(:,1:end-1)/3E8); % cono de luz guia;
%
% dphi(Iz)=0;
%
% q0 = 0;
% f0 = 2e12;
% q1 = 6e5;
% f1 = 10e12;
% penA = (f1-f0)/(q1-q0);
% const = f0;
%
% recta = penA*Q(:,1:end-1) + const;
%
% Izr = find(F(:,1:end-1) > recta);
%
% dphi(Izr)=0;
%
%
% surf(Q(:,1:end-1),F(:,1:end-1),abs(dphi),'FaceAlpha',1,'EdgeColor','none');view(2);shading interp;colorbar;
%
%
% I = find(abs(dphi)>= valmax);
%
%  Msave = [Q(I) F(I)];
%
% save -ascii 'Real_kx_vs_f_modo1_R4m_dg05_EF02.dat' Msave
