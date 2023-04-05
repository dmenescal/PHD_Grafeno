%% DBA 01022021 en este programa se grafica la parte imaginaria de la reflexion para ver la disperiosn de plasmons
%% para un sistema de 3 medios


eps1 = 1;
eps2 = (3.42)^2 + 0*i;
eps3 = (1.98)^2;

c = 3e8;

d = 15E-6;

Ef = 0.4;
hbar = 6.582119569E-16; %eV.s

fmin = 0.01;
fmax = 10;

vf = (fmin:0.05:fmax)*1e12;

t = sqrt(eps2);

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

phi = angle(den);

dphi = diff(arg(den),1,2)/pq;
#dphi = diff(phi,1,2)./diff(real(kx),1,2);

#surf(kx,f,log(imrp13));view(2);shading interp;
#surf(kx,f,(abs(dphi)));view(2);shading interp;

#surf(real(kx(:,1:end-1)),f(:,1:end-1), dphi,'FaceAlpha',1,'EdgeColor','none'); colorbar;view(2);
#plot(real(kx(:,1:end-1)), phi);


# Obtém os valores máximos, para cada frequencia (indice e valor)
#[pico, ic] = max(max(dphi));

ic = 106;

# Obtém apenas os valores máximos do dphi (relativos aos possíveis modos guiados)
# Nao acontece pois cada coluna tem um valo maximo, logo esse filtro nao tem mudanca
Vwdphi = dphi(ic, :);

# Obtém apenas as frequências que geraram picos de dphi
valwf = vf(ic);

# Obtém apenas os kx que geraram picos de dphi
Vwq = vkx(ic);

vsup = real(0.98*sqrt(eps2)*2*pi*valwf/c);
vinf = sqrt(eps3)*2*pi*valwf/c;

M_vwq = max(Vwq);

Iinf = find(Vwq <= vinf);
Isup = find(Vwq >= vsup);

Vwdphi(Iinf) = 0;
Vwdphi(Isup) = 0;

%% recta para eliminar los otros modos para alta frequencia

q0 = 0;
f0 = 2e12;
q1 = 6e5;
f1 = 9e12;
penA = (f1-f0)/(q1-q0);
const = f0;

recta = (valwf - const)/penA;
Izr = find(Vwq(1,:) <= recta);
Izr = [Izr];

Zeros = find(Vwq(1,:) > recta);
Zeros = [Zeros];
New_dphi = dphi;
New_dphi(Zeros) = 1000;

hold on;
surf(real(kx(:,1:end-1)),f(:,1:end-1), New_dphi,'FaceAlpha',1,'EdgeColor','none'); colorbar;view(2);
reta = (f(:,1:end-1) - const)/penA;
plot(reta, f(:,1:end-1),'w-');
hold off;


####tam = length(Izr);
####
####izr = zeros(tam);
####
####for j = 1:tam
####  if (Vwq(j) < recta)
####    izr(j) =  j;
####  endif
####end
####
####result = izr(:,1);
##Vwdphi(Izr) = 0;
##
##hold on;
##
##[dphimax posm] = max(Vwdphi);
##[dphimin posmin] = min(Vwdphi);
##
##zero = 0;
##Inz = find(Vwdphi > zero);
##NVwdphi = Vwdphi(Inz);
##NVwq = Vwq(Inz);
##
##[yprime1 params1 resnorm1 residual1 jacobain1] = lorentzfit(NVwq,NVwdphi);
##
##plot(NVwq/(2*pi*valwf/c),NVwdphi,'*k',NVwq/(2*pi*valwf/c),yprime1);xlim([0 4])
##
##P1 = params1(1);
##P2 = params1(2);
##P3 = params1(3);
##C = params1(4);
##
##kxm = P2;
##
##cay = 0.5*(P1/P3 - C);
##
##imkx = sqrt(P1/cay - P3);
##Izr;
##
##%surf(real(kx(:,1:end-1)),abs(rp13(:,1:end-1)),f(:,1:end-1),'FaceAlpha',1,'EdgeColor','none'); colorbar;view(2);
##
##%surf(kx,f,log(abs(rp13)));shading interp;
##%%% truco para plotar cuando tiene 0's abs'newM = imrp13 + exp(-50);
##
##%x = kx(:,1:end-1);
##%y = dphi;
##
##%[yprime1 params1 resnorm1 residual1] = lorentzfit(x,y);
##
##%[pico, idx] = max(max(yprime1)); %Obtem o valor onde yprime atinge um pico.
##
##%kx_pico_real = max(x(:,idx)); %Obtem a altura do pico indicada pelo yprime.
##
##%kx_img;
##
##%surf(x, y, yprime1);colorbar;view(2);
##%figure; %plot(x,y,'b.','LineWidth',2);
##%hold on; plot(x,yprime1,'r-','LineWidth',2);
