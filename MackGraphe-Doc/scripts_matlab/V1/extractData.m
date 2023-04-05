function extractData(eps1,eps2,eps3)
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

dphi = diff(arg(den),1,2)/pq;

[nf nc] = size(dphi);

% % selecionar dependendo da posicao no vetro de frequancias vf onde começa o modo 1. f = 2.4600e+12
valini = 50;
ic = valini;
ickx = 1;
imkx = 0;


while (ic <= nf)

  Vwdphi = dphi(ic,:);
  valwf = vf(ic);
  Vwq = Q(ic,1:end-1);

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

  %% se elimina valores negativos para hacer fit

  Inz = find(Vwdphi > 0);
  NVwdphi = Vwdphi(Inz);
  NVwq = Vwq(Inz);

  [dphimax posm] = max(NVwdphi);

  kxm = NVwq(posm);

  [ff cf] = size(Inz);

  alpha_loss = -1;
  if (cf > 1) % si no entra se grava igual el imkx

      [yprime1 params1 resnorm1 residual1 jacobain1] = lorentzfit(NVwq,NVwdphi);

      P1 = params1(1);
      P2 = params1(2);
      P3 = params1(3);
      C = params1(4);

      cay = abs(0.5*(P1/P3 - C));
      imkx = sqrt(abs(P1/cay - P3));
      alpha_loss = 20*(imkx/log(10));

  end


  [x y] = size(valwf);
  if(y == 0)
    Imkxf(ickx,1) = 0;
  else
    Imkxf(ickx,1) = valwf;
  endif

  Imkxf(ickx,2) = imkx;
  Imkxf(ickx,3) = alpha_loss;

  clear Inz;
  clear NVwq;
  clear NVwdphi;

  ickx = ickx + 1;
  ic = ic + 1;
end

######################################################
# Extração de Dados: Base de Imagens - Cone de Luz
surf(real(Q(:,1:end-1)),F(:,1:end-1), dphi,'FaceAlpha',1,'EdgeColor','none');view(2);

imgfilename = 'imgs/cone_luz_sem_grafeno_n2_(';
n2_real = real(eps2);
n2_img = imag(eps2);
n2_real_str = num2str(n2_real);
n2_img_str = num2str(n2_img);
n2_full_str = [n2_real_str '_+_'];
n2_full_str = [n2_full_str n2_img_str];
n2_full_str_ext = [n2_full_str 'i).jpg'];
imgfilename = [imgfilename n2_full_str_ext]

x0=100;
y0=100;
width=300;
height=300;
set(gca,'XTick',[], 'YTick', []);
set(gcf,'position',[x0,y0,width,height]);
saveas(gcf,imgfilename);

############################################################
# Extração de Dados: Base de Dados - Perdas por frequencia

csvfilename = 'losses/cone_luz_sem_grafeno_n2_(';
n2_full_str_ext = [n2_full_str ').csv'];
csvfilename = [csvfilename n2_full_str_ext]

csvwrite(csvfilename, Imkxf);

endfunction
