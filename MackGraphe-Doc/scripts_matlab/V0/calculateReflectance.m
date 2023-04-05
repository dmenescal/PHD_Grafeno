function [kx, ref, imref, dphi] = calculateImgKx(eps1,eps2,eps3)
  c = 3e8;
  d = 5E-6;
  Ef = 0.4;
  hbar = 6.582119569E-16; %eV.s
  fmin = 0.01;
  fmax = 10;
  vf = (fmin:0.05:fmax)*1e12; % Intervalos de frequencia
  kmax = sqrt(eps2)*2*pi*vf(end)/c;
  pq = 0.005;
  vkx = (0.001:pq:1.5)*kmax; % Intervalos de Kx
  [kx f] = meshgrid(vkx,vf);
  w = 2*pi*f;
  kz1 = sqrt( eps1.*(w/c).^2 - kx.^2);
  kz2 = sqrt( eps2.*(w/c).^2 - kx.^2);
  kz3 = sqrt( eps3.*(w/c).^2 - kx.^2);
  sigma = 0;
  rp12 = crp(kz1,kz2,eps1,eps2,w,sigma);
  rp23 = crp(kz2,kz3,eps2,eps3,w,0);
  fase = kz2*d;
  num = (rp12 + rp23.*exp(i*2*fase));
  den = (1 + rp12.*rp23.*exp(i*2*fase));

  ref = num./den;
  imref = abs(imag(ref));
  dphi = diff(arg(den),1,2)/pq;
end
