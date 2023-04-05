function  rp = crp(kz1,kz2,eps1,eps2,w,sigma)

rp = -((eps2.*kz1 - eps1.*kz2 +((kz1.*kz2).*(sigma./(w.*8.8419e-12))))./(eps2.*kz1+ eps1.*kz2+((kz1.*kz2).*(sigma./(w.*8.8419e-12)))));

end
