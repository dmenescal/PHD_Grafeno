function rp = func_rp_x(n1,n2,x)

rp = (n2.^2.*sqrt(n1.^2 - x.^2) - n1.^2.*sqrt(n2.^2 - x.^2))./(n2.^2.*sqrt(n1.^2 - x.^2) + n1.^2.*sqrt(n2.^2 - x.^2));

end
