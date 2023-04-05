    function rp = func_rs_x(n1,n2,x)

rp = (sqrt(n1.^2 - x.^2) - sqrt(n2.^2 - x.^2))./(sqrt(n1.^2 - x.^2) + sqrt(n2.^2 - x.^2));

end
