i=sqrt(-1);

n1 = 1;
n2 = 1.77;
n3 = 1.45;

lamb= 1e-6;
d = 1e-6;

px = 0.0005;
x = (0:px:1.77);

r12 = func_rp_x(n1,n2,x);
r23 = func_rp_x(n2,n3,x);
phi = (2*pi/lamb)*sqrt(n2.^2-x.^2)*d;

num = (r12 + r23.*exp(i*2*phi));
den = (1 + r12.*r23.*exp(i*2*phi));

r123 = num./den;

%plot(x,arg(r123));
%plot(x,arg(den)/pi);

dphi = diff(arg(den))/px;

vxd = [x(:,1:end-1)' dphi'];

[nvxd si] = sortrows(vxd,2);

modo1 = nvxd(end,1)
modo2 = nvxd(end-1,1)



plot(x(:,1:end-1),dphi);xlim([n3 n2]);

xlabel ("x");
ylabel ("fase");

%plot(x,r123.*conj(r123));
%xlim([0 90]);
%xlabel ("Ângulo de indeciência");
%ylabel ("R");
%legend ("TM","TE",'location','north');
%title("Coeficiente de refletividade:Ar-vidro")
