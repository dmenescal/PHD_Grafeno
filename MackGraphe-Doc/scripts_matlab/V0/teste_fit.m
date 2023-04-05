x = -16:0.1:35;
y = 19.4./((x - 7).^2 + 15.8) + randn(size(x))./10;
[yprime1 params1 resnorm1 residual1] = lorentzfit(x,y);
figure; plot(x,y,'b.','LineWidth',2)
hold on; plot(x,yprime1,'r-','LineWidth',2)
