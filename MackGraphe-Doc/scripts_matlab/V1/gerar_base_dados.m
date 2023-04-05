rng('default') % For reproducibility
i = sqrt(-1);

eps2_real = (3.42)^2;
delta = 0.1;

eps2_range = [0:delta:2*eps2_real];%parte imaginaria para o eps2

v_kx = zeros(length(eps2_range));
v_alpha = zeros(length(eps2_range));
for n = 1:length(eps2_range)
  [imkx dphi alpha_loss] = calculateImgKx(1,eps2_real + eps2_range(n), (1.98)^2);
  v_kx(n) = imkx;
  v_alpha(n) = alpha_loss;
  surf(log10(abs(dphi)),'FaceAlpha',1,'EdgeColor','none');
  view(2);
  x0=100;
  y0=100;
  width=300;
  height=300;
  set(gca,'XTick',[], 'YTick', []);
  set(gcf,'position',[x0,y0,width,height]);
  file_name = 'imagens/waveguide_complex_';
  file_index = num2str(n);
  file_index = [file_index '_'];
  n2_complex_i = num2str(eps2_range(n));
  file_name = [file_name file_index];
  file_name = [file_name 'n2_']
  file_name = [file_name n2_complex_i];
  file_name = [file_name 'i_complex(kx)_'];
  kx_complex_i = num2str(imkx);
  file_name = [file_name kx_complex_i];
  loss = num2str(alpha_loss);
  file_name = [file_name 'i_loss(alpha)_'];
  file_name = [file_name loss];
  file_name = [file_name '.png'];
  file_name
  saveas(gcf,file_name)
end




