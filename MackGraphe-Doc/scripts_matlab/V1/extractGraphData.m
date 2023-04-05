function extractGraphData(n2, n3)
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

  kz2 = sqrt(n2.*(w/c).^2 - Q.^2);
  kz3 = sqrt(n3.*(w/c).^2 - Q.^2);

  rp23 = crp(kz2,kz3,n2,n3,w,0);

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

  alpha_loss = -99.99;
  while (ic <= nf)
    Vwdphi =   dphi(ic,:);
    valwf = vf(ic);
    Vwq = Q(ic,1:end-1);

    %% se eleimina cono de luz guia

    vsup = 0.98*sqrt(n2)*2*pi*valwf/c;
    vinf = sqrt(n3)*2*pi*valwf/c;

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
        imkx = sqrt(abs(Dkxsq/(Cay - 1)));
        alpha_loss = 20*(imkx/log(10));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [x y] = size(valwf);
    if(y == 0)
      Imkxf(ickx,1) = 0;
    else
      Imkxf(ickx,1) = valwf;
    endif
    Imkxf(ickx,2) = imkx;
    #Imkxf(ickx,3) = kxm;
    Imkxf(ickx,4) = alpha_loss;

    ickx = ickx + 1;
    ic = ic + 1;
  end
  ######################################################
  # Extração de Dados: Base de Imagens - Cone de Luz
  surf(real(Q(:,1:end-1)),F(:,1:end-1), dphi,'FaceAlpha',1,'EdgeColor','none');view(2);

  imgfilename = 'imgs_grafeno/cone_luz_grafeno_n2_(';
  n2_real = real(n2);
  n2_img = imag(n2);
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
  csvfilename = 'losses_grafeno/cone_luz_grafeno_n2_(';
  n2_full_str_ext = [n2_full_str ').csv'];
  csvfilename = [csvfilename n2_full_str_ext]

  csvwrite(csvfilename, Imkxf);

endfunction
