fprintf('\n')
fprintf('**************************************************************\n')
fprintf('Begin imaginary axis calculation.\n')

wn = pi/beta;
numwi = round(Ncut*nuscale/wn);
% numwi = 400 

fprintf('  T = %g K, w(n=1) = %12.8f eV, number of wn>0 numwi =%6d\n',T,wn,numwi)
fprintf('  initial gap = %g eV\n',gap)
%define the freqency grid
WN = pi*(2*(-numwi:(numwi-1))+1)/beta;
%setup the arrays to store the self-energy on the imaginary axis,
%length(WN) = 2*numwi
Z = zeros(2*nk,2*nk,length(WN));
X = zeros(2*nk,2*nk,length(WN));
P = zeros(2*nk,2*nk,length(WN));
Zold = Z;
Xold = X;
Pold = P;

%temporary matrices to hold the arrays on the imaginary axis.
mat1 = zeros(2*nk,2*nk,length(WN));
mat2 = zeros(2*nk,2*nk,length(WN));
mat3 = zeros(2*nk,2*nk,length(WN));
mat4 = zeros(2*nk,2*nk,length(WN));
mat5 = zeros(2*nk,2*nk,length(WN));
mat6 = zeros(2*nk,2*nk,length(WN));

%initialize self-energy
for nn = 1:length(WN)
    Z(:,:,nn) = WN(nn);
    X(:,:,nn) = 0;
    P(:,:,nn) = gap;
end

loop = 1;
iter = 0;
%
if delta_a2f_ph == 1
    fprintf('  NU is summed at NU =%7.3f eV.\n',wph)
end
if mixing == 1
    fprintf('  mixing_beta = %6.3f\n',mixing_beta)
end
fprintf(['  Iter.  _______dZ_______  _______dX_______  _______dP_______',...
    '  _______mu_______  _______<n>______\n'])
while(loop == 1)
    iter = iter + 1;
    %save the old self-energy
    if (mixing == 1) && (iter>1)
        Zold = mixing_beta*Z + (1-mixing_beta)*Zold;
        Xold = mixing_beta*X + (1-mixing_beta)*Xold;
        Pold = mixing_beta*P + (1-mixing_beta)*Pold;
    else
        Zold = Z;
        Xold = X;
        Pold = P;
    end
 
      %update ek at new mu for the fixed filling filling0
    ek = energy(KX,KY,mu);  % ********
    
 %********************** impurity ******************************

    for nn = 1:length(WN)
        gee0 = 1/pi/dos0/(4*nk*nk)* Zold(:,:,nn)./  ... 
            ( Zold(:,:,nn).^2 + (ek+Xold(:,:,nn)).^2 + Pold(:,:,nn).^2 );
        gee1 = 1/pi/dos0/(4*nk*nk)*  Pold(:,:,nn)./ ... 
             ( Zold(:,:,nn).^2 + (ek+Xold(:,:,nn)).^2 + Pold(:,:,nn).^2 );
         g0 = sum(sum(gee0));
         g1 = sum(sum(gee1));
    end
    %initialize with zeros or constants
    for nn = 1:length(WN)
        Z(:,:,nn) = WN(nn) + nimp*g0/cee/cee ;
        X(:,:,nn) = nimp/cee;    
        P(:,:,nn) = -nimp*g1/cee/cee;
    end
 %***************** end of impurity *******************************
 
 
 
  
    %Begin Matsubara frequency w-sum and k-sum. k-sum is done with FFT.

    %First, do the k-sum (FFT convolution). Build some look up tables to
    %store the convolved components of the Green's function.
    for mm = 1:length(WN)
        Den = Zold(:,:,mm).^2 + (ek(:,:) + Xold(:,:,mm)).^2 + Pold(:,:,mm).^2;
        gzee = Zold(:,:,mm)./Den;
        gchi = (ek(:,:) + Xold(:,:,mm))./Den;
        gphi = Pold(:,:,mm)./Den;
        Gzee = fft2(gzee);            %Store the Fourier transforms.
        Gchi = fft2(gchi);
        Gphi = fft2(gphi);
        mat1(:,:,mm) = fftshift(ifft2(Qsf.*Gzee));
        mat2(:,:,mm) = fftshift(ifft2(Qph.*Gzee));
        mat3(:,:,mm) = fftshift(ifft2(Qsf.*Gchi));
        mat4(:,:,mm) = fftshift(ifft2(Qph.*Gchi));
        mat5(:,:,mm) = fftshift(ifft2(Qsf.*Gphi));
        mat6(:,:,mm) = fftshift(ifft2(Qph.*Gphi));
    end
    %Second, do the Matsubara frequency wm-sum at each frequency wn
    for idx = 1:numwi
        nn = idx + numwi;
        wn = WN(nn);
        for mm = 1:length(WN)
            wm = WN(mm);
            %evaluate the nu part of the integral
            if delta_a2f_ph == 1
                % trapezoidal rule will smear out a delta-function like
                % phonon spectrum a2f_ph. Use right Riemann sum.
                Dsf = 0;    %because mat1, mat3, mat5 = 0
                Dph = wph * wph /(wph^2 + (wn-wm)^2);
            else
                Fsf = 2*NU.*a2f_sf./(NU.^2 + (wn-wm)^2);
                Dsf = trapz(NU,Fsf);
                Fph = 2*NU.*a2f_ph./(NU.^2 + (wn-wm)^2);
                Dph = trapz(NU,Fph);
            end
            Z(:,:,nn) = Z(:,:,nn) + Dsf*mat1(:,:,mm)/(4*nk*nk*beta) ...
                + Dph*mat2(:,:,mm)/(4*nk*nk*beta);
            X(:,:,nn) = X(:,:,nn) - Dsf*mat3(:,:,mm)/(4*nk*nk*beta) ...
                - Dph*mat4(:,:,mm)/(4*nk*nk*beta);
            P(:,:,nn) = P(:,:,nn) + Dsf*mat5(:,:,mm)/(4*nk*nk*beta) ...
                + Dph*mat6(:,:,mm)/(4*nk*nk*beta);
        end
        Z(:,:,numwi-idx+1) = -Z(:,:,nn);
        X(:,:,numwi-idx+1) = X(:,:,nn);
        P(:,:,numwi-idx+1) = P(:,:,nn);
    end

    Nk = [ nk  nk ];
    %Done k-sum and \nu-sum. Update chemical potential
    mu = get_mu(Z,X,P,WN,KX,KY,beta,filling0,1);  %  ****************
   % mu = mu0;  
    
    filling = filling_Imaxis(Z,X,P,WN,KX,KY,beta,mu);

    diffz = 0;
    diffp = 0;
    diffx = 0;
    for nn = 1:numwi
        diffz = max([diffz, max(max(abs(Z(:,:,nn)-Zold(:,:,nn))))]);
        diffx = max([diffx, max(max(abs(X(:,:,nn)-Xold(:,:,nn))))]);
        diffp = max([diffp, max(max(abs(P(:,:,nn)-Pold(:,:,nn))))]);
    end
    if (diffz < eps && diffx < eps && diffp < eps), loop = 0; end;
    %Output results in this loop
    fprintf('  %5d  %16.12f  %16.12f  %16.12f',iter,diffz,diffx,diffp)
    fprintf('  %16.12f  %16.12f\n',mu,filling)
    fprintf('       P = %16.12f\n',max(max(real(P(:,:,numwi+1)))))    
    %Decide if we need to exit the loop.
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check convergence!\n',maxiter)
    end    
end