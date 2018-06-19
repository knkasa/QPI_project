% Code that calculates electron self-energies in the superconducting state
% that arise from electron-phonon coupling.  Output can be used for the QPI
% calculation.

clc;
clear;

nk = 32  ;     %  (nk=32) # of k points  (-nk:nk-1)  
eps = 1e-6;    % convergence criteria
T = 10;     % temp *********
q0 = 0.3;   % forward scattering criteria, smaller the more forward 
lambda_ph = 0.8;  % lambda from phonon (forward scattering)
lambda_sf = 0.0;  % lambda from spin fluctuation

%define the frequency grid W on the real axis(make sure it covers bandwidth)
dw = 0.002  ;  % 0.001 default **********************
numw = 500;    % 600 default (note numw>numnu) ********************

%define the frequency grid NU on the real axis (subset of W)
dnu = dw;
numnu = 100;  % (default 110) make sure it covers all boson modes ******

timp = 0.003;   % smoothen DOS pick 0.001 for better result******

wph = 0.03;    % 13 phonon energy peak in a2F forward scattering (80mev)****
wsf = 0.0;    % sf energy [eV] peak in a2F (if lambda_sf =/= 0)

% impurity parameter 
cee = 10;  % 10 default scattering rates
nimp = 0.0;  % 0.1 default n_imp (set this to zero for no impurity)
dos0 = 1;    % dos at fermi level

%required constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxiter = 200;         % maximum iteration
kb = 8.617e-5;          %Boltzmann's constant [eV/K]
mu0 = -0.235     %(-1.485 for qpi-paper) initial chemical potential(default=-0.235)

sig = 0.02;      % 0.015 default fermi-surface broadening
Ncut = 6;               %Matsubara frequency cut-off = Ncut*\Omega_c
gap = 0.007;            %initial guess of the gap [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input parameters

%------------- pick ph scattering option for of your chice --------------
 vertexinput = {'exp',q0};  % default (forward scattering) ***********
% vertexinput = {'gauss',q0};
% vertexinput = {'kd',nk};
% vertexinput = {'window',nk,0.3};
 %vertexinput = {'const',nk};  % this is for uniform scattering

%--------------------------------------------------------------------------------
 
gphinput = 0.04;    %  [eV]
lambda_mass = 0.2;  %  0.2;

beta = 1/(kb*T);    
mixing = 1;
mixing_beta = 0.7;

delta_a2f_ph = 1 ;   % choose this for default **************
 %delta_a2f_ph = 1;  % choosing this to make a2F a delta function 

save_real_axis = 1;
save_imag_axis = 0;
Fixlambda = 'lambdaph';
% Fixlambda = 'gph';
% Fixlambda = 'lambdamass';
% gap = 0;

%define the momentum grid
K = (-nk:(nk-1))*pi/nk;
[KX, KY] = meshgrid(K);

W = (-numw:numw)*dw;
wscale = dw*numw;

NU = (1:numnu)*dnu;
nuscale = dnu*numnu;
if nuscale > wscale
    warning('nuscale is too large!')
end
%Note: Band width 8t = 8*0.075 eV = 0.6 eV. If the band crosses Fermi level
%(metal), -8t <= W_min < W_max = wscale <= 8t. dnu = dw to simplify the
%calculation on real axis.

%set energy grid and dk = delta function \delta(\epsilon_k)
ek = energy(KX,KY,mu0);
dk = gauss(ek,sig);
DK = fft2(dk);          %Fourier transform of dk
Nf = sum(sum(dk))/(4*nk*nk);            %DOS at the Fermi level.
%define the boson spectra: phonons and spin-fluctuations, normalized to
%\int d\nu (2/\nu) a2f(\nu) = 1

%*********  spin fluctuation spectral density  *******************
a2f_sf = (NU./(NU.^2 + wsf^2));  % spin fluctuation MMP model
%del = 0.001;
%x1 = 0.011;
%x2 = 0.021;
% two phonon modes.  uniform scattering
%a2f_sf = del/pi./( (NU-x1).^2 + del*del ) + del/pi./( (NU-x2).^2 + del*del );

a2f_sf = a2f_sf/trapz(NU,a2f_sf*2./NU); %normalizing to one
%*****************************************************************


if delta_a2f_ph == 1
    a2f_ph = (wph/2)*(NU==wph)/dnu;
else
    a2f_ph = gauss(NU-wph,0.001); %***  phonon peak width (0.001)*****
    a2f_ph = a2f_ph/trapz(NU,a2f_ph*2./NU);
end

%set up coupling vertex (wph/2)*|g(q)|^2 and take its Fourier transform
gph = 1;
%qph = gph*exp(-2*sqrt(KX.*KX+KY.*KY)/q0);
qph = gph*vertexq(KX,KY,vertexinput{:});
Qph = fft2(qph);
%calculate lambda_sf; use the fact that lambda_sf is defined in part by a
%convolution of the two delta functions.
mat = dk.*fftshift(ifft2(Qph.*DK))/((4*nk*nk)^2);
rtmp1 = sum(sum(mat))/Nf;
%set coupling vertex gph
% Electron-phonon coupling vertex |g(q)|^2 = g_0^2 f(q). f(q) is normalized
% to unity in full q space or in the Brillouin zone.
% g_0^2 = gph*(wph/2) is determined in three cases
% 'lambdaph'    : double-delta function averaged lambda_ph is fixed
% 'gph'         : coupling vertex g_0 itself is fixed
% 'lambdamass'  : \lambda_m = -Re\partial\Sigma(\omega,k)/\partial\omega is
% fixed. Approximate formula \lambda_m = g_0^2/wph^2 = gph/(2*wph) is used
% to determin gph.


if strcmp(Fixlambda,'lambdaph')
    gph = lambda_ph/rtmp1;
elseif strcmp(Fixlambda,'gph')
    gph = gphinput;
    lambda_ph = gph*rtmp1;
elseif strcmp(Fixlambda,'lambdamass')
    gph = lambda_mass*(2*wph);
    lambda_ph = gph*rtmp1;
else
    error(['Fixlambda method ',Fixlambda,' is not implemented.'])
end
%adjust coupling vertex to fix lambda and Fourier transform
qph = gph*vertexq(KX,KY,vertexinput{:});
Qph = fft2(qph);

%
%--------- phonon parameter --------------------------------------
   
   % qph = ones(size(KX))    ; % phonon coupling constant a(q)
    % qph = ( sin(KX/2).^2 + sin(KY/2).^2 ) ;  % breathing mode
     %qph = 1 + 0.5*(cos(KX) - cos(KY) ) ;
   qph = ( cos(KX/2).^2 + cos(KY/2).^2 ) ;  % buckling mode
    % qph = 2*pi/q0/q0* exp( -2*sqrt(KX.*KX+KY.*KY)/q0 ) ;  % forward scatering
   
   Qph = fft2(qph);

   gfourier = dk.*fftshift(ifft2(Qph.*DK))   ;
   norm_fac = sum(sum(gfourier)) /Nf /((4*nk*nk)^2)   ;
   
   qph = lambda_ph*qph/norm_fac   ;  
   Qph = fft2(qph) ; 
   
%-------------------------------------------------------------------
%}




%plot( linspace(-pi,pi,256) , qph(:,128) )  % for nk=128 (for example)



%set up coupling vertex (wsf/2)*|g(q)|^2 and take its Fourier transform
gsf = 1;
qsf = gsf*ones(size(KX));
Qsf = fft2(qsf);
%calculate lambda_sf
mat = dk.*fftshift(ifft2(Qsf.*DK))/((4*nk*nk)^2);
rtmp2 = sum(sum(mat))/Nf;
%set coupling vertex gsf
gsf = lambda_sf/rtmp2;
%adjust coupling vertex to fix lambda and Fourier transform
qsf = gsf*ones(size(KX));
Qsf = fft2(qsf);

%get the initial filling
filling0 = 2*sum(sum(fermi(energy(KX,KY,mu0),beta)))/(4*nk*nk);
mu = mu0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileDir = './';  % change this to your directory*********
filnamestr = ['_gph=' num2str(gph) '_lamph=' num2str(lambda_ph) ...
    '_q0=' num2str(q0) '_T=' num2str(T) '_nk=' num2str(nk) '.dat'];
filAkpath = ['Akpath' filnamestr];
filAkGamm = ['AkGamm' filnamestr];
filDOS = ['DOS' filnamestr];
filout = ['out' filnamestr];
fid = fopen([fileDir,filout],'w');
fclose(fid);
%
diary([fileDir, filout])
%print input parameters
fprintf('\n')
fprintf('spectral function calculation\n')
fprintf('  grid nk = %4d, convergence criterion = %.2e, smearing = %g eV\n',nk,eps,sig)
fprintf('  wph = %12.8f eV,    wsf = %12.8f eV\n',wph,wsf)
fprintf('  gph = %12.8f eV,    gsf = %12.8f eV\n',gph,gsf)
fprintf('  lambda_ph = %12.8f, lambda_sf = %12.8f\n',lambda_ph,lambda_sf)
fprintf('  q0 = %g, vertex f(q) = ''%s''\n',q0,vertexinput{1})
fprintf('  mu0= %g eV, N(E_f) = %12.8f (1/eV/per spin/Volume)\n',mu0,Nf)
fprintf('  mu = %g eV, filling0 = %12.8f\n',mu,filling0)
fprintf('  numw =%6d, numnu =%6d, maxiter =%6d\n',numw,numnu,maxiter)
if save_real_axis == 1
    filraxis = ['real_axis'  filnamestr(1:(end-4)) '.mat'];
    var_raxis = {'nk','T','mu','W','ek','Z','X','P'};
    fprintf('  Real axis self-energy will be saved.\n')
end
if save_imag_axis == 1
    filraxis = ['imag_axis'  filnamestr(1:(end-4)) '.mat'];
    var_iaxis = {'nk','T','mu','WN','ek','Z','X','P'};
    fprintf('  Imaginary axis self-energy will be saved.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Imaginary axis calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
calculate_im_axis;
t(1) = toc;
fprintf('Done imaginary axis calculation. Time = %.2f s\n',t(1));
if save_imag_axis == 1
    save([fileDir, filraxis],var_iaxis{:},'-mat')
    fprintf(['  ',fileDir,filraxis,' saved.\n'])
end


wm = pi/beta;
gap0 = max(max(wm*P(:,:,numwi+1)./Z(:,:,numwi+1)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calculation of the driving term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
calculate_driving_terms;
t(2) = toc;
fprintf('Done driving term calculation. Time = %.2f s\n',t(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Calculation of the real-axis term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
calculate_real_axis;
t(3) = toc;
fprintf('Done real axis calculation. Time = %.2f s\n',t(3));
if save_real_axis == 1
    save([fileDir, filraxis],var_raxis{:},'-mat')
    fprintf(['  ',fileDir,filraxis,' saved.\n'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate lambdam_avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ek = energy(KX,KY,mu);  %at final mu
dk = gauss(ek,sig);
Nf = sum(sum(dk))/(4*nk*nk);
iw0 = find(W==0);
zk = real(Z(:,:,iw0+1)-Z(:,:,iw0-1))/(W(iw0+1)-W(iw0-1))-1;
lambdam_avg = sum(sum(zk.*dk))/sum(sum(dk));
fprintf('lambdam_avg = %12.8f\n',lambdam_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the final results to the various data files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
output_results;
t(4) = toc;
fprintf('Done save data. Time = %.2f s\n',t(4));
fprintf('Done spectral function calculation. Total Time = %.2f s\n',sum(t));
diary off

%}


