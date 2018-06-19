function out = energy(kx,ky,mu)


t = 0.38  ;   % t=0.38(for qpi-paper) default positive (0.075) ******
 %t = 0.075 ;

  out = -2*t*(cos(kx)+cos(ky))-mu  ;   

%{
t1=   0.36;  t2=-0.1;  t3=0.035;  t4=0.01;    tz=0.036;  tbi=0.11;  a0=0.4;
%t1 = 0.126;  t2=-0.036;  t3=0.015;  t4=0.015;  tz=0.005;  tbi=0.03; a0=0.4;  
Az = 4*tz*cos(0.5*kx).*cos(0.5*ky) ;
Tz = -sqrt( tbi^2 + Az.^2 ) ;
Ez = -Tz.*(  0.25*(cos(kx)-cos(ky)).^2 + a0 ) ;
out = -2*t1*( cos(kx)+cos(ky) ) - 4*t2*cos(kx).*cos(ky) - 2*t3*( cos(2*kx) + ...
    cos(2*ky) ) - 4*t4*( cos(2*kx).*cos(ky) + cos(2*ky).*cos(kx) ) + Ez +0.0 ;
%}

%  t1 = -0.5908; t2 = 0.0962; t3 = -0.1306; t4 = -0.0507; t5 = 0.0939;    %units: eV
%  eps0 = 0.1;
%  out =   0.5*t1*(cos(kx)+cos(ky)) + t2*cos(kx).*cos(ky) ...
%   + 0.5*t3*(cos(2*kx)+cos(2*ky)) + 0.5*t4*(cos(2*kx).*cos(ky) + cos(2*ky).*cos(kx)) ...
%    + t5*cos(2*kx).*cos(2*ky) + eps0;



