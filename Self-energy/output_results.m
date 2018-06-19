%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output A(k,w) along cut (-pi,0)--(+pi,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxomg =  dw*numw;   % max omega for plot (-maxomg:maxomg) **********

fid = fopen([fileDir,filAkpath],'w');
indW = sort(find(W>=-maxomg & W<=maxomg),'ascend');  % -0.6~0.5 default*********
[KK WW] = meshgrid(K,W(indW));
%Akw = zeros(size(KK));
%nky = find(K == 0);
for nkx = 1:numel(K)
     kx = K(nkx);
    for nky = 1:numel(K)
     ky = K(nky);
   % ky = 0;
   % ek = energy(kx,ky,mu);
    for nw = indW
        ek = energy(kx,ky,mu);
        
        w = W(nw);
        zee = Z(nky,nkx,nw);
        chi = X(nky,nkx,nw);
        phi = P(nky,nkx,nw);
        %fermi = 1/(exp(beta*w)+1);
       % Denom = zee^2 - (ek+chi)^2 - phi^2;
      %  Green = (zee + ek + chi)/Denom;
       % Akw(nw,nkx) = -2*imag(Green)*fermi/pi;
        fprintf(fid,[repmat('%16.10f ',1,9), '\n'],...
            kx, ky, w,  ...
            real(zee),imag(zee),real(chi),imag(chi),real(phi),imag(phi));
    end
    end
end
fclose(fid);
fprintf(['  ',fileDir,filAkpath,' saved.\n'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output A(k = 0,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([fileDir,filAkGamm],'w');
indW = sort(find(W>=-maxomg & W<=maxomg),'ascend');   % -0.6~0.5 default*********
nk = find(K == 0);
kx = 0;
ky = 0;
ek = energy(kx,ky,mu);
Akw = zeros(size(W(indW)));
for nw = indW
    w = W(nw);
    zee = Z(nk,nk,nw);
    chi = X(nk,nk,nw);
    phi = P(nk,nk,nw);
    fermi = 1/(exp(beta*w)+1);
    Denom = zee^2 - (ek+chi)^2 - phi^2;
    Green = (zee + ek + chi)/Denom;
    Akw(nw) = -2*imag(Green)*fermi/pi;
    fprintf(fid,[repmat('%16.10f ',1,11), '\n'],...
        kx, ky, w, Akw(nw), ek,...
        real(zee),imag(zee),real(chi),imag(chi),real(phi),imag(phi));
end
fclose(fid);
fprintf(['  ',fileDir,filAkGamm,' saved.\n'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output the Density of States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([fileDir,filDOS],'w');
ek = energy(KX,KY,mu);
mu  % chemical potential
timp

const = complex(0,ones(64,64)) ;

for nw = 1:numel(W)
    w = W(nw);
    zee = Z(:,:,nw) ;
    chi = X(:,:,nw);
    phi = P(:,:,nw);
    Denom = ((zee-complex(0,0.000))).^2 - (ek+chi).^2 - phi.^2;
    Green = (zee + ek + chi)./Denom;
    DOS = -2*sum(sum(imag(Green)))/(4*nk*nk*pi);
    fprintf(fid,'%16.10f %16.10f\n', w,DOS);
end
fclose(fid);
fprintf(['  ',fileDir,filDOS,' saved.\n'])