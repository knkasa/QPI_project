fprintf('\n')
fprintf('**************************************************************\n')
fprintf('Begin driving term calculation.\n')
%
if delta_a2f_ph == 1
    fprintf('  NU is summed at NU =%7.3f eV.\n',wph)
end

%set arrays to store the driving terms
Cz = zeros(2*nk,2*nk,length(W));
Cx = zeros(2*nk,2*nk,length(W));
Cp = zeros(2*nk,2*nk,length(W));


i = sqrt(-1);
smear = timp;
%small positive imaginary part to avoid the poles on real axis

%temporary matrices to hold the arrays on the imaginary axis.
%size = [2*nk, 2*nk, length(WN)], length(WN) < length(W)
mat1 = zeros(size(Z));
mat2 = zeros(size(Z));
mat3 = zeros(size(Z));
mat4 = zeros(size(Z));
mat5 = zeros(size(Z));
mat6 = zeros(size(Z));

for mm = 1:length(WN)
    Den = Z(:,:,mm).^2 + (ek + X(:,:,mm)).^2 + P(:,:,mm).^2;
    gzee = Z(:,:,mm)./Den;
    gchi = (ek(:,:)+X(:,:,mm))./Den;
    gphi = P(:,:,mm)./Den;
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

for idx = 0:numw
    nw = idx + numw + 1;
    w = W(nw);
    % Cz(:,:,nw) = w + i*smear;
    
    %small array to speed up the loop. Cz is large when numw is large
    %%{
    tmpCz = (w + i*smear)*ones(2*nk, 2*nk);
    tmpCx = zeros(2*nk, 2*nk); 
    tmpCp = zeros(2*nk, 2*nk);
    %}
    for mm = 1:length(WN)
        iwm = i*WN(mm);
        wm = WN(mm);
        %evaluate the nu part of the integral
        if delta_a2f_ph == 1
            clam_ph = wph * wph /(wph^2 - (w - iwm)^2);
            clam_sf = 0;    %because mat1, mat3, mat5 = 0
        else
            F = 2*a2f_ph.*NU./(NU.*NU - (w - iwm).^2);
            clam_ph = trapz(NU,F);
            F = 2*a2f_sf.*NU./(NU.*NU - (w - iwm).^2);
            clam_sf = trapz(NU,F);
        end
        %{
        Cz(:,:,nw) = Cz(:,:,nw) + i*clam_sf*mat1(:,:,mm)/(4*nk*nk*beta) ...
            + i*clam_ph*mat2(:,:,mm)/(4*nk*nk*beta);
        Cx(:,:,nw) = Cx(:,:,nw) - clam_sf*mat3(:,:,mm)/(4*nk*nk*beta) ...
            - clam_ph*mat4(:,:,mm)/(4*nk*nk*beta);
        Cp(:,:,nw) = Cp(:,:,nw) + clam_sf*mat5(:,:,mm)/(4*nk*nk*beta) ...
            + clam_ph*mat6(:,:,mm)/(4*nk*nk*beta);
        %}
        
        %%{
        tmpCz = tmpCz + i*clam_sf*mat1(:,:,mm)/(4*nk*nk*beta) ...
            + i*clam_ph*mat2(:,:,mm)/(4*nk*nk*beta);
        tmpCx = tmpCx - clam_sf*mat3(:,:,mm)/(4*nk*nk*beta) ...
            - clam_ph*mat4(:,:,mm)/(4*nk*nk*beta);
        tmpCp = tmpCp + clam_sf*mat5(:,:,mm)/(4*nk*nk*beta) ...
            + clam_ph*mat6(:,:,mm)/(4*nk*nk*beta);
        %}
    end
    %%{
    Cz(:,:,nw) = tmpCz;
    Cx(:,:,nw) = tmpCx;
    Cp(:,:,nw) = tmpCp;
    %}
    Cz(:,:,numw+1-idx) = -conj(Cz(:,:,nw));
    Cx(:,:,numw+1-idx) = conj(Cx(:,:,nw));
    Cp(:,:,numw+1-idx) = conj(Cp(:,:,nw));
end
%Note: Cx, Cp should be real, imag(Cz) = smear
%conj to ensure the frequency symmetry is correct even when Cx, Cp are
%not exact real due to very small numerical errors.