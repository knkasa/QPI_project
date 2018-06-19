function out = filling_Imaxis(Z,X,P,WN,KX,KY,beta,mu);
ek = energy(KX,KY,mu);
green = zeros(size(ek));
den = zeros(size(ek));
N = size(ek);
i = sqrt(-1);
fill = 0;
fp = fermi(ek,beta);

fill = sum(sum(fp));
for nn = 1:length(WN),
    wn = WN(nn);
    den = (i*wn - ek(:,:)).*(ek(:,:) - i*Z(:,:,nn) + X(:,:,nn));    
    green = (i*wn - i*Z(:,:,nn) + X(:,:,nn))./den;
    fill = fill - sum(sum(green/beta));
end;

clear green;
clear den;

%factor of two is for spin.

out = real(2*fill/(N(1)*N(2)));