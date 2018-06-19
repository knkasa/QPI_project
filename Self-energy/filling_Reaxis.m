function out = filling_Reaxis(Z,X,P,W,KX,KY,beta,mu)
ek = energy(KX,KY,mu);
Den = zeros(size(KX));
dw = W(2)-W(1);

fill = 0;
for nw = 1:length(W),
  w = W(nw);
  fermi = 1/(exp(beta*w)+1);
  Den = (Z(:,:,nw).^2 - (ek(:,:)+X(:,:,nw)).^2 - P(:,:,nw).^2);
  Green = (Z(:,:,nw) + ek(:,:) + X(:,:,nw))./Den;
  DOS =-2*sum(sum(imag(Green)))/(length(KX)*length(KX)*pi);
  fill = fill + fermi*dw*DOS;
end;

out = fill;