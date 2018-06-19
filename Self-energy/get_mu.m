function out = get_mu(Z,X,P,WN,KX,KY,beta,fill,Imaxis)
%set the initial values of the chemical potential.
muL = energy(0,0,0)   ;       
muR = energy(pi,pi,0) ;
muM = (muR+muL)/2; 

%Get the fillings
if(Imaxis == 1),
  fillL = filling_Imaxis(Z,X,P,WN,KX,KY,beta,muL);
  fillR = filling_Imaxis(Z,X,P,WN,KX,KY,beta,muR);
  fillM = filling_Imaxis(Z,X,P,WN,KX,KY,beta,muM);
else
  fillL = filling_Reaxis(Z,X,P,WN,KX,KY,beta,muL);
  fillR = filling_Reaxis(Z,X,P,WN,KX,KY,beta,muR);
  fillM = filling_Reaxis(Z,X,P,WN,KX,KY,beta,muM);
end; 

while(abs(fillM-fill)>1e-6),
  if(fillM > fill),
      fillR = fillM;
      muR = muM;
  elseif(fillM < fill)
      fillL = fillM;
      muL = muM;
  end;
  muM = (muR+muL)/2; 
  %get the new filling at the mid point.
  if(Imaxis == 1),
    fillM = filling_Imaxis(Z,X,P,WN,KX,KY,beta,muM);
  else
    fillM = filling_Reaxis(Z,X,P,WN,KX,KY,beta,muM);
  end;
end;


out = muM;