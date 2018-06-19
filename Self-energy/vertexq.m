function vq = vertexq(KX,KY,varargin)
%
funtype = varargin{1};
if strcmp(funtype,'exp')
    q0 = varargin{2}/2;
    vq = 4*(pi^2)*exp(-sqrt(KX.*KX+KY.*KY)/q0)/(2*pi*q0*q0);
elseif strcmp(funtype,'gauss')
    q0 = varargin{2};
    vq = 4*(pi^2)*exp(-(KX.*KX + KY.*KY)/(2*q0*q0))/(2*pi*q0*q0);
elseif strcmp(funtype,'kd')
    nk = varargin{2};
    vq = ((KX.*KX+KY.*KY)==0) * ((2*nk)^2);
elseif strcmp(funtype,'window')
    nk = varargin{2};
    delq = varargin{3};
    vqmat = (KX.*KX+KY.*KY)<=(delq^2);
    Nq = numel(find(vqmat));
    vq = vqmat*((2*nk)^2)/Nq;
    
elseif strcmp(funtype,'const') % uniform scattering ***************
    nk = varargin{2};
    vq = ones(2*nk,2*nk);  %  *(4*nk*nk);



else
    error([funtype, 'is not implemented!'])
end