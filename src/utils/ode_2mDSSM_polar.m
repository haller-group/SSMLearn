function y = ode_2mDSSM_polar(z, p, data)
% ODE_2MDSSM_POLAR This function presents vectorized implementation of vector
% field of reduced dynamics on 2m-dimensional SSMs. Here z is a 
% 2m-dimensinoal state vector and p is parameter vector for excitation
% frequency and amplitude. All other info such as eigenvalues of master
% subspace, coefficients of nonlinear terms is included in the structure
% data. The state vector here is in the form of polar coordinate.
%
% See also: ODE_2MDSSM_CARTESIAN

assert(~isempty(data), 'Structure data in ode_2mDSSM_polar is empty');
% extract data fields
beta   = data.beta;
kappa  = data.kappa;
lamdRe = data.lamdRe;
lamdIm = data.lamdIm;
mFreqs = data.mFreqs;
iNonauto = data.iNonauto;
rNonauto = data.rNonauto;

% rename state and parameter
rho = z(1:2:end-1,:);
th  = z(2:2:end,:);
om   = p(1,:);
epsf = p(2,:);

% autonomous linear part
yrho = lamdRe.*rho;
yth  = lamdIm-mFreqs(:)*om;

% autonomous nonlinear part
m  = numel(mFreqs);
em = eye(m);
for i=1:m
    kappai = kappa{i};
    kappai = full(kappai);
    betai  = beta{i};
    nka = size(kappai,1);
    nbe = numel(betai);
    assert(nka==nbe, 'Size of kappa%d and beta%d does not match',i,i);
    ei = em(i,:);
    for k=1:nka
        ka = kappai(k,:);
        be = betai(k);
        l = ka(1:2:end-1);
        j = ka(2:2:end);
        ang = (l-j-ei)*th;
        rhopower = rho.^((l+j)');
        pdrho = prod(rhopower,1);
        yrho(i,:) = yrho(i,:)+pdrho.*(real(be)*cos(ang)-imag(be)*sin(ang));
        yth(i,:)  = yth(i,:)+pdrho./rho(i,:).*(real(be)*sin(ang)+imag(be)*cos(ang));
    end
end

% nonautonomous leading part
for i=1:numel(iNonauto)
    id = iNonauto(i);
    r  = rNonauto(i);
    r  = epsf*r; 
    rRe = real(r);
    rIm = imag(r);
    yrho(id,:) = yrho(id,:)+rRe.*cos(th(id,:))+rIm.*sin(th(id,:));
    yth(id,:)  = yth(id,:)-rRe.*sin(th(id,:))./rho(id,:)+rIm.*cos(th(id,:))./rho(id,:);
end

nt = numel(om);
y  = zeros(2*m,nt);
y(1:2:end-1,:) = yrho;
y(2:2:end,:)   = yth;

end