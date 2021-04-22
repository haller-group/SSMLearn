function [DS, S, mfd] = getSSM(M, C, K, fnl, SSMDim, varargin)

mfdOrder = 9;
if ~isempty(varargin)
    mfdOrder = varargin{:};
end

DS = DynamicalSystem();
set(DS, 'M', M, 'C', C, 'K', K, 'fnl', fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
S.choose_E([1:SSMDim]);
mfd = S.compute_whisker(mfdOrder);