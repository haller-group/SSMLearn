function [V,D,d_cont] = eigSorted(A,varargin)
% Sorting eigenvalues as cc with positive imaginary parts, reals and cc 
% with negative imaginary parts. Frrequencies are sorted from slowest to
% fastest. varargin is an optional input for discrete dynamical systems,
% which indicates the sampling time as the ordering is based on the
% continuous eigenvalues.

k = size(A,1);
[V,D] = eig(A);
if isempty(varargin) == 1
   d_cont = diag(D);
else
   Dt = [varargin{:}]; d_cont = log(diag(D))/Dt;
end
n_re_eig = sum(imag(d_cont)==0);
n_cc_eig = (k-n_re_eig)/2;
[~,pos] = sort(imag(d_cont),'descend'); pos = transpose(pos);
if n_re_eig == 0
    pos_sort = [fliplr(pos(1:n_cc_eig)) pos(n_cc_eig+[1:n_cc_eig])];
else
    pos_sort = [fliplr(pos(1:n_cc_eig)) pos(n_cc_eig+[1:n_re_eig]) ...
                                      pos(n_cc_eig+n_re_eig+[1:n_cc_eig])];
end
d = diag(D); d = d(pos_sort); D = diag(d); V = V(:,pos_sort);
d_cont = d_cont(pos_sort);
end