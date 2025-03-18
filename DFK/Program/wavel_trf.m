function [psi_jk] = wavel_trf_v2(j, t_k, t)
%create wavelet transform

d = 2^(j/2); % dilation

psi_jk = basic_wavel( d*d*t - t_k );
psi_jk = d*psi_jk;

end