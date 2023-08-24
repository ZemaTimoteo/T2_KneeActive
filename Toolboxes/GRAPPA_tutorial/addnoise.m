function [out_mat] = addnoise(in_mat,level)

% addnoise
% This function adds complex gaussian noise to a matrix.  The scale value
% scales the noise.
%
% Felix Breuer 28.06.02


rand_real = randn(size(in_mat));
rand_imag = randn(size(in_mat));
out_mat = double(in_mat) + level*(rand_real + i*rand_imag);


