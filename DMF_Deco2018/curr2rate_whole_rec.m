function f_rate = curr2rate_whole_rec(x,wgain,g,I,c,receptors)
% Computes transfer function from unit current to unit firing rate by a
% nonlinear function
%
% From Deco et al 2014.
% g is d
% I is b
% c is a

y=bsxfun(@times,(c.*x-I),(1+receptors*wgain)');
if y~=0
    f_rate = y./(1-exp(-g.*y));
else
    f_rate=0;
end
end