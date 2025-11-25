
function [W, Dm] = dd_w_fun(x1, x2, L, A, sig_d, d0)

if nargin < 6; d0 = 0; end

if isrow(x1); x1 = x1'; end
if isrow(x2); x2 = x2'; end

n_sh = round(d0 / L *length(x1));

D = pdist2(x1, x2);

Dm = D;
Dm(Dm>(L/2)) = L - Dm(Dm>(L/2));

W = exp(-Dm/sig_d);

W = W / mean(sum(W)) * A;

W = W .* rand(size(W));

W = circshift(W, -n_sh);

end
