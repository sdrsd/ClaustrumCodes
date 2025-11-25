
function [r_all] = simulate_network(I, W, T, dt, tau, th)

if nargin < 6; th = 0; end

N = size(I,1);

r = zeros(1,N);
r_all = zeros(N,length(T));
for i = 1:length(T)-1
    inp_rec = r*W;
    inp_ffw = I(:,i)';
    
    dr = dt/tau * (-r + rectify(inp_rec+inp_ffw -th));
            
    %size(dr)
    
    r = r + dr;
    
    r_all(:,i+1) = r;
end
    
end

