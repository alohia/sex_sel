function [ diff ] = my_phi_dtax( v,al,p,x,y,theta )
%phi function
%   Gives the value for phi(v) - phi to find the value for v
phi_val = log(y - (2/theta)*(my_psi(v,al)) + (x/theta)) + log((1-al)*(my_psi(v,al)));
diff = phi_val - p;
end