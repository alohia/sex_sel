function [ diff ] = my_phi_girls( v,al,p,x,y,tr )
%phi function
%   Gives the value for phi(v) - phi to find the value for v
phi_val = log(y - 2*(my_psi(v,al)) + x) + log((1-al)*(my_psi(v,al)) + tr);
diff = phi_val - p;
end