function [ Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S,dgiven ] = solve_model( X, Y, params )
%solve_model Function that takes X, Y and the parameters and returns H(k*)
%   Solving the model recursively starting from the topmost income class

a = params(1);
al = params(2);
n = length(X); %Number of income classes

%H distribution for decision of sex selection
H = makedist('Uniform','lower',0,'upper',a);
%H = makedist('Normal','mu',0.5,'sigma',0.2);

x = n;
y = n;
i = 1;
d = (0.4)*Y(y,1);
t = (X(x,1)+d)/2;
Cx = X(x,1) + d - t;
cx = al*t;
vin = log(Cx) + log(cx);

v_high = log((9*al/16)*X(x,1)*X(x,1));
v_low = log((1/4)*X(x,1)*X(x,1));
% disp(v_low)
% disp(vin)
% disp(v_high)

vxx = zeros(n,1);
vxx(n,i) = vin;
i = i+1;
vxx(n,i) = vin;

while(1)
%     vin = vxx(n,i-1);
%     temp = @(vi) real(set_vin(n,H,al,vi,X,Y,i,vxx));
    temp = @(vi) set_vin(n,H,al,vi,X,Y,i,vxx);
    %disp(strcat('i=',num2str(i)))
    options = optimset('Display','off');
%     ind = 1;
%     for i2 = v_low:0.001:v_high
%         vals(ind) = (temp(i2));
%         ind = ind+1;
%     end
%     vdom = v_low:0.001:v_high;
%     plot(vdom,vals)
%     pause;
%     [vinit,fval] = fzero(temp,vin,options);
    [vinit,fval] = fzero(temp,[v_low,v_high],options);
    if (fval > 10^-2)
        disp('boundary condition not met')
        disp(fval)
    end
    if (vinit > v_high || vinit < v_low)
        disp('wrong solution')
    end
    
    [~,xend,vxx,Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,dy,S,dgiven ] = set_vin(n,H,al,vinit,X,Y,i,vxx);
    while (xend ~= 0)
        vxx(xend,i) = S(xend,1);
        CXx(xend,2) = my_psi(S(xend,1),al);
        cXx(xend,2) = al*my_psi(S(xend,1),al);
        d = 2*(my_psi(S(xend,1),1)) - X(xend,1);
        dy(xend,2) = d;
        xend = xend-1;
    end
    i = i+1;
    if (i>2 & (vxx(:,i-1) == vxx(:,i-2) | i > 35))
        break
    end
end
Hky = Hky(:,2);
end