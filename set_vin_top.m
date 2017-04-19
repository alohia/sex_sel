function [ diff, xend, vxx, Hky, ks, CYy, cYy, CXx, cXx, ky, vyy, dy, S ] = set_vin_top( n, H, al, vin, X, Y, i, vxx, tr_in )
%set_vin Function to set initial value of v(x)
%   This function solves the model to get a value of v(x) for the last
%   matched boy and the income xbar of this boy.

S(:,1) = (log((X(:,1).^2)/4));
CXx(:,1) = X(:,1);
cXx(:,1) = X(:,1);
CYy(:,1) = Y(:,1);
cYy(:,1) = Y(:,1);
ky(:,1) = Y(:,1);
Hky(:,1) = Y(:,1);
vyy(:,1) = Y(:,1);
dy(:,1) = X(:,1);
Hks = zeros(n,n);
ks = zeros(n,n,1);


xi = 2;
x = n;
y = n;
xchng = 0;
del = (vin-S(1,1))/(n-1);
if (i==2)
    vxx(x,i-1) = vin;
    for b = n-1:-1:1
        %vxx(b,i-1) = vxx(b+1,i-1) -  my_vp(vxx(b+1,i-1),X(b+1,1),Y(b,1),sig,bet,al)*(X(b+1,1) - X(b,1));
        vxx(b,i-1) = vxx(b+1,i-1)-del;
    end
end
vxx(x,i) = vin;
d = 2*(my_psi(vxx(x,i),al)) - X(x,1);
dy(x,2) = d;
t = (X(x,1)+d)/2;
Cx = X(x,1) + d - t;
cx = al*t;
CXx(x,2) = Cx;
cXx(x,2) = cx;

Cy = Y(y,1) - d;
cy = (1-al)*t;
CYy(y,xi) = Cy;
cYy(y,xi) = cy;
vyy(y,2) = log(Cy) + log(cy);

exitflag = 0;

while Y(y,2)>=0
    y = y+1;
        
    while (X(x,2)>0) % Going down to the next girl as long as boys are left
        y = y-1;
        if (y<=2)
            tr=tr_in;
        else
            tr=0;
        end
        if (y==0)
            break
        end
            
        if (~xchng) % No change in x means we go down for the next match
            Cy = Y(y,1) + tr - d;
            cy = (1-al)*t;
            vyy(y,2) = log(Cy) + log(cy);
            kstr = vxx(y,i-1) - vyy(y,2);
            Hk = cdf(H,kstr);
            X(y,2) = X(y,2) + Hk*Y(y,2);
            Y(y,2) = Y(y,2) * (1 - Hk);
        else
            xchng = 0;
        end
            
        CYy(y,xi) = Cy;
        cYy(y,xi) = cy;
        ky(y,xi) = kstr;
        Hky(y,xi) = Hk;
        ks(y,x,i) = kstr;
        Hks(y,x) = Hk;
        bleft = X(x,2);
        X(x,2) = X(x,2) - Y(y,2);
        if (X(x,2) > 0)
            xi = 2;
        end
    end

    if (y==0)
        vend = log(Cx)+log(cx);
        xend = x;
        break
    end
    
    if (y < n)
        xi=xi+1;
        Y(y,2) = Y(y,2) - bleft;       
        xchng = 1;
        x = x-1;
        if (x ~= n)
            v_high = log((X(x,1)+Y(y,1))/2) + log(al*((X(x,1)+Y(y,1))/2));
            v_low = log(X(x,1)/2) + log((X(x,1)/2));
            temp = @(v) real(my_phi(v,al,vyy(y,2),X(x,1),(Y(y,1)+tr)));
            opts1=  optimset('display','off');
%             disp('temp v_low in setvin')
%             temp(v_low)
%             disp('temp v_high in setvin')
%             temp(v_high)
%             if (temp(v_high-0.0001)*temp(v_low+0.0001) > 0)
%                 exitflag = 1;
%                 break
%             end
%             vnew = lsqnonlin(temp,vxx(x+1,i),v_low,v_high,opts1);
%             vdom = v_low:0.001:v_high;
%             vals = temp(vdom);
%             plot(vdom,vals)
            vnew = fsolve(temp,vxx(x+1,i),opts1);
            %vnew = fzero(temp,[v_low+0.0001,v_high-0.0001],opts1);
%             syms v_val;
%             phi_eqn = vyy(y,2) == log(Y(y,1) - 2*(my_psi(v_val,al)) + X(x,1)) + log((1-al)*(my_psi(v_val,al)));
%             vnew = vpasolve(phi_eqn,v_val)
            %if (vnew > v_high || vnew < v_low)
                %disp('wrong solution in setvin')
%                 vdom = v_low:0.001:v_high;
%                 vals = temp(vdom);
%                 plot(vdom,vals)
%                 pause;
           % end
            vxx(x,i) = vnew;
            %disp(vxx(x,i) < log((al*(X(x,1)+Y(y,1)))/4));
            %disp(vxx(x,i) > log((al*(X(x,1)+Y(y,1)))/16));
            d = 2*(my_psi(vxx(x,i),al)) - X(x,1);
            %disp((Y(y,1)-d)/((Y(y,1)-d)-(X(x,1)+d)) < 0);
            dy(x,2) = d;
            t = (X(x,1)+d)/2;
            Cx = X(x,1) + d - t;
            cx = al*t;
            CXx(x,2) = Cx;
            cXx(x,2) = cx;
        end
    end
end
% if exitflag == 1
%     diff = -1000;
%     return 
% end
diff = vend - S(xend);
end