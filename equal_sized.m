function [ X ] = equal_sized( X,nclass,pop )
%equal_sized function
%   Gives equal sized income classes for a given income distribution
n = length(X);
Y = zeros(nclass,1);
Y(:,1) = 0;
Y(:,2) = pop/nclass;
tmp = 0;
reqd = 0;
ind = 1;
l2=1;
while (ind<nclass+1)
    reqd = pop/nclass - tmp;
    tmp = tmp+X(l2,2);
    if (tmp < pop/nclass)
        Y(ind,1) = Y(ind,1) + (X(l2,2)*X(l2,1));
    elseif (tmp >= pop/nclass)
        Y(ind,1) = Y(ind,1) + (reqd * X(l2,1));
        X(l2,2) = X(l2,2) - reqd;
        Y(ind,1) = Y(ind,1)/(pop/nclass);
        ind = ind+1;
        tmp = 0;
        reqd = 0;
        l2 = l2-1;
    end
    l2=l2+1;
end
X=Y;
end