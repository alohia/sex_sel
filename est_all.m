function [ SSR ] = est_all( castes,params )
%Returns the sum of squared errors across all castes
%   This function solves the model for the given parameters for each caste
%   and gives the sum of squared errors across all castes

%num_castes = numel(castes(:,1));
num_castes = numel(castes);
SSR = 0;
for i = 1:num_castes
    [X, Hky ,n] = getdata(castes(i));
    %[X, Hky ,n] = getdata(castes(i,1));
    %X(:,1) = X(:,1)*10;
    Y = X;
    SSR = SSR + sumsqr(Hky - solve_model(X,Y,params));
%     SSR = SSR + castes(i,2)*sumsqr(Hky - solve_model(X,Y,params));
end
end

