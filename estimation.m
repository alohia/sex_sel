clc;
%clearvars;
%castes = [1,5,7,10,14,17,23,28,33,38,44];
% anom_castes = [4,5,7,44];

% castes = [1 .0984789; 4 0.0090223; 5 0.018784; 7 0.0151011; 10 .0221519; 14 .0589121; 17 .0112403; 23 .0183048; 28 .4044578; 33 .0585972; 28 .2693966; 44 .0155529];
% caste = 23;

% [X, Hky, n] = getdata(caste);
% %X = equal_sized(X,n,sum(X(:,2)));
% X(:,1) = X(:,1)*10;
% Y = X;

% Parameters guess - 10,0.75 are the correct params for caste 100
a = 12.; %upper limit for H distribution
al = 0.6; %exogenous alpha - husband's share
params= [a,al];
% index = 1;

pguess = params;
% options = optimoptions(@lsqnonlin,'Display','iter','OptimalityTolerance',1e-10);

% temp3 = @(pars) (Hky - solve_model(X,Y,pars));
% [par,res] = lsqnonlin(temp3, pguess, [3,0.3], [60,1], options)
% par = fminsearch(temp, pguess, options)

% estimates(:,1) = castes;

delete(gcp('nocreate')) % Delete any existing parallel pool
parpool % Start a parallel pool


% for caste = castes
%     disp(caste)
%     [X, Hky, n] = getdata(caste);
%     Y = X;
%     temp4 = @(pars) (Hky - solve_model(X,Y,pars));
%     problem = createOptimProblem('lsqnonlin','x0',pguess,'objective',temp4,'lb',[3,0.2],'ub',[85,1]);
%     ms = MultiStart('UseParallel',true);
%     [xmin,fmin,flag,outpt,allmins] = run(ms,problem,30);
%     disp(xmin)
%     estimates(index,2) = xmin(1);
%     estimates(index,3) = xmin(2);
%     index = index+1;
% end


% estimates
% temp4 = @(pars) (Hky - solve_model(X,Y,pars));
% problem = createOptimProblem('lsqnonlin','x0',pguess,'objective',temp4,'lb',[3,0.2],'ub',[100,1]);
temp5 = @(pars) est_all(castes,pars);
problem = createOptimProblem('fmincon','x0',pguess,'objective',temp5,'lb',[3,0.2],'ub',[100,1]);
ms = MultiStart('PlotFcns',@gsplotbestf,'UseParallel',true);
[xmin,fmin,flag,outpt,allmins] = run(ms,problem,40)
delete(gcp('nocreate'))
save('result.mat');
