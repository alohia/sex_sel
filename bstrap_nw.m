clc;

sample_num = castes(end);
castes = castes(1:end-1);
load(strcat('./bs_nw/result',num2str(sample_num)),'xmin');
a = xmin(1);
al = xmin(2);
params = [a,al];

num_castes = numel(castes);
Hk_pred = zeros(num_castes,8);
for i = 1:num_castes
    [X, Hky ,n] = getdata(castes(i));
    Y = X;
    Hk_pred(i,:) = (solve_model(X,Y,params))';
end
H = mean(Hk_pred);
save(strcat('./bs_nw/hk_result',num2str(sample_num),'.mat'));