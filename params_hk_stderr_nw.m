clc;

nsamples = 100;
Hks = zeros(nsamples,8);
params = zeros(nsamples, 2);
i=1;

% for i=1:nsamples
%     load(strcat('./bs_nw/result',num2str(i)),'xmin');
%     params(i,:) = xmin;
%     disp(xmin)
% end
% 
% ests = mean(params);
% se = std(params);

for i=1:nsamples
    load(strcat('./bs_nw/hk_result',num2str(i)),'H');
    Hks(i,:) = H;
end

pred_csr = mean(100*((1+Hks)./(1-Hks)))';
se_csr = std(100*((1+Hks)./(1-Hks)))';

pred_PG = mean((1-Hks)./2)'; %with substitution
se_PG = std((1-Hks)./2)';

