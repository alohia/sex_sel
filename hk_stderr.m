clc;
clearvars;

nsamples = 100;
Hks = zeros(nsamples-1,8);
i=1;

for i=1:nsamples
    if i==28
        continue;
    end
    load(strcat('./bs_results/hk_result',num2str(i)),'H');
    Hks(i,:) = H;
end

pred_csr = mean(100*((1+Hks)./(1-Hks)))';
se_csr = std(100*((1+Hks)./(1-Hks)))';