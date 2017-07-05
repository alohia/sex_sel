clc;
clearvars;
figurepath = strcat('./Figures/');
jlim=1;
l=0;
tic
for j=1:jlim
    
    pop = 72943;
    nclass = 10;
    n = 100; %Number of income classes
    X = zeros(n,1);
    xin = 0.5595853 - 0*(0.01*l) + 0.6;
    xl = 3.488218 + 0*(0.01*l);
    m = 1.971046 + 0*0.02*l; %mean
    vr = 0.4205418 + 0*(0.03*l); %variance
    %m = 8;
    %vr = 0.27;
    mu = log(m) - (0.5*log((vr/(m*m))+1));
    sig = sqrt(log((vr/(m*m))+1));
    %Xdist = makedist('Normal','mu',m,'sigma',vr);
    %Xdist = makedist('Lognormal','mu',mu,'sigma',sig);
    Xdist = makedist('Uniform','lower',xin,'upper',xl);
    X(:,1) = linspace(xin,xl,n); %Income classes
    X(:,2) = pdf(Xdist,X(:,1));
    %cdf_dom = cdf(Xdist,xl)-cdf(Xdist,xin);
    %X(:,2) = pdf(Xdist,X(:,1))/cdf_dom; %No. of boys/girls in each income class
    %delx = (xl-xin)/(n-1);
    %delx = 1;
    %X(:,2) = 15055.5*(X(:,2)/sum(delx*X(:,2)));
    X(:,2) = pop*(X(:,2)/sum(X(:,2)));
    %del(j)=delx;
    
    %[X, Hky_data, n] = getdata(28);
    %X = equal_sized(X,nclass,pop);
    %n = nclass;
    %X = [X(:,1)*100 X(:,2)]
    %X(:,1) = X(:,1)*10;
    Y = X;
    
    % Parameters for caste 28
    a = 12.6887;
    al = 0.6121;
    
    % Final estimates
    a = 16.4475;
    al = 0.6090;
    
    % Estimates with weights
    
    a = 12.0153;
    al = 0.6111;
    
    % Parameters for caste 1
    a = 10.1179; %upper limit for H distribution
    al = 0.6043; %exogenous alpha - husband's share
    
    a = 12 + 0*8*l;
    al = 0.6 - 0*0.2*l;
    
    %a = 12.9751;
    %al = 0.6128;
    params = [a,al];

    
    [Hky,ks,CYy,cYy,CXx,cXx,ky,phi,vxx,i,dy,S,dgiven] = solve_model(X,Y,params);
    %PG(:,j) = (1-Hky)./(2-Hky); %without substitution
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = (1+Hky)./(1-Hky);
    %PG_data(:,j) = (1-Hky_data)./2;
    H(:,j) = Hky;
    x(:,j) = X(:,2);
    domx(:,j) = X(:,1);
    domy(:,j) = Y(:,1);
    d(:,j) = dy(:,2);
    d(d(:,j)==0,j) = NaN;
    dgiv(:,j) = dgiven(:,2);
    dgiv(dgiv(:,j)==0,j) = NaN;
    alph(j) = al;
    as(j) = a;
    kvals(:,:,j) = ks(:,:,end);
    %var(j) = vr;
    l=l+1;
end
toc
%%
%X-Tick labels for matching plot go in labmy
labmy = {};
r=1;
newv = linspace(xin,xl,(10)+1);
while (r<=11)
    labmy{r} = num2str(newv(r));
    r = r + 1;
end
%%
figure(1)
set(figure(1),'defaulttextinterpreter','latex');
plot(CYy(:,1),CYy(:,2),cYy(:,1),cYy(:,2),CXx(:,1),CXx(:,2),cXx(:,1),cXx(:,2))
xlabel('$x,y$','FontSize',14)
ylabel('$C_x, c_x, C_y, c_y$','FontSize',14)
title(strcat('Consumptions, $\alpha = ',num2str(al),'$'),'FontSize',14)
legend('C_y','c_y','C_x','c_x','Location','southeast')
print('-dpdf', strcat(figurepath, 'Consumptions.pdf'));
hold off
close all

figure(2)
set(figure(2),'defaulttextinterpreter','latex');
scatter(ky(:,1),ky(:,2));
%scatter(linspace(1,n),ky(:,2));
xlabel('$y$','FontSize',14)
ylabel('$k*$','FontSize',14)
title(strcat('k* vs y, $\alpha = ',num2str(al),'$'),'FontSize',14)
print('-dpdf', strcat(figurepath, 'kstar.pdf'));
hold off
close

figure(3)
set(figure(3),'defaulttextinterpreter','latex');
%gin='';
%slp='';
lvar={};
col=hsv(jlim);
hold on
for j=1:jlim
    %s=polyfit(domy(:,j),H(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),H(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),H(:,j),'color',col(j,:))
    %gin = strcat(gin,'gini',num2str(j),'=',num2str(ginicoeff(x(:,j),domx(:,j))),', ');
    %lvar{j} = strcat('x',num2str(j),':var = ',num2str(var(j)));
    lvar{j} = strcat('\alpha = ',num2str(alph(j)));
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
%xlabel(strcat('$y$,',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('$H(k^{\star})$','FontSize',14)
legend(lvar,'Location','southeast')
%title(strcat('H(k*) vs y, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('H(k*) vs y, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('H(k*) vs Income-class rank','FontSize',14)
print('-dpdf', strcat(figurepath, 'H(kstar).pdf'));
hold off
close

figure(4)
set(figure(4),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),PG(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),PG(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),PG(:,j),'color',col(j,:))
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% scatter(linspace(1,n,n),PG_data(:,j))
% lsline
%xlabel(strcat('$y$,',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('Pr(Girl)','FontSize',14)
legend(lvar,'Location','northeast')
%title(strcat('Pr(Girl) vs y, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Pr(Girl) vs y, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('Pr(Girl) vs Income-class rank','FontSize',14)
print('-dpdf', strcat(figurepath, 'Pr(Girl).pdf'));
hold off
close

figure(5)
set(figure(5),'defaulttextinterpreter','latex');
plot(domx(:,jlim),S,domx(:,jlim),vxx)
xlabel('$x$','FontSize',14)
ylabel('$v(x)$','FontSize',14)
title(strcat('v(x), $\alpha = ',num2str(al),'$, Iterations = ',num2str(i)),'FontSize',14)
legend('S','v1','v2','v3','v4','Location','southeast')
print('-dpdf', strcat(figurepath, 'v(x).pdf'));
close

figure(6)
set(figure(6),'defaulttextinterpreter','latex');
hold on
kvl = (kvals(:,:,1)>0)';
[kvlx, kvly] = find(kvl);
plot(kvlx,kvly,'.-r','LineWidth',0.5,'MarkerSize',2)
xlim([0 n])
ylim([0 n])
%spy(kvals(:,:,1),'r',7)
%spy(kvals(:,:,2),'g')
%spy(kvals(:,:,3),'b')
grid on
set(gca,'YDir','normal', 'GridLineStyle', ':')
refline(1,0)
xlabel('male wealth class','FontSize',14)
ylabel('female wealth class','FontSize',14)
%set(gca,'XTickLabel',labmy)
%set(gca,'YTickLabel',labmy)
%legend(lvar,'Location','southeast')
%title('Matching','FontSize',14)
print('-dpdf', strcat(figurepath, 'Match.pdf'));
close

figure(7)
set(figure(7),'defaulttextinterpreter','latex');
plot(phi(:,1),S,phi(:,1),phi(:,2))
xlabel('$y$','FontSize',14)
ylabel('$\Phi(y)$','FontSize',14)
title(strcat('$\Phi(y), \alpha = ',num2str(al),'$'),'FontSize',14)
legend('S','\Phi(y)','Location','southeast')
print('-dpdf', strcat(figurepath, 'Phi(y).pdf'));
close

figure(8)
%lvar={};
set(figure(8),'defaulttextinterpreter','latex');
hold on
for j=1:jlim
    %lvar{j} = strcat('x',num2str(j),':var = ',num2str(var(j)));
    plot(domx(:,j),x(:,j),'color',col(j,:))
end
xlabel('$x$','FontSize',14)
ylabel('$No. of people$','FontSize',14)
%title(strcat('$X, Total = ',num2str(delx*sum(X(:,2))),'$'),'FontSize',14)
title(strcat('$X, Total = ',num2str(sum(X(:,2))),'$'),'FontSize',14)
%legend(lvar,'Location','southeast')
print('-dpdf', strcat(figurepath, 'X.pdf'));
hold off
close

figure(9)
set(figure(9),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),d(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),d(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),d(:,j),'color',col(j,:))
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% xlabel(strcat('$y$, ',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('$d$','FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('Dowry demanded vs Income-class rank','FontSize',14)
legend(lvar,'Location','southeast')
print('-dpdf', strcat(figurepath, 'd.pdf'));
hold off
close

figure(10)
j=1;
set(figure(10),'defaulttextinterpreter','latex');
%title(strcat('Dowry/Wealth and Pr(Girl)',' | $\alpha=$',num2str(al)),'FontSize',14)
xlabel('wealth class','FontSize',14)
lvar2={};
hold on
for j=1:jlim
    yyaxis right
    plot(linspace(1,n,n),d(:,j)./X(:,1))
    yyaxis left
    plot(linspace(1,n,n),PG(:,j))
    lvar2{j} = strcat('\alpha = ',num2str(alph(j)));
end
yyaxis right
ylabel('dowry/wealth','FontSize',14)
yyaxis left
ylabel('Pr(Girl)','FontSize',14)
% set(get(ax(1),'Ylabel'),'String','$d$','FontSize',14)
% set(get(ax(2),'Ylabel'),'String','Pr(Girl)','FontSize',14)
% ylabel(gca(1),'$d$','FontSize',14)
% ylabel(gca(2),'Pr(Girl)','FontSize',14)
legend(lvar2,'Location','southwest')
legend('boxoff')
print('-dpdf', strcat(figurepath, 'd+Pr.pdf'));
hold off
close

figure(11)
set(figure(11),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),d(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),d(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),(d(:,j)./X(:,1)),'color',col(j,:))
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% xlabel(strcat('$y$, ',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('$d/x$','FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('Dowry/Income vs Income-class rank','FontSize',14)
legend(lvar,'Location','southeast')
print('-dpdf', strcat(figurepath, 'dstat.pdf'));
hold off
close

figure(12)
set(figure(12),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),PG(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),PG(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),csr(:,j),'color',col(j,:))
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% scatter(linspace(1,n,n),PG_data(:,j))
% lsline
%xlabel(strcat('$y$,',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('Child sex ratio','FontSize',14)
legend(lvar,'Location','southeast')
%title(strcat('Pr(Girl) vs y, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Pr(Girl) vs y, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('CSR vs Income-class rank','FontSize',14)
print('-dpdf', strcat(figurepath, 'CSR.pdf'));
hold off
close

figure(13)
set(figure(13),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),d(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),d(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),(d(:,j)./X(:,1)),'color',col(j,:))
    plot(linspace(1,n,n),(dgiv(:,j)./Y(:,1)),'color','blue')
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% xlabel(strcat('$y$, ',gin),'FontSize',14)
xlabel('wealth class','FontSize',14)
ylabel('dowry/wealth','FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),'$'),'FontSize',14)
%title('Dowry/Income vs Income-class rank','FontSize',14)
legend('Received','Given','Location','southeast')
print('-dpdf', strcat(figurepath, 'dstat_both.pdf'));
hold off
close

figure(14)
set(figure(14),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),d(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),d(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),(d(:,j)),'color',col(j,:))
    plot(linspace(1,n,n),(dgiv(:,j)),'color','blue')
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% xlabel(strcat('$y$, ',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('$d$','FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('Dowry vs Income-class rank','FontSize',14)
legend('Received','Given','Location','southeast')
print('-dpdf', strcat(figurepath, 'd_both.pdf'));
hold off
close

figure(15)
j=1;
set(figure(15),'defaulttextinterpreter','latex');
%title(strcat('Dowry/Wealth and Pr(Girl)',' | $\alpha=$',num2str(al)),'FontSize',14)
xlabel('wealth class','FontSize',14)
hold on
for j=1:jlim
    yyaxis left
    p1 = plot(linspace(1,n,n),PG(:,j));
    yyaxis right
    p2 = plot(linspace(1,n,n),(d(:,j)./X(:,1)));
    p3 = plot(linspace(1,n,n),(dgiv(:,j)./Y(:,1)));
end
yyaxis left
ylabel('Pr(Girl)','FontSize',14)
yyaxis right
ylabel('dowry/wealth','FontSize',14)
% set(get(ax(1),'Ylabel'),'String','$d$','FontSize',14)
% set(get(ax(2),'Ylabel'),'String','Pr(Girl)','FontSize',14)
% ylabel(gca(1),'$d$','FontSize',14)
% ylabel(gca(2),'Pr(Girl)','FontSize',14)
legend([p2 p3],{'Received','Given'},'Location','southwest')
legend('boxoff')
print('-dpdf', strcat(figurepath, 'dstat+Pr.pdf'));
hold off
close


close all
%%


%{
xin = 7.45 + 0.2;
xl = 8 + 0.2;
m = 7.639474 + 0.2; %mean
vr = 0.0102265 + 0; %variance
mu = log(m) - (0.5*log((vr/(m*m))+1));
sig = sqrt(log((vr/(m*m))+1));
%Xdist = makedist('Normal','mu',m,'sigma',vr);
Ydist = makedist('Lognormal','mu',mu,'sigma',sig);
%Xdist = makedist('Uniform','lower',xin,'upper',xl);
Y(:,1) = linspace(xin,xl,n); %Income classes
Y(:,2) = pdf(Ydist,Y(:,1)); %No. of boys/girls in each income class
Y(:,2) = 15055.5*(Y(:,2)/sum(Y(:,2)));
%}

castes = [1 .0984789; 4 0.0090223; 5 0.018784; 7 0.0151011; 10 .0221519; 14 .0589121; 17 .0112403; 23 .0183048; 28 .4044578; 33 .0585972; 38 .2693966; 44 .0155529];
[X, Hky_data, n] = getdata(28);

figure(16)
set(figure(16),'defaulttextinterpreter','latex');
hold on
xlabel('wealth class','FontSize',14)
ylabel('child sex ratio','FontSize',14)
%title(strcat('Counterfactual Simulations, $\alpha = ',num2str(al),'$'),'FontSize',14)
csr = zeros(length(X),1);
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    [Hky,ks,CYy,cYy,CXx,cXx,ky,phi,vxx,~,dy,S,dgiven] = solve_model(X,Y,params);
    %PG(:,j) = (1-Hky)./(2-Hky); %without substitution
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
end
%plot(linspace(1,n,n),100*csr(:,j),'color','b', 'LineWidth', 1)
plot(linspace(1,n,n),100*csr(:,j),'color','b')
%plot(domy(:,j),PG(:,j))


%[X, Hky_data, n] = getdata(28);
%X = equal_sized(X,nclass,pop);
%n = nclass;
%X = [X(:,1)*100 X(:,2)]
%X(:,1) = X(:,1)*10;
csr = zeros(length(X),1);
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    tr = 0.2 * ((X(1,1)+X(2,1))/2);
    [Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,~,dy,S] = solve_model_top(X,Y,params,tr);
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
    csr
end
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
%plot(linspace(1,n,n),100*csr(:,j),'color','r', 'LineWidth', 1.5)
plot(linspace(1,n,n),100*csr(:,j),'color','r')
%plot(domy(:,j),PG(:,j),'color','g')


csr = zeros(length(X),1);
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    tr = (0.2 * ((X(1,1)+X(2,1))/2))/4;
    [Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,~,dy,S] = solve_model_all(X,Y,params,tr);
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
    csr
end
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
%plot(linspace(1,n,n),100*csr(:,j),'color','r','LineStyle','--', 'LineWidth', 1)
plot(linspace(1,n,n),100*csr(:,j),'color','r','LineStyle','--')

csr = zeros(length(X),1);
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    tr = (0.2 * ((X(1,1)+X(2,1))/2))/4;
    [Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,~,dy,S] = solve_model_girls(X,Y,params,tr);
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
    csr
end
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
%plot(linspace(1,n,n),100*csr(:,j),'color',[0 .5 0])
%plot(linspace(1,n,n),100*csr(:,j),'color','k', 'LineStyle', '--', 'LineWidth', 1)
plot(linspace(1,n,n),100*csr(:,j),'color','k', 'LineStyle', '--')
%plot(linspace(1,n,n),100*csr(:,j),'color','k', 'LineWidth', 1.5)


legend('Benchmark','Transfer to parents of poor girls','Transfer to parents of all girls','Transfer to all girls','Location','southeast')
print('-dpdf', strcat(figurepath, 'counter_black_normal.pdf'));
hold off
close all
%%

[X, Hky_data, n] = getdata(28);

figure(17)
set(figure(17),'defaulttextinterpreter','latex');
hold on
xlabel('wealth class','FontSize',14)
ylabel('child sex ratio','FontSize',14)
%title(strcat('Counterfactual Simulations, $\alpha = ',num2str(al),'$'),'FontSize',14)

csr = zeros(length(X),1);
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    [Hky,ks,CYy,cYy,CXx,cXx,ky,phi,vxx,~,dy,S,dgiven] = solve_model(X,Y,params);
    %PG(:,j) = (1-Hky)./(2-Hky); %without substitution
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
end
plot(linspace(1,n,n),100*csr(:,j),'color','b', 'LineWidth', 1.5)
lvar={};
lvar{1} = 'Benchmark';
%plot(domy(:,j),PG(:,j))

csr = zeros(length(X),1);
theta = 0.9;
for i = 1:length(castes)
    [X, Hky_data, n] = getdata(castes(i));
    Y = X;
    j = 1;
    [Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,~,dy,S,dgiven] = solve_model_dtax(X,Y,params,theta);
    PG(:,j) = (1-Hky)./2; %with substitution
    csr(:,j) = csr(:,j) + castes(i,2).*((1+Hky)./(1-Hky));
end
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
dgiv(:,j) = dgiven(:,2);
dgiv(dgiv(:,j)==0,j) = NaN;
%var(j) = vr;
plot(linspace(1,n,n),100*csr(:,j),'color','r','LineStyle','--', 'LineWidth', 1.5)
lvar{2}=strcat('\theta=',num2str(theta));
%plot(domy(:,j),PG(:,j),'color','g')


[X, Hky_data, n] = getdata(28);
Y = X;
j = 1;
theta = 1.1;
[Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S] = solve_model_dtax(X,Y,params,theta);
PG(:,j) = (1-Hky)./2; %with substitution
csr(:,j) = (1+Hky)./(1-Hky);
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
plot(linspace(1,n,n),100*csr(:,j),'color','r')
lvar{3} = strcat('\theta=',num2str(theta));

legend(lvar,'Location','southeast')
print('-dpdf', strcat(figurepath, 'counter_tax.pdf'));
hold off
close all

figure(18)
set(figure(18),'defaulttextinterpreter','latex');
hold on
%slp='';
for j=1:jlim
    %s=polyfit(domy(:,j),d(:,j),1);
    %s=s(1)*100;
    %plot(domy(:,j),d(:,j),'color',rand(1,3))
    plot(linspace(1,n,n),(d(:,j)./X(:,1)),'color',col(j,:))
    plot(linspace(1,n,n),(dgiv(:,j)./Y(:,1)),'color','blue')
    %slp = strcat(slp,'slope',num2str(j),'=',num2str(s),', ');
end
% xlabel(strcat('$y$, ',gin),'FontSize',14)
xlabel('$income-class rank$','FontSize',14)
ylabel('$d/x$','FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),', ',slp,'$'),'FontSize',14)
%title(strcat('Dowry, $\alpha = ',num2str(al),'$'),'FontSize',14)
title('Dowry/Income vs Income-class rank','FontSize',14)
legend(lvar,'Location','southeast')
%print('-dpdf', strcat(figurepath, 'dstat.pdf'));
hold off
close