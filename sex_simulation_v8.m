clc;
clearvars;
figurepath = strcat('./Figures/');
jlim=1;
l=0;
tic
for j=1:jlim
    
    pop = 15055.5;
    nclass = 100;
    n = 100; %Number of income classes
    X = zeros(n,1);
    xin = 7.45 - 0*(0.01*l);
    xl = 8.5 + 0*(0.01*l);
    m = 7.639474 + 0*0.02*l; %mean
    vr = 0.0102265 + 0*(0.03*l); %variance
    %m = 8;
    %vr = 0.27;
    mu = log(m) - (0.5*log((vr/(m*m))+1));
    sig = sqrt(log((vr/(m*m))+1));
    Xdist = makedist('Normal','mu',m,'sigma',vr);
    %Xdist = makedist('Lognormal','mu',mu,'sigma',sig);
    %Xdist = makedist('Uniform','lower',xin,'upper',xl);
    X(:,1) = linspace(xin,xl,n); %Income classes
    X(:,2) = pdf(Xdist,X(:,1));
    %cdf_dom = cdf(Xdist,xl)-cdf(Xdist,xin);
    %X(:,2) = pdf(Xdist,X(:,1))/cdf_dom; %No. of boys/girls in each income class
    %delx = (xl-xin)/(n-1);
    %delx = 1;
    %X(:,2) = 15055.5*(X(:,2)/sum(delx*X(:,2)));
    X(:,2) = pop*(X(:,2)/sum(X(:,2)));
    %del(j)=delx;
    
    %[X, Hky_data, n] = getdata(1);
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
    
    a = 5;
    al = 0.65 + 0.2*l;
    
    params = [a,al];

    
    [Hky,ks,CYy,cYy,CXx,cXx,ky,phi,vxx,i,dy,S] = solve_model(X,Y,params);
    %PG(:,j) = (1-Hky)./(2-Hky); %without substitution
    PG(:,j) = (1-Hky)./2; %with substitution
    %PG_data(:,j) = (1-Hky_data)./2;
    H(:,j) = Hky;
    x(:,j) = X(:,2);
    domx(:,j) = X(:,1);
    domy(:,j) = Y(:,1);
    d(:,j) = dy(:,2);
    alph(j) = al;
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
legend(lvar,'Location','southeast')
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
spy(kvals(:,:,1),'r')
%spy(kvals(:,:,2),'g')
%spy(kvals(:,:,3),'b')
grid on
set(gca,'YDir','normal', 'GridLineStyle', ':')
refline(1,0)
xlabel('$x$','FontSize',14)
ylabel('$y$','FontSize',14)
%set(gca,'XTickLabel',labmy)
%set(gca,'YTickLabel',labmy)
legend(lvar,'Location','southeast')
title('Matching','FontSize',14)
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
hold on
[ax,h1,h2]=plotyy(domy(:,j),d(:,j),domy(:,j),PG(:,j));
xlabel('$y$','FontSize',14)
set(get(ax(1),'Ylabel'),'String','$d$','FontSize',14)
set(get(ax(2),'Ylabel'),'String','Pr(Girl)','FontSize',14)
% ylabel(gca(1),'$d$','FontSize',14)
% ylabel(gca(2),'Pr(Girl)','FontSize',14)
title(strcat('Dowry and Pr(Girl), $\alpha = ',num2str(al),'$'),'FontSize',14)
%legend('Dowry','Pr(Girl)','Location','southeast')
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

close all
%%
%{

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
figure(12)
set(figure(12),'defaulttextinterpreter','latex');
hold on
xlabel('$y$','FontSize',14)
ylabel('Pr(Girl)','FontSize',14)
title(strcat('Counterfactual Simulations, $\alpha = ',num2str(al),'$'),'FontSize',14)
plot(linspace(1,n,n),PG(:,j))
%plot(domy(:,j),PG(:,j))

[X, Hky_data, n] = getdata(28);
%X = equal_sized(X,nclass,pop);
%n = nclass;
%X = [X(:,1)*100 X(:,2)]
%X(:,1) = X(:,1)*10;
Y = X;
j=1;
tr = 0.2 * ((X(1,1)+X(2,1))/2);
% [Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S] = solve_model(X,Y,params);
% PG(:,j) = (1-Hky)./2; %with substitution
% H(:,j) = Hky;
% x(:,j) = X(:,2);
% domx(:,j) = X(:,1);
% domy(:,j) = Y(:,1);
% d(:,j) = dy(:,2);
% var(j) = vr;
% plot(linspace(1,n),PG(:,j),'color','r')
% plot(domy(:,j),PG(:,j),'color','r')
[Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S] = solve_model_top(X,Y,params,tr);
PG(:,j) = (1-Hky)./2; %with substitution
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
plot(linspace(1,n,n),PG(:,j),'color','g')
%plot(domy(:,j),PG(:,j),'color','g')

tr = tr/4;
[Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S] = solve_model_all(X,Y,params,tr);
PG(:,j) = (1-Hky)./2; %with substitution
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
plot(linspace(1,n,n),PG(:,j),'color','r')

[Hky,ks,CYy,cYy,CXx,cXx,ky,vyy,vxx,i,dy,S] = solve_model_girls(X,Y,params,tr);
PG(:,j) = (1-Hky)./2; %with substitution
H(:,j) = Hky;
x(:,j) = X(:,2);
domx(:,j) = X(:,1);
domy(:,j) = Y(:,1);
d(:,j) = dy(:,2);
%var(j) = vr;
plot(linspace(1,n,n),PG(:,j),'color','black')

legend('Benchmark','Tr to fathers of poor girls','Tr to fathers of all girls','Tr to all girls','Location','northeast')
print('-dpdf', strcat(figurepath, 'counter1.pdf'));
hold off
close all
%}