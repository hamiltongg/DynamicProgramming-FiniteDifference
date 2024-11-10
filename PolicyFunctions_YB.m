%{
Application: Dynamic Programming approach to solve 
two-agent heterogeneous economy: Wang(1996)'s paper
-------------------------------------------------
I solve the PDE of Ak using a function: A1, A2
Then, I find "S" based on equilibrium conditions

----------------------------
Author: Hamilton Galindo Gil
Date:   2023 (July), 2024 (Nov)
Paper base: Wang(1996) 
----------------------------
%}
%=========================================

%% Solution of the Model
close all;
clear all;
clc;
gamma1 = 1.5;
gamma2 = 0.5;

gammakg1 = gamma1;
gammakg2 = gamma2;

mu = 0.05;
sigma = 0.3;
lambda = 1/3;

[c1g1,c2g1,Yg1,A1,psig1,rg1, N1, N2] = Ak_YB_fun(gamma1,gamma2,gammakg1,lambda, mu, sigma); % Sol of Agent 1:gammakg1 
[c1g2,c2g2,Yg2,A2,psig2,rg2, N1, N2] = Ak_YB_fun(gamma1,gamma2,gammakg2,lambda, mu, sigma); % Sol of Agent 2:gammakg2

%% Variables in Equilibrium
% State var
Y = Yg1';
r = rg1';
psi = psig1;

% Risky asset price: S
S = c1g1'.*A1.^(1/gamma1) + c2g1'.*A2.^(1/gamma2);

% Volatility of S: sigmat
deltaY = Y(2)-Y(1);
Sy =  (S(2:end) - S(1:end-1))/deltaY;
    sigma0 = sigma*(Y(2:end)./S(2:end)).*Sy; %500 points
sigmat = [0;sigma0]; %501 points

% Price-Dividend Ratio
pd = S./Y;

%Expected Rate of Return
beta = psi.*sigmat + r;

%Stochastic Discount Factor
m = lambda*c1g1.^(-gamma1);

% Wealth
W1 = c1g1'.*A1.^(1/gamma1); % agent1
W2 = c2g1'.*A2.^(1/gamma2); % agent2

% Optimal Portfolio
dA1Y = (A1(2:end) - A1(1:end-1))/deltaY; % starts at i=2 (500 points)
dA2Y = (A2(2:end) - A2(1:end-1))/deltaY; % starts at i=2 (500 points)

    %risky asset
    w110 = ( (sigma*Y(2:end))./sigmat(2:end) ).*dA1Y./(gamma1*A1(2:end)) + psi(2:end)'./(gamma1*sigmat(2:end));
    w210 = ( (sigma*Y(2:end))./sigmat(2:end) ).*dA2Y./(gamma2*A2(2:end)) + psi(2:end)'./(gamma2*sigmat(2:end));

    %riskless asset
    w120 = 1- w110;
    w220 = 1 -w210;
    
    w11 = [0; w110]; %agent1: risky asset
    w21 = [0; w210]; %agent2: risky asset
    w12 = [0; w120]; %agent1: riskless asset
    w22 = [0; w220]; %agent2: riskless asset

% Risky Asset Shares
N11 = w11.*W1./S;
N21 = w21.*W2./S;

% Money invested in riskless asset
NB1 = w12.*W1;
NB2 = w22.*W2;

%% Graphs 2 (policy functions)
tinit = 2;
tend = 101;

% Plot (Fig 1: Optimal consumption)
figname1 = strcat('Fig1: Optimal consumption', ' (\lambda=',num2str(lambda),')');

figure('Name',figname1)
subplot(2,2,1)
    plot(Y(tinit:tend),c1g1(tinit:tend),'r:',...
         Y(tinit:tend),c2g1(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    %titlestr = strcat('Endowment and Optimal Consumption','($\lambda$=',num2str(lambda),')');
    %title(titlestr,'interpreter','latex')
    title('Optimal Consumption (c_k)')
    %l1 = strcat('c_1','(\gamma_1=',num2str(gamma1),')');
    %l2 = strcat('c_2','(\gamma_2=',num2str(gamma2),')');
    %legend(l1, l2,Location='best')
    legend('c_1', 'c_2',Location='best')
    grid;

subplot(2,2,2)
    plot(Y(tinit:tend),S(tinit:tend),'k',...
         Y(tinit:tend),W1(tinit:tend),'r:',...
         Y(tinit:tend),W2(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Stock Price (S) and Wealth (W_k)')
    legend('Stock Price (S)','Wealth of Agent 1 (W_1)','Wealth of Agent 2 (W_2)',Location='best')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P2','_Fig1.pdf'));    

%--------------------------------------
% Plot (Fig 2: Optimal Portfolio)
figname2 = strcat('Fig2: Optimal Portfolio', ' (\lambda=',num2str(lambda),')');

figure('Name',figname2)
subplot(2,2,1)
    plot(Y(tinit:tend), w11(tinit:tend),'r:',......
         Y(tinit:tend), w21(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Optimal Portfolio: risky asset (\omega_k^{(1)})')
    legend('\omega_1^{(1)}','\omega_2^{(1)}',Location='best')
    grid;

subplot(2,2,2)
    plot(Y(tinit:tend), w12(tinit:tend),'r:',......
         Y(tinit:tend), w22(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Optimal Portfolio: riskless asset (\omega_k^{(2)})')
    legend('\omega_1^{(2)}','\omega_2^{(2)}',Location='best')
    grid;

subplot(2,2,3)
    plot(Y(tinit:tend),N11(tinit:tend),'r:',...
         Y(tinit:tend),N21(tinit:tend),'b--',...
        'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    %titlestr1 = strcat('Riksy Asset Shares', '($\lambda$=',num2str(lambda),')');
    %title(titlestr1,'Interpreter','latex')
    title('Risky Asset Shares (N_k^{(1)})')
    legend('N_1^{(1)}', 'N_2^{(1)}',Location='best')
    grid;
    
subplot(2,2,4)
    plot(Y(tinit:tend),NB1(tinit:tend),'r:',...
         Y(tinit:tend),NB2(tinit:tend),'b--',...
        'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Money invested in riskless asset (B*N_k^{(2)})')
    legend('B*N_1^{(2)}', 'B*N_2^{(2)}',Location='best')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P2','_Fig2.pdf'));    

%--------------------------------------
% Plot (Fig 3: Asset Prices)
figname3 = strcat('Fig2: Asset Prices', ' (\lambda=',num2str(lambda),')');

figure('Name',figname3)
subplot(2,2,1)
    plot(Y(tinit:tend), pd(tinit:tend),'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Price-Dividend Ratio (S/Y)')
    grid;
    
subplot(2,2,2)
    plot(Y(tinit:tend),r(tinit:tend),'LineWidth',1.5)
    title('Interest Rate (r)')
    xlabel('Endowment (Y)')
    grid;

subplot(2,2,3)
    plot(Y(tinit:tend),psi(tinit:tend),'r',...
         Y(tinit:tend),beta(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Prices I')
    legend('Price of Risk (\psi)', 'Expected Rate of Return (\beta)',Location='best')
    grid;

subplot(2,2,4)
    plot(Y(tinit:tend),m(tinit:tend),'r',...
         Y(tinit:tend),sigmat(tinit:tend),'k--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Prices II')
    legend('Stochastic Discount Factor (m)', 'Stock Volatility (\nu)',Location='best')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P2','_Fig3.pdf'));    

%% Consumption-Wealth ratio

cw1 = c1g1'./W1;
cw2 = c2g1'./W2;

figure('Name','C-W ratio')
plot(Y,cw1,'--',Y,cw2,...
     Y, N1*ones(501,1),'k:',...
     Y, N2*ones(501,1),'k--','LineWidth', 1.5)
    l1 = strcat('c_1/W_1','(RRA=',num2str(gamma1),')');
    l2 = strcat('c_2/W_2','(RRA=',num2str(gamma2),')');
    legend(l1, l2,'Only Agent 1', 'Only Agent 2')
%legend('Agent 1 (RRA=0.8)', 'Agent 2 (RRA=0.5)', 'Only Agent 1', 'Only Agent 2')
title('Consumption-Wealth Ratio (c/W)')
grid

%% Consumption-share

cs1 = c1g1'./Y;
cs2 = c2g1'./Y;

figure('Name','C-W ratio')
plot(Y(tinit:tend),cs1(tinit:tend),'--',...
     Y(tinit:tend),cs2(tinit:tend),'LineWidth', 1.5)
    l1 = strcat('c_1/Y','(RRA=',num2str(gamma1),')');
    l2 = strcat('c_2/Y','(RRA=',num2str(gamma2),')');
    legend(l1, l2)
%legend('Agent 1 (RRA=0.8)', 'Agent 2 (RRA=0.5)', 'Only Agent 1', 'Only Agent 2')
title('Consumption-share(c/Y)')
grid

%% Stats
corr_sigmat_Y = corrcoef(sigmat(2:501), Y(2:501)); %-0.26