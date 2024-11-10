%{
Inspecting Risk-Aversion Heterogeneity:
    - Increasing \gamma_1
----------------------------
Author: Hamilton Galindo Gil
Date:   2024 (Nov)
Paper base: Wang(1996) 
----------------------------
%}

close all;
clear all;
clc;

% Same parameters for both simulations
mu = 0.05;
sigma = 0.3;
lambda = 1/3;

%--------------
% Simulation 1
%--------------
gamma1 = 0.8;
gamma2 = 0.5;

gammakg1 = gamma1;
gammakg2 = gamma2;

[c1g1,c2g1,Yg1,A1,psig1,rg1, N1, N2] = Ak_YB_fun(gamma1,gamma2,gammakg1,lambda, mu, sigma); % Sol of Agent 1:gammakg1 
[c1g2,c2g2,Yg2,A2,psig2,rg2, N1, N2] = Ak_YB_fun(gamma1,gamma2,gammakg2,lambda, mu, sigma); % Sol of Agent 2:gammakg2

%--------------
% Simulation 2
%--------------
gamma10 = 1.6;
gamma20 = 0.5;

gammakg10 = gamma10;
gammakg20 = gamma20;

[c1g10,c2g10,Yg10,A10,psig10,rg10, N10, N20] = Ak_YB_fun(gamma10,gamma20,gammakg20,lambda, mu, sigma); % Sol of Agent 2:gammakg2
[c1g20,c2g20,Yg20,A20,psig20,rg20, N10, N20] = Ak_YB_fun(gamma1,gamma2,gammakg2,lambda, mu, sigma); % Sol of Agent 2:gammakg2

%--------------
% Variables (Simulation 1)
%--------------
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

% Wealth
W1 = c1g1'.*A1.^(1/gamma1); % agent1
W2 = c2g1'.*A2.^(1/gamma2); % agent2

% Optimal Portfolio
dA1Y = (A1(2:end) - A1(1:end-1))/deltaY; % starts at i=2 (500 points)
dA2Y = (A2(2:end) - A2(1:end-1))/deltaY; % starts at i=2 (500 points)

    %risky asset
    w110 = ( (sigma*Y(2:end))./sigmat(2:end) ).*dA1Y./(gamma1*A1(2:end)) + psi(2:end)'./(gamma1*sigmat(2:end));
    w210 = ( (sigma*Y(2:end))./sigmat(2:end) ).*dA2Y./(gamma2*A2(2:end)) + psi(2:end)'./(gamma2*sigmat(2:end));
    
    w11 = [0; w110]; %agent1: risky asset
    w21 = [0; w210]; %agent2: risky asset

%--------------
% Variables (Simulation 2)
%--------------
% State var
Y0 = Yg10';
r0 = rg10';
psi0 = psig10;

% Risky asset price: S
S0 = c1g10'.*A10.^(1/gamma10) + c2g10'.*A20.^(1/gamma20);

% Volatility of S: sigmat
deltaY0 = Y0(2)-Y0(1);
Sy0 =  (S0(2:end) - S0(1:end-1))/deltaY0;
    sigma00 = sigma*(Y0(2:end)./S0(2:end)).*Sy0; %500 points
sigmat0 = [0;sigma00]; %501 points

% Price-Dividend Ratio
pd0 = S0./Y0;

% Wealth
W10 = c1g10'.*A10.^(1/gamma10); % agent1
W20 = c2g10'.*A20.^(1/gamma20); % agent2

% Optimal Portfolio
dA1Y0 = (A10(2:end) - A10(1:end-1))/deltaY0; % starts at i=2 (500 points)
dA2Y0 = (A20(2:end) - A20(1:end-1))/deltaY0; % starts at i=2 (500 points)

    %risky asset
    w1100 = ( (sigma*Y0(2:end))./sigmat0(2:end) ).*dA1Y0./(gamma10*A10(2:end)) + psi0(2:end)'./(gamma10*sigmat0(2:end));
    w2100 = ( (sigma*Y0(2:end))./sigmat0(2:end) ).*dA2Y0./(gamma20*A20(2:end)) + psi0(2:end)'./(gamma20*sigmat0(2:end));
    
    w110 = [0; w1100]; %agent1: risky asset
    w210 = [0; w2100]; %agent2: risky asset

%--------------
%% Graphs
%--------------
tinit = 2;
tend = 101;

% Plot (Fig 1: consumption | wealth | risky portfolio)
figname1 = strcat('Fig1: Optimal consumption', ' (\lambda=',num2str(lambda),')');

figure('Name',figname1)
subplot(2,3,1)
    plot(Y(tinit:tend),c1g1(tinit:tend),'r',...
         Y(tinit:tend),c2g1(tinit:tend),'b',...
         Y0(tinit:tend),c1g10(tinit:tend),'r--',...
         Y0(tinit:tend),c2g10(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    %titlestr = strcat('Endowment and Optimal Consumption','($\lambda$=',num2str(lambda),')');
    %title(titlestr,'interpreter','latex')
    title('Optimal Consumption (c_k)')
    %l1 = strcat('c_1','(\gamma_1=',num2str(gamma1),')');
    %l2 = strcat('c_2','(\gamma_2=',num2str(gamma2),')');
    %legend(l1, l2,Location='best')
    legend('c_1 (\gamma_1=0.8)', 'c_2 (\gamma_2=0.5)',...
           'c_1 (\gamma_1=1.6)', 'c_2 (\gamma_2=0.5)',Location='best')
    grid;

subplot(2,3,2)
    plot(Y(tinit:tend),W1(tinit:tend),'r',...
         Y(tinit:tend),W2(tinit:tend),'b',...
         Y0(tinit:tend),W10(tinit:tend),'r--',...
         Y0(tinit:tend),W20(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Wealth (W_k)')
    legend('W_1 (\gamma_1=0.8)', 'W_2 (\gamma_2=0.5)',...
           'W_1 (\gamma_1=1.6)', 'W_2 (\gamma_2=0.5)',Location='best')
      grid;

subplot(2,3,3)
    plot(Y(tinit:tend), w11(tinit:tend),'r',...
         Y(tinit:tend), w21(tinit:tend),'b',...
         Y0(tinit:tend), w110(tinit:tend),'r--',...
         Y0(tinit:tend), w210(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Optimal Portfolio: risky asset (\omega_k^{(1)})')
    legend('\omega_1^{(1)} (\gamma_1=0.8)', '\omega_2^{(1)} (\gamma_2=0.5)',...
           '\omega_1^{(1)} (\gamma_1=1.6)', '\omega_2^{(1)} (\gamma_2=0.5)',Location='northeast')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P2','_RRA_Fig1.pdf'));    

%--------------------------------------
% Plot (Fig 2: PD | r | stock volatility)
figname3 = strcat('Fig2: Asset Prices', ' (\lambda=',num2str(lambda),')');

figure('Name',figname3)
subplot(2,3,1)
    plot(Y(tinit:tend), pd(tinit:tend),'b',...
         Y0(tinit:tend), pd0(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Price-Dividend Ratio (S/Y)')
    legend('S/Y (\gamma_1=0.8, \gamma_2=0.5)','S/Y (\gamma_1=1.6, \gamma_2=0.5)',Location='best')
    grid;
    
subplot(2,3,2)
    plot(Y(tinit:tend),r(tinit:tend), 'b', ...
         Y0(tinit:tend),r0(tinit:tend), 'b--','LineWidth',1.5)
    title('Interest Rate (r)')
    xlabel('Endowment (Y)')
    legend('r (\gamma_1=0.8, \gamma_2=0.5)','r (\gamma_1=1.6, \gamma_2=0.5)',Location='best')
    grid;

subplot(2,3,3)
    plot(Y(tinit:tend),sigmat(tinit:tend),'b',...
         Y0(tinit:tend),sigmat0(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Stock Return Volatility (\nu)')
    legend('\nu (\gamma_1=0.8, \gamma_2=0.5)','\nu (\gamma_1=1.6, \gamma_2=0.5)',Location='best')
    grid;


% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P2','_RRA_Fig2.pdf'));   