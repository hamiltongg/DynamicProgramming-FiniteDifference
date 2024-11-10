function [c1,c2,Y,Aknew,psi,r, N1, N2,n] = Ak_YB_fun(gamma1,gamma2,gammak,lambda0,mu0,sigma0)
%{
The Implicit Method
-
==========================================
Application: Dynamic Programming approach to solve 
two-agent heterogeneous economy: Wang(1996)'s paper
I want to solve: 1 PDE:
- PDE of Ak (value function agent k)

-I use forward and backward approx
- Boundaries in S(n+1): the right way to do that is to substitute Eq1 and EqI+1 by
  Boundaries conditions

----------------------------
Author: Hamilton Galindo Gil
Date:   2023 (May, June, July), 2024 (Nov)
Paper base: Wang(1996) 
----------------------------
%}
%=========================================
%clear; clc;
%close all;
tic;
%% STEP 1: Parameters
% A. Preferences
rho    = 0.1;   % impatience rate in discount factor e^(-rho*t)
g1 = gamma1; %0.8;   % g1=gamma1 = RRA of agent 1 (more risk averse)
g2 = gamma2; %0.5;   % g2=gamma2 = RRA of agent 2 (less risk averse)

gk = gammak; %g2;    % gk=gammak

% B. Exogeneous State Variable Dynamic (Y)
mu      = mu0; %0.05; %E[dY/Y]
sigma   = sigma0; %0.3;  %Volatility term of dY/Y

% C. Weights in Utility Function of Social Planner (U)
lambda = lambda0; % weight of utility function of agent 1 (gamma_1) in U 
a = (1-lambda)/lambda;
                                              
%% STEP 2: Discretization 
% A. State space: structured grid
Ymax = 100;   
Ymin = 1;
I = 500;                % N of points in the grid: I + 1
deltaY = (Ymax-Ymin)/I; % the distance between grid points
  Y = Ymin:deltaY:Ymax; %the vector of the state variable (grid)

% B. Boundaries of Ak
r1 = rho + g1*mu - (1/2)*g1*(1+g1)*sigma^2; % "r" in RA economy with only \gamma1 agent
r2 = rho + g2*mu - (1/2)*g2*(1+g2)*sigma^2; % "r" in RA economy with only \gamma2 agent

N1 = r1 - mu + g1*sigma^2; % c1/W1 ratio = A1^{-1/g1}
N2 = r2 - mu + g2*sigma^2; % c2/W2 ratio = A2^{-1/g2}

% I use: theta1 \in (0,1) to find this boundaries 
Akmin = N1^(-gk);  %lower bound of Ak
Akmax = N2^(-gk);  %upper bound of Ak

    %A1min = N1^(-g1);  %lower bound of A1
    %A1max = N2^(-g1);  %upper bound of A1
    %A2min = N1^(-g2);  %lower bound of A2
    %A2max = N2^(-g2);  %upper bound of A2  

%% STEP 3: Preliminaries for iteration of Ak
% A. Storage of Ak for every iteration
Akmatrix = [];   

% B. Variables from Social Planner solution: 
% Consumption
c10 = 0.1; %initial point of agent1's consumption
for i=1:I+1
    c1(i) = fsolve(@(c1) c_op(c1,a,g1,g2,Y(i)),c10); %c1
end
c2 = Y - c1; %c2

% Interest Rate (r) and Price of Risk (psi)
a1 = c2*g1 + c1*g2; % aux1
a2 = -c1*(g2^2)*(1+g1) - c2*(g1^2)*(1+g2); % aux2

r   = rho + (mu*Y)*(g1*g2)./a1 + (( (sigma*Y).^2 )/2).*(g1*g2*a2./a1.^3);
psi = (sigma*Y)*(g1*g2)./a1;

%% STEP 4: INITIAL GUESS of Ak (for every point of the state var)
% A. Initial guess of "Ak"
    % Ak = [Ak_1, Ak_2, ..., Ak_I], k={1,2}
    ak0 = Y.^0.5; % Value function of agent k   
    ak = ak0;     % row-vector of (I+1) columns
    
%% STEP 5: Iteration of Ak
maxit= 1000;
crit=10^(-6);  %the criterion to stop iteration and
               %get the solution of "Ak"
deltat = 1000; %time length (from Achdou et al (2022))

for n=1:maxit
    %% STEP-5.1: Initial point of Ak
    Ak = ak;
    Akmatrix = [Akmatrix; Ak]; %We save the initial Ak of every iteration

    %% STEP-5.2: Finite Difference (Forward/Backward Diff Approx & central)
        % A. Forward and Backward Difference
          % forward difference (A1Y, A2Y, SY)
            dAkf  = [(Ak(2:end) - Ak(1:end-1))/deltaY 0];
            % Boundary nodes (Ymax): dAkf(I+1,:)= 0
                             % it will never be used
                             % because at Ymax we use backward
                             % ghost node: (i = I+2)
          
          % backward difference (A1Y, A2Y, SY)                        
            dAkb  = [0 (Ak(2:end) - Ak(1:end-1))/deltaY];
            % Boundary nodes (Ymin): dAkb(1,:)= 0
                             % it will never be used
                             % because at Ymin we use forward
                             % ghost node: (i = 0)

          % Central difference (SYY)
            ddAkYY = [(Ak(2) - Ak(1))/deltaY^2,...
                    (Ak(3:end) - 2*Ak(2:end-1) + Ak(1:end-2))/deltaY^2,...
                    (-Ak(end) + Ak(end-1))/deltaY^2 ];
                          
    %% STEP-5.3: Upwind scheme                     
     % (A) aki: coefficient
        % Agent k
        akcoef_f = ( ((sigma*Y).^2)/(2*gk) ).*( dAkf./Ak ) + mu*Y/(1-gk) + sigma*Y.*psi/gk;    
        akcoef_b = ( ((sigma*Y).^2)/(2*gk) ).*( dAkb./Ak ) + mu*Y/(1-gk) + sigma*Y.*psi/gk;    

     % (B) Indicator Functions
        % dAkY_upwind makes a choice of forward or backward differences based on
        % based on the sign of the drift (AkY):
            % Agent k
        Ifak = akcoef_f > 0; %positive drift --> forward difference
        Ibak = akcoef_b < 0; %negative drift --> backward difference:
                      %Ifak is a logic vector: zeros and ones: 
                      %1 means "true"
        I0ak = (1-Ifak-Ibak); %when a1=0 
        
        %check: NaN or Inf at "i=1" and "i=I+1"
        check1 = [akcoef_f', akcoef_b', Ifak', Ibak', I0ak'];
                 
     % (C) Boundaries conditions   
      % To be sure that in Ymin we will use Forward
        % consistent with X1=0 
         %If(1)=1; Ib(1)=0; I0(1)=0;
      % To be sure that in Y(I+1) we will use Backward
        % consistent with Z(I+1)=0
         %If(end)=0; Ib(end)=1; I0(end)=0;
      % Already taken care of automatically
              
     % (D) The first derivative with Upwind scheme            
        AkY_Upwind = dAkf.*Ifak + dAkb.*Ibak + dAkf.*I0ak;
            %AkY_Upwind(I+1) = 0 due to dAkf(I+1)=0
 
    %% STEP-5.4: Discretization of PDE system (Ak)
       % Implicit method
       % A. Coefficients (column vectors)
        alphak = (psi.^2)/(2*gk) + r;
        
        % Agent k
        Xk = - min(akcoef_b',0)/deltaY + (1/(2*(1-gk)))*( (sigma.*Y').^2 )/deltaY^2;
        Hk =  alphak' - max(akcoef_f',0)/deltaY + min(akcoef_b',0)/deltaY...
               - (1/(1-gk))*((sigma.*Y').^2)/deltaY^2; 
        Zk =  max(akcoef_f',0)/deltaY + (1/(2*(1-gk)))*((sigma.*Y').^2)/deltaY^2;            

       % B. Matrix of coefficients: "Mk"        
        % Up Diagonal (Zk)
            updiagMk = [ 0; 0; Zk(2:end-1)];
        % Central Diagonal (Hk)
            centdiagMk = [ -1; Hk(2:end-1); -1];
        % Down Diagonal (Xk)
            lowdiagMk = [ Xk(2:end-1); 0; 0];
        % Mk
            Mk = spdiags([lowdiagMk centdiagMk updiagMk], -1:1, I+1, I+1);
            
        % See the diagonal matrix XHZ
            spy(Mk)
        
            %{ 
            look at this example to undertand how "spdiags" works
            ZZL = [1 2 3]'
            ZZC = [-1 -1 -1]'
            ZZU = [4 5 6]'
            
            ZZZ = spdiags([ZZL ZZC ZZU], -1:1, 3,3)
            %}
       
       % C. Vectors: Pk, Qk, Yk_tilde
        Itilde = diag([0 ones(1,I-1) 0]);
        Pk = 1/((1-gk)*deltat)*Itilde;
        Qk = rho/(1-gk)*Itilde;
        Yk_tilde = [N1^(-gk); (gk/(1-gk))*Ak(2:I)'.^(1 - 1/gk); N2^(-gk)];

       % D. Left-hand matrix: "Bn" 
        Bn = Pk + Qk - Mk;

       % E. Right-hand matrix: "bn" (column vector)       
        bn = Yk_tilde + Pk*Ak';
        
       % F. Solve the system of equations: finding V^n+1        
        Aknew = Bn\bn; % column vector 

    %% STEP-5.5: Update of the value function        
        Akchange = Aknew - Ak';  % since we have "Ak", we calculate "Akchange"
        ak = Aknew';            % the "new initial guess" (row vector)

    %% STEP-5.6: The optimal value function
        %We use the "Absolute-value norm" 
        %We can use others: e.g., Euclidean norm        
        dist(n) = max(abs(Akchange));
        if dist(n)<crit %crit=10^(-6)
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end
        
    % To know in what "iteration" we are
    disp(n)                 
end
toc;

end