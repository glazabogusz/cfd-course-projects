%% *****************************************************************************
% File: NS2dValidateHfunction.m
%   Script for validation of implemented function for computation of H1- and H2-
%   functions, for simulation of translating lid flow in square cavity described
%   by nondimensional incompressible unsteady 2D Navier-Stokes equations.
%   Finite Volume Method on staggered cartesian mesh with uniform square cells.
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator).
%   H1- and H2-expressions:
%        H1 = 1/dx*(1/Re*dudx-u*u)|_w^e + 1/dy*(1/Re*dudy-u*v)|_s^n
%        H2 = 1/dx*(1/Re*dvdx-u*v)|_w^e + 1/dy*(1/Re*dvdy-v*v)|_s^n
%   Test functions for u,v:
%        u  = sin(2*pi*Xu).*sin(2*pi*Yu);
%        v  = sin(2*pi*Xv).*sin(2*pi*Yv);
%   yields expressions for H1,H2, defined in the velocity cell centers, which
%   can be used for a convergence study -- confirm that you obtain the
%   expected rate of convergence.
%
% Input   :
%   Nvec  :  Vector with number of cells in x,y, 
%            e.g. Nvec = [11 21 41 81 161 321 641];
%   Re    :  Global Reynolds number 
%   K     :  Vector of wave numbers for velocity test functions
%
% Output  :
%   Plot of convergence of truncation errors for H1- and H2-functions
%
%   Author      : <<author>>
%   Date        : <<date>>
%   Version     : <<version#>>
%*******************************************************************************

%% Clear
clear all
close all
clc

%% Input
Nvec = [11 21 41 81 161 321 641];   % vector with number of cells in x,y
Re   = 3;                           % global Reynolds number
K    = pi*(1:4);                    % test-fun. wave numbers, K=[k1,k2,k3,k4]

%% Function handles for u,v,H1exact,H2exact
% << Evaluate exact relations for H1,H2 based on u,v test functions, and then >>
% << implement the relations in the function handles >>
U   = @(x,y,k1,k2) sin(k1*x).*sin(k2*y);     % hor. vel. test function
V   = @(x,y,k3,k4) sin(k3*x).*sin(k4*y);     % vert. vel. test function
H1e = @(x,y,k1,k2,k3,k4) NaN*x;  % INSERT function % exact H1 expression
H2e = @(x,y,k1,k2,k3,k4) NaN*x;  % INSERT function % exact H2 expression
% << NOTE: Insert your own expressions for H1e and H2e >>

%% Initialize arrays
ERROR = zeros(2,length(Nvec));      % alloc. array for inf-norm of trunc. errors
TL    = zeros(2,length(Nvec));      % alloc. array for inf-norm of slope lines

%% Compute truncation error on H1,H2 terms for different number of cells
% << Loop over the different number of cells in Nvec. For each n compute: >>
% <<  - the necessary mesh arrays with StaggeredMesh2dSquare >>
% <<  - the test velocities u,v and the exact relations h1e,h2e >>
% <<  - the FV-approximation h1,h2 using your own NS2dHfunctions >>
% <<  - the maximum truncation errors on both h1,h2 compared to h1e,h2e >>
%
% << Hints: >>
% <<  - use max(max(abs(error))) to compute max. trunc. error, or >>
% <<    be sure to use a induced vector such as norm(error(:),inf) (!!) >>

%% Compute the rate of convergence for the FV-scheme
% << Compute the asymptotic rates of convergence for both h1,h2 >>
% << See guidelines and hints in the template files to the previous reports >>

%% Plot
% << Plot the truncation errors and the estimated rate of convergence >>
% << Plot h1,h2,h1e,h2e for one particular n to visualize the results >>
% << Remember: axis lables, title, grid lines, legends with rates of conv. ...>>

%% Initialize arrays
ERROR = zeros(2,length(Nvec));      % alloc. array for inf-norm of trunc. errors
count=1;

for ii=Nvec
    n = ii;
    [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi]=StaggeredMesh2dSquare(n);
    dy = dx;
    
    %% Function handles for u,v,H1exact,H2exact
    U   = @(x,y,k1,k2) sin(k1*x).*sin(k2*y);
    V   = @(x,y,k3,k4) sin(k3*x).*sin(k4*y);
    U_fun = U(Xu, Yu, K(1), K(2));
    V_fun = V(Xv, Yv, K(3), K(4));
    
    H1e = @(x,y,k1,k2,k3,k4) -(k1^2 + k2^2)/Re.*sin(k1*x).*sin(k2*y) ...
        - 2*k1.*(sin(k2*y).^2).*sin(k1*x).*cos(k1*x) ...
        - (sin(k1*x).*sin(k3*x)).*(k2.*cos(k2*y).*sin(k4*y) + k4.*cos(k4*y).*sin(k2*y));
    H2e = @(x,y,k1,k2,k3,k4) -(k3^2 + k4^2)/Re.*sin(k3*x).*sin(k4*y) ...
        - 2*k4.*(sin(k3*x).^2).*sin(k4*y).*cos(k4*y) ...
        - (sin(k4*y).*sin(k2*y)).*(k1.*cos(k1*x).*sin(k3*x) + k3.*cos(k3*x).*sin(k1*x));
    
    H1_exact = H1e(Xu, Yu, K(1), K(2), K(3), K(4));
    H1_exact = H1e(Xu(2:end-1, :), Yu(2:end-1, :), K(1), K(2), K(3), K(4));
    H2_exact = H2e(Xv, Yv, K(1), K(2), K(3), K(4));
    H2_exact = H2e(Xv(:, 2:end-1), Yv(:, 2:end-1), K(1), K(2), K(3), K(4));
    
    [H1,H2] = NS2dHfunctions(n,dx,Re,U_fun,V_fun);
    
    ERROR(1, count) = max(abs(H1(:) - H1_exact(:)));
    ERROR(2, count) = max(abs(H2(:) - H2_exact(:)));
    count = count+1;
end

p1 = polyfit(log10(Nvec), log10(ERROR(1, :)), 1);
p2 = polyfit(log10(Nvec), log10(ERROR(2, :)), 1);
slope1 = p1(1);
slope2 = p2(1);

loglog(Nvec, ERROR(1, :), 'ro-')
hold on
loglog(Nvec, ERROR(2, :),'bo-')
grid on
legend('H1 error', 'H2 error', 'Interpreter', 'LaTex')
xlabel('Grid points [-]', 'Interpreter', 'LaTex', 'FontSize', 14)
ylabel('Error [-]', 'Interpreter', 'LaTex', 'FontSize', 14)
%title('H validation')
saveas(gcf, 'plots/H_validation.png');