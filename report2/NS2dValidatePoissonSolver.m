% Input   :
%   Nvec  :  Vector with number of cells in x,y, 
%               e.g. Nvec = [11 21 41 81 161 321];
%   K     :  Wave number for H1,H2 test functions and exact P solution
%
% Output  :
%   Plot of convergence of truncation errors for pressure 

%% Clear
clear all
close all
clc

folderName = 'plots';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

%% Input
Nvec = [11 21 41 81 161 321];     
Re = 3;
K    = 2*pi;   
%% Function handles for H1,H2,Pexact
H1  = @(x,y,k) -k*sin(k*x).*cos(k*y); % H1 test function
H2  = @(x,y,k) -k*cos(k*x).*sin(k*y); % H2 test function
Pe  = @(x,y,k) cos(k*x).*cos(k*y);    % pressure field exact solution

%% Initialize arrays
ERROR = zeros(1,length(Nvec));        % alloc. array for inf-norm of trunc. errors

%% Compute truncation error on pressure for difference number of cells
% << Loop over the different number of cells in Nvec. For each n compute: >>
for i=1:length(Nvec)
    % <<  - the necessary mesh arrays with StaggeredMesh2dSquare >>
    n = Nvec(i);
    [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi]=StaggeredMesh2dSquare(n);
    dy = dx;

    % <<  - the test functions h1,h2 and the exact pressure pe >>

    H1_exact = H1(Xu, Yu, K);
    H2_exact = H2(Xv, Yv, K);
    P_exact = Pe(Xp, Yp, K);

    % <<  - the FV Laplacian operator matrix A with your own NS2dLaplaceMatrix >>
    A = NS2dLaplaceMatrix(n);

    % <<  - the FV source array s incl. homogeneous Neumann BC >>
    H1e = H1_exact(2:n+1,2:n+1);
    H1w = H1_exact(2:n+1,1:n);
    
    H2n = H2_exact(2:n+1,2:n+1);
    H2s = H2_exact(1:n,2:n+1);

    % <<  - the discrete pressure field by solving A*p = s(:); >>
    s = (H1e-H1w)*dy + (H2n - H2s)*dx;
    P = zeros(n);
    P(:) = A\s(:);

    % set reference level for the solution
    P_ref = P(round(n/2),round(n/2)); % pressure @ the center of the domain
    P_relative = P-P_ref+1; % normalize so P at the center = 1

    % <<  - the maximum truncation errors on the pressure p compared to pe >>
    error = P_relative-P_exact;
    ERROR(1,i) = max(max(abs(error)));

    % << Plot P,Pe,H1,H2,etc for one particular n to visualize the results >>
    if i == 3
        figure(101)
        % --- P ---
        mesh(Xp, Yp, P_relative)
        hold on
        contour3(Xp, Yp, P_relative, 10, 'k') % Contour lines on base
        hold off
        view(3)
        grid on
        xlabel('x'), ylabel('y'), zlabel('$P$', 'Interpreter', 'latex')
        saveas(gcf, 'plots/P_validation.eps'); 
        
        % --- P_error ---
        figure(102)
        mesh(Xp, Yp, error)
        hold on
        contour3(Xp, Yp, error, 10, 'k') % Contour lines on base
        hold off
        view(3)
        grid on
        xlabel('x'), ylabel('y'), zlabel('$P_{error}$', 'Interpreter', 'latex')
        saveas(gcf, 'plots/error_validation.eps'); 
    end
end
%
% << Hints: >>
% <<  - use max(max(abs(error))) to compute max. trunc. error, or >>
% <<    be sure to use a induced vector such as norm(error(:),inf) !! >>
% <<  - set the pressure in one cell to establish a reference level  >>
% <<    since the hom. Neumann BCs do not provide this. E.g. fix the  >>
% <<    pressure at the center to 1 by setting p=p-p(ceil(n/2),ceil(n/2))+1;  >>

%% Compute rate of convergence for FV-scheme
% << Compute the asymptotic rate of convergence for pressure p >>
% << See guidelines and hints in the template files to the previous reports >>
p = polyfit(log10(Nvec), log10(ERROR), 1);
ROC = p(1);
TL = polyval(p, log10(n));


%% Plot
% << Plot the truncation error and rate of convergence >>

figure(1)
p = polyfit(log10(Nvec), log10(ERROR), 1);
slope = p(1);
loglog(Nvec, ERROR, 'o-')
grid on
xlabel('N', 'Interpreter', 'LaTex', 'FontSize', 14)
ylabel('Error [-]', 'Interpreter', 'LaTex', 'FontSize', 14)
saveas(gcf, 'plots/pressure_validation.eps'); 

% << Remember: axis lables, title, grid lines, legends with rate of conv. ...>>