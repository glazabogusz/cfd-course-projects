%% *****************************************************************************
% File: FiniteVolumeConvectionDiffusion2D.m
%   Script for simulation of steady convection-diffusion heat problem on two-
%   dimensional square domain using the Finite Volume Method on a cartesian,
%   structured mesh with square cells. Central difference fluxes applied for
%   the diffusive terms, and either central or upwinded difference fluxes
%   applied for the convective terms.
%   Implementation does not include source terms and is limited to the
%   following two sets of convective velocity fields and boundary conditions:
%     1) uniform flow [u,v] = [1,0], homogeneous Neumann BC. at north and south
%     walls (dTn/dy=dTs/dy=0), homogeneous Dirichlet BC at west wall (Tw=0),
%     inhomogeneous Dirichlet BC at east wall (Te=1).
%     2) corner flow [u,v] = [-x,y], homogeneous Neumann BC. at north wall
%     (dTn/dy=0), homogeneous Dirichlet BC at west wall (Tw=0), inhomogeneous
%     Dirichlet BC at east wall (Te=1), inhomogeneous Dirichlet BC
%     at south wall (Ts=x)
%   Linear system of equations solved either directly using Matlab's backslash
%   operator, or iteratively using either Jacobi, Gauss-Seidel, or SOR stationary
%   iterative methods.
%
% Input         :
%   n           :  Number of cells along x,y-axis
%   L           :  Size of square in x,y-direction
%   Pe          :  Global Peclet number
%   problem     :  Problem #: 1 or 2, selects case of convective field and BCs
%   fvscheme    :  Finite volume scheme for convection-diffusion,
%                  either 'cds' or 'uds'
%
% Output        :
%   T           :   Temperature at cell nodes, T(1:n,1:n)
%   A           :   Convection-diffusion system matrix, A(1:n^2,1:n^2)
%   s           :   Source array with BC contributions, s(1:n,1:n)
%   TT          :   Temperature field extrapolated to walls, TT(1:n+2,1:n+2)
%   CF,DF       :   Conv. and diff. fluxes through walls, CF=[CFw,CFe,CFs,CFn]
%   GHC         :   Global heat conservation, scalar (computed from wall fluxes)
%   solve_time  :   CPU-time for solving system of linear equations A*T(:)=s(:)
%   Plots of the temperature field and convective velocity field
%
% Author        : <<author>>
% Date          : <<date>>
% Version       : <<version#>>
%*******************************************************************************

%% Clear
clear all      % Clears all variables from workspace
close all      % Closes all open figure windows
clc            % Clears the command window

%% =================================== INPUT ===================================
% Define problem parameters
n = [5, 20];                 % Number of cells along x,y-axis
L         = 1;                   % Size of the square domain in x,y-direction
Pe        = 10;                  % Global Peclet number (dimensionless)
problem   = 2;                   % Problem number to solve (1 or 2 - see problem description)

folderName = 'plots';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Loop for different cell numbers (n) to examine convergence
n = [10,20,50,100,200,500];
% n = linspace(10,100,10);
N = n;   % Store values of n to loop over them

for scheme = 1:2
    if scheme == 1
        fvscheme = 'cds';
        count = 1;
    else
        fvscheme = 'uds';
        count = 1;
    end

    for ii = N
        n = ii   % Assign current value of n to variable for each iteration
    
        %% 1) Coordinate arrays
        % Define grid spacing and cell center coordinates
        dx = L/n;                       % Cell size in x-direction
        dy = L/n;                        % Cell size in y-direction
        xf = 0:dx:L;                     % Cell face coordinates along x-axis
        yf = 0:dy:L;                     % Cell face coordinates along y-axis
        xc = dx/2:dx:L-dx/2;             % Cell center coordinates along x-axis
        yc = dy/2:dy:L-dy/2;             % Cell center coordinates along y-axis

        P(count) = Pe*dx;
    
        %% 2) Generate convective velocity field at cell faces
        % Depending on the selected problem, set up velocity fields
    
        if problem == 1  % Problem 1 - Uniform flow
            UFw = ones(n);  % Western face velocity, u = 1
            UFe = ones(n);  % Eastern face velocity, u = 1
            VFs = zeros(n); % Southern face velocity, v = 0
            VFn = zeros(n); % Northern face velocity, v = 0
    
        elseif problem == 2  % Problem 2 - Corner flow
            UFw = repmat(-xf(1:n),n,1);    % Western face velocity, u = -x
            UFe = repmat(-xf(2:n+1),n,1);  % Eastern face velocity, u = -x
            VFs = repmat(xf(1:n)',1,n);    % Southern face velocity, v = y
            VFn = repmat(xf(2:n+1)',1,n);  % Northern face velocity, v = y
        else
            disp('problem not implemented') % Display error if problem number is invalid
        end
    
        clearvars aW aS aE aP aN   % Clear previously defined variables for fluxes
    
        % Initialize source array for boundary conditions
        s = zeros(n);   % Source array
    
        %% 3) Initialize fluxes and finite volume coefficients
        % Coefficients for diffusive terms (central differences) - see page 121
        Dw = -dy/dx;
        De = -dy/dx;
        Ds = -dx/dy;
        Dn = -dx/dy;
    
        % Compute fluxes based on Peclet number and velocity fields
        Fw = Pe * UFw * dy;
        Fe = Pe * UFe * dy;
        Fs = Pe * VFs * dx;
        Fn = Pe * VFn * dx;
    
        if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
            % << compute coefficients for CDS scheme here >>
            % << check whether the cell Peclet number dx*Pe>=2 - then print warning >>
            aW(1:n, 1:n) = Dw - Fw/2;
            aS(1:n, 1:n) = Ds - Fs/2;
            aN(1:n, 1:n) = Dn + Fn/2;
            aE(1:n, 1:n) = De + Fe/2;
            aP(1:n, 1:n) = -(aW + aE + aS + aN) + Fe - Fw + Fn - Fs;
        elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
    
            % Apply UDS (Upwind Difference Scheme) for convective terms
            aW(1:n, 1:n) = Dw + min(0, -Fw);    % West flux term
            aS(1:n, 1:n) = Ds + min(0, -Fs);    % South flux term
            aN(1:n, 1:n) = Dn + min(0, Fn);     % North flux term
            aE(1:n, 1:n) = De + min(0, Fe);     % East flux term
            aP(1:n, 1:n) = -(aW + aE + aS + aN) + Fe - Fw + Fn - Fs; % Coefficients for the main diagonal
        else                              % fvscheme unknown
            disp('fvscheme not implemented')
        end
    
    
        %% 4) Apply boundary conditions to the system
        % West wall (Dirichlet condition)
        Ta = 0;
        aP(:, 1) = aP(:, 1) - aW(:, 1);     % Modify main diagonal for west boundary
        s(:, 1) = s(:, 1) - 2*Ta*aW(:, 1);  % Source term for west boundary
        aW(:, 1) = 0;                       % Set west flux term to 0
    
        % East wall (Dirichlet condition)
        Tb = 1;
        aP(:, n) = aP(:, n) - aE(:, n);     % Modify main diagonal for east boundary
        s(:, n) = s(:, n) - 2*Tb*aE(:, n);  % Source term for east boundary
        aE(:, n) = 0;                       % Set east flux term to 0
    
        % North wall (Neumann condition)
        Da_n = 0;
        aP(n, :) = aP(n, :) + aN(n, :);     % Modify main diagonal for north boundary
        s(n, :) = s(n, :) - aN(n, :)*dy*Da_n; % Source term for north boundary
        aN(n, :) = 0;                       % Set north flux term to 0
    
        % South wall (depends on problem type)
        if problem == 1
            Da_s = 0;           %No flux out or in the domain - (Neumann)
            aP(1, :) = aP(1, :) + aS(1, :);  % Modify main diagonal for south boundary
            s(1, :) = s(1, :) + aS(1, :)*dy*Da_s; % Source term for south boundary
            aS(1, :) = 0;                    % Set south flux term to 0
        elseif problem == 2
            Ts = xc;                         % Variable temperature at south wall
            aP(1, :) = aP(1, :) - aS(1, :);   % Modify main diagonal for south boundary
            s(1, :) = s(1, :) - 2*Ts.*aS(1, :); % Source term for south boundary
            aS(1, :) = 0;                    % Set south flux term to 0
        else
            disp('problem not implemented')  % Display error if problem number is invalid
        end
    
        % Build sparse system matrix A
        b1 = [aW(n+1:end), zeros(n, 1)'];
        b2 = [aS(2:end), 0];
        b3 = aP(:)';
        b4 = [0, aN(1:end-1)];
        b5 = [zeros(n, 1)', aE(1:end-n)];
        A = spdiags([b1; b2; b3; b4; b5]', [-n, -1:1, n], n^2, n^2);  % Sparse matrix A
    
        % Solve the system of equations for temperature
        if exist('A','var') && exist('s','var')
            T = zeros(n);                   % Initialize temperature field
            solve_time = tic;               % Start timing for solution
            T(:) = A\s(:);                  % Solve for temperature field
            solve_time = toc(solve_time);   % Stop timing for solution
        end
    
        %% ============================== Post-processing ==============================
    
    
    
        %% 1) Coordinate arrays
        % Input  : n, L
        % Output : dx cell size
        %          xf cell face coordinate vector (of length n+1)
        %          xc cell center coordinate vector (of length n)
        dx = L/(n);                       % Cell size in x-direction
        dy = L/(n);                       % Cell size in y-direction
        xf = 0:dx:L;                     % Cell face coordinates along x-axis
        yf = 0:dy:L;                     % Cell face coordinates along y-axis
        xc = dx/2:dx:L-dx/2;             % Cell center coordinates along x-axis
        yc = dy/2:dy:L-dy/2;             % Cell center coordinates along y-axis
        
        %% 2) Generate convective velocity field
        % Input  : n, xf
        % Output : UFw, UFe, VFs, VFn cell faces velocities (arrays of size n*n)
        
        
        %% 3) Extrapolate temperature field to walls
        % << Input  : n, dx, fvscheme, problem, T >>
        % << Output : TT extended temperature field array (of size (n+2)*(n+2) ) >>
        % 
        % << hints  : >>
        % <<  - the temperature in the corners is not uniquely defined from the >>
        % <<    extrapolation expressions, however, as the corner temperatures >>
        % <<    are not used, just leave the corners equal to zero (plotting the >>
        % <<    temperatur field with 'surf' will not show this) >>
        % <<  - see notes section 5.3 >>
        TT = zeros(n+2);                 % Initialize extended temperature array
        TT(2:n+1, 2:n+1) = T;            % Insert internal temperature values into TT
        
        % Horizontal is Dirichlet for each of the problems
        if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
            TT(:, 1) = Ta;                   % Set west boundary value
            TT(:, n+2) = Tb;                 % Set east boundary value
        
        elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
            TT(:, 1) = Ta;        % Set west boundary value
            TT(:, n+2) = Tb;                 % Set east boundary value
        end
        
        if problem == 1                 % problem 1 - uniform flow
            if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
                TT(1,2:n+1) = T(1,:);    % South - Neumann
                TT(n+2,2:n+1) = T(n,:);   % North - Neumann
    
                
                Ts_wall = TT(2,2:n+1);    % South - Neumann
                Tn_wall = TT(n+1,2:n+1) + 0.5*dy*Da_n;   % North - Neumann
                Tw_wall = Ta;
                Te_wall = Tb;
        
            elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
                TT(1,2:n+1) = T(1,:);    % South - Neumann
                TT(n+2,2:n+1) = T(n,:);   % North - Neumann
    
                Ts_wall = TT(2,2:n+1);    % South - Neumann
                Tn_wall = TT(n+1,2:n+1);   % North - Neumann
                Tw_wall = 2*Ta - TT(2:n+1,2);
                Te_wall = TT(2:n+1,n+1);
            end

        elseif problem == 2             % problem 2 - corner flow
            if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
                TT(1,2:n+1) = Ts;    % South - Diricht
                TT(n+2,2:n+1) = TT(n+1,2:n+1);
                TT(:,n+2) = Tb;
                
                Ts_wall = Ts;    % South - Diricht
                Tn_wall = TT(n+1,2:n+1) + 0.5*dy*Da_n;   % North - Neumann
                Tw_wall = Ta;
                Te_wall = Tb;
        
            elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
                TT(1,2:n+1) = Ts;
                TT(n+2,:) = TT(n+1,:);
                TT(:, n+2) = Tb;
    
                Ts_wall = 2*Ts - TT(2,2:n+1);    % South - Diricht
                Tn_wall = TT(n+1,2:n+1);   % North - Neumann
                Tw_wall = TT(2:n+1,2);
                Te_wall = 2*TT(2:n+1,n+2) - TT(2:n+1,n+1);
    
            end
        end
        
        %% 4) Compute temperature gradients at walls
        % << Input  : n, dx, problem, T
        % << Output : DTw, DTe, DTs, DTn temperature wall gradients (vectors of length n) >>
        %
        % << hints: >>
        % <<  - see notes section 5.3 >>
        
        DTw = 2/dx * (TT(2:n+1,2)-Ta);
        DTe = 2/dx * (Tb - T(:,n));
        
        if problem == 1                 % problem 1 - uniform flow
                DTn = 0; 
                DTs = 0;  
        
        elseif problem == 2             % problem 2 - corner flow
                DTs = 2/dy * (TT(2, 2:n+1) - Ts);     % South - Diricht
                DTn = Da_n;    % North - Neumann
        end
        
        %% 5) Convective and diffusive wall fluxes, and global conservation of heat
        % << Input  : n, dx, TT, DTw, DTe, DTs, DTn, Pe, UFw, UFe, VFs, VFn >>
        % << Output : CF=[CFw,CFe,CFs,CFn] convective fluxes summed up over walls >>
        % <<          DF=[DFw,DFe,DFs,DFn] diffusive fluxes summed up over walls >>
        % <<          GHC global heat conservation (scalar) >>
        %
        % << hints: >>
        % <<  - see notes section 5.4 >>
    
        CFw = (Pe .* UFw(:,1) .* Tw_wall)';
        CFe = (Pe .* UFe(:,n) .* Te_wall)';
        CFs = Pe .* VFs(1,:) .* Ts_wall;
        CFn = Pe .* VFn(n,:) .* Tn_wall;
        
        
        %CF=[CFw,CFe,CFs,CFn];
        
        DFw = -DTw';
        DFe = -DTe';
        DFs = -DTs;
        DFn = -DTn;
        
        DF=[DFw,DFe,DFs,DFn];
        
        GHC(scheme,count) = dy * sum((CFe+DFe)-(CFw+DFw)) + dx * sum((CFn+DFn)-(CFs+DFs));
    
        % Create meshgrid for plotting
        XC =[0, xc, L];
        YC = [0, yc, L];
        [XX, YY] = meshgrid(XC, YC);
    
        % Plot the temperature field
        if ii == 20
            figure(1)
            surf(XX, YY, TT)
            grid on
            xlabel('X', 'Interpreter', 'LaTex', 'FontSize', 14)
            ylabel('Y', 'Interpreter', 'LaTex', 'FontSize', 14)
            zlabel('T', 'Interpreter', 'LaTex', 'FontSize', 14)
            % title([fvscheme,'Temperature field, N = 20, Problem:',num2str(problem)])
            colorbar
            %colormap(cool)
            
            % if problem == 1
            % % Set the figure filename based on the value of 'scheme'
            %     if scheme == 1
            %         print('-depsc', 'plots/figure_2_1_cds.eps')
            %     elseif scheme == 2
            %         print('-depsc', 'plots/figure_2_1_uds.eps')
            %     end
            % elseif problem == 2
            %     if scheme == 1
            %         print('-depsc', 'plots/figure_2_2_cds.eps')
            %     elseif scheme == 2
            %         print('-depsc', 'plots/figure_2_2_uds.eps')
            %     end
            % end
        end


        % Exact soultion for problem 1
        if problem == 1
            T_exact = meshgrid((exp(Pe*xc) - 1)/(exp(Pe) - 1)); % 2-D solution should be uniform in the vertical direction and can be compared to the exact 1-D example discussed in the textbook. 
            e(scheme,count) = norm(TT(2:n+1,2:n+1) - T_exact, inf) / norm(T_exact, inf);
        end

        
        % convergence of heat flux along the western wall
        conv_heat_west(scheme,count) = sum(DFw)*dy;
        count = count + 1;
    end
end
%% plotting

if problem == 1
    figure(10)
    loglog(N, e(1,:), 'o-', 'LineWidth', 2) % CDS
    hold on
    loglog(N, e(2,:), 's-', 'LineWidth', 2) % UDS
    
    % --- Add slope lines ---
    p1 = 2; % Expected convergence rate for CDS (example)
    p2 = 1; % Expected convergence rate for UDS (example)
    
    % Pick anchor points for reference slopes
    x_ref = [N(1), N(end)];
    y1_ref = e(1,end) * (x_ref / N(end)).^(-p1);
    y2_ref = e(2,end) * (x_ref / N(end)).^(-p2);
    
    loglog(x_ref, y1_ref, 'o--', 'LineWidth', 1.5)
    loglog(x_ref, y2_ref, 's-.', 'LineWidth', 1.5)
    xlim([10,500]);
   
    
    hold off
    % title('Convergence of the solutions', 'FontSize', 14, 'FontWeight', 'bold')
    legend('CDS', 'UDS', ['Slope -' num2str(p1)], ['Slope -' num2str(p2)], 'Location', 'Best', 'FontSize', 12)
    xlabel('N', 'FontSize', 12)
    ylabel('T Error', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_1_convergence.eps')

end

if problem == 1
    % plotting convergence
    figure(20)
    plot(N, conv_heat_west(1,:), 'LineWidth', 2)
    hold on
    plot(N, conv_heat_west(2,:), 'LineWidth', 2)
    hold off
    % title('Convergence of the Diffusive Heat Flux along Western Wall', 'FontSize', 14, 'FontWeight', 'bold')
    legend('CDS', 'UDS', 'Location', 'Best', 'FontSize', 12)
    xlabel('N', 'FontSize', 12)
    ylabel('Diffusive heat flux, western wall', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_1_DTw.eps')

    figure(30)
    loglog(N, abs(GHC(1,:)), 'LineWidth', 2)
    hold on
    loglog(N, abs(GHC(2,:)), 'LineWidth', 2)
    hold off
    % title('Convergence of the Diffusive Heat Flux along Western Wall', 'FontSize', 14, 'FontWeight', 'bold')
    legend('CDS', 'UDS', 'Location', 'Best', 'FontSize', 12)
    xlabel('N', 'FontSize', 12)
    ylabel('|GHC|', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_1_GHC.eps')
end

if problem == 2
    % plotting convergence
    figure(20)
    plot(N, conv_heat_west(1,:), 'LineWidth', 2)
    hold on
    plot(N, conv_heat_west(2,:), 'LineWidth', 2)
    hold off
    % title('Convergence of the Diffusive Heat Flux along Western Wall', 'FontSize', 14, 'FontWeight', 'bold')
    legend('CDS', 'UDS', 'Location', 'Best', 'FontSize', 12)
    xlabel('N', 'FontSize', 12)
    ylabel('Diffusive heat flux, western wall', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_2_DTw.eps')

    figure(30)
    plot(N, GHC(1,:), 'LineWidth', 2)
    hold on
    plot(N, GHC(2,:), 'LineWidth', 2)
    hold off
    % title('Convergence of the Diffusive Heat Flux along Western Wall', 'FontSize', 14, 'FontWeight', 'bold')
    legend('CDS', 'UDS', 'Location', 'Best', 'FontSize', 12)
    xlabel('N', 'FontSize', 12)
    ylabel('GHC', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_2_GHC.eps')
end


