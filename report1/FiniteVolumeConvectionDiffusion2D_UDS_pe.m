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

%% Clear
clear all      % Clears all variables from workspace
close all      % Closes all open figure windows
clc            % Clears the command window

%% =================================== INPUT ===================================
% Define problem parameters
n         = 10;                 % Number of cells along x,y-axis
L         = 1;                   % Size of the square domain in x,y-direction
Pe_values = 10.^[1:8];           % Pe = 10^1, 10^2, ..., 10^8                % Global Peclet number (dimensionless)
problem   = 2;                   % Problem number to solve (1 or 2 - see problem description)
fvscheme = 'uds';


% Loop for different cell numbers (n) to examine convergence
n = [10,50,100,200];
N = n;   % Store values of n to loop over them

for p = 1:length(Pe_values)
    Pe = Pe_values(p);
    count = 1;

    for ii = N
        n = ii;   % Assign current value of n to variable for each iteration
    
        %% 1) Coordinate arrays
        % Define grid spacing and cell center coordinates
        dx = L/n;                       % Cell size in x-direction
        dy = L/n;                        % Cell size in y-direction
        xf = 0:dx:L;                     % Cell face coordinates along x-axis
        yf = 0:dy:L;                     % Cell face coordinates along y-axis
        xc = dx/2:dx:L-dx/2;             % Cell center coordinates along x-axis
        yc = dy/2:dy:L-dy/2;             % Cell center coordinates along y-axis
    
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
        
        DTw = 2/dx * (T(:,1)-Ta);
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
        
        GHC = dy * sum((CFe+DFe)-(CFw+DFw)) + dx * sum((CFn+DFn)-(CFs+DFs));
    
        % Create meshgrid for plotting
        XC =[0, xc, L];
        YC = [0, yc, L];
        [XX, YY] = meshgrid(XC, YC);
    
        % Plot the temperature field
        if ii == 20
            figure(3)
            surf(XX, YY, TT)
            grid on
            xlabel('X', 'Interpreter', 'LaTex', 'FontSize', 14)
            ylabel('Y', 'Interpreter', 'LaTex', 'FontSize', 14)
            zlabel('T', 'Interpreter', 'LaTex', 'FontSize', 14)
            title([fvscheme,'Temperature field, N = 20, Problem:',num2str(problem)])
            colorbar
            colormap(cool)
        end

        % Exact soultion for problem 1
        if problem == 1
            T_exact = meshgrid((exp(Pe*xc) - 1)/(exp(Pe) - 1)); % 2-D solution should be uniform in the vertical direction and can be compared to the exact 1-D example discussed in the textbook. 
            e(scheme,count) = norm(TT(2:n+1,2:n+1) - T_exact, inf) / norm(T_exact, inf);
        end

        if problem == 2
        % convergence of heat flux along the western wall
            conv_heat_west(p,count) = sum(DFw)*dy;
        end
        [Pe_values(p),n]
        count = count + 1;
    end
    
end

%% plotting
folderName = 'plots';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
    
if problem == 2
    % plotting convergence
    figure(20)
    loglog(Pe_values, -conv_heat_west(:,1), 'LineWidth', 2)
    hold on
    loglog(Pe_values, -conv_heat_west(:,2), 'LineWidth', 2)
    loglog(Pe_values, -conv_heat_west(:,3), 'LineWidth', 2)
    loglog(Pe_values, -conv_heat_west(:,4), 'LineWidth', 2)
    hold off
    % title('Diffusive heat flux for different Pe and resolution', 'FontSize', 14, 'FontWeight', 'bold')
    legend('N = 10', 'N = 50','N = 100','N = 200', 'Location', 'Best', 'FontSize', 12)
    xlabel('Pe', 'FontSize', 12)
    ylabel('-Diffusive heat flux, western wall', 'FontSize', 12)
    grid on
    print('-depsc', 'plots/figure_2_2_Pe.eps')
end


