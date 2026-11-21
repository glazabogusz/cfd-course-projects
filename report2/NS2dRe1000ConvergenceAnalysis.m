%% *****************************************************************************
% File: NS2dRe1000ConvergenceAnalysis.m
%   Driver for analysis of convergence of velocity and pressure profiles in 
%   centerlines to benchmark results from Bruneau & Saad [2006]. Results obtained 
%   for simulation of translating lid flow in square cavity described by 
%   incompressible unsteady 2D Navier-Stokes equations in nondimensional form.
%   Finite Volume Method on staggered Cartesian mesh with uniform square cells.
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator).
%   NS equations integrated in time using explicit Euler time stepping.
%   SIMPLE scheme used to derive Poisson equation for pressure. Solved directly
%   by LU-factorizing Laplacian matrix before time-stepping, then forward- and
%   backward-substituting using Matlab's backslash operator in each time step.
%
% Input        :
%   Nvec       :  Vector with number of cells along x,y-axis,
%                   e.g. Nvec = [11 21 41 81 161];
%   Re         :  Global Reynolds number number
%   Data files with (u,v,p)-fields corresponding to Nvec,Re
%
% Output       :
%   Plots:        i)  plot of velocity and pressure profiles in centerlines, 
%                 ii) plot of max errors, including convergence rates
%
%   Author      : <<author>>
%   Date        : <<date>>
%   Version     : <<version#>>
%*******************************************************************************

%% Clear
clear all
close all
clc

%% Bruneau and Saad, 2006 benchmark results

% Vert.[ y        u        p       ]% results in vertical centerline
vert = [ 1.0000 -1.00000  0.052971
     0.9688 -0.58031  0.051493
     0.9531 -0.47239  0.050314
     0.7344 -0.18861  0.012113
     0.5000  0.06205  0.000000
     0.2813  0.28040  0.040381
     0.1016  0.30029  0.104416
     0.0625  0.20227  0.109160
     0.0000  0.00000  0.110560 ];
% Hor. [ x        v        p        ]% results in horizontal centerline
hor  = [ 0.0000  0.00000  0.077429
     0.0391 -0.29330  0.078658
     0.0547 -0.41018  0.077128
     0.1406 -0.42634  0.049004
     0.5000  0.02580  0.000000
     0.7734  0.33398  0.047259
     0.9062  0.33290  0.084369
     0.9297  0.29627  0.087625
     1.0000  0.00000  0.090448 ]; 

%% Loop inputs
% Nvec = [11 21 41 81 161];       % vector with number of cells in x,y
Nvec = [11 21 41 81 161];
Re         = 1e3;                   % global Reynolds number
Ulid       = -1;                % lid velocity (+-1)
maxstep    = 1e6;               % maximum number of explicit Euler time steps 
steadytol  = 1e-6;              % tolerance on dUmax/dt to reach steady solution
statstride = 1e2;               % stride for writting status to screen

correction = 1;                 % correction not to include round-off error
step = 0;

%% Preallocate storage
u_interp_all = cell(1, length(Nvec));
v_interp_all = cell(1, length(Nvec));
px_interp_all = cell(1, length(Nvec));
py_interp_all = cell(1, length(Nvec));
ERROR = zeros(4, length(Nvec));
x_all = cell(1, length(Nvec));
y_all = cell(1, length(Nvec));
u_all = cell(1, length(Nvec));
v_all = cell(1, length(Nvec));
p_all = cell(1, length(Nvec));

ERROR = zeros(4,length(Nvec));      % allocate for errors on u,v,p
ROC   = zeros(4,1);                 % allocate for polyfit coeff.
TC    = zeros(4,length(Nvec));      % allocate for polyval tend. lines
%% << Obtain simulation results >>
% << Use your own NavierStokes2dMovingLidSquareCavityFlow script to obtain >>
% << simulation results for (u,v,p) for each n in Nvec using Re=1e3 and >>
% << Ulid=-1. Store each of the simulations in a datafile with n in the >>
% << filename (so it can be accessed in a loop). >>

%% Loop over Nvec to extract centerline results and plot profiles
for case_num = 1:length(Nvec)
    n = Nvec(case_num);
    display(n)
%% Grid arrays
% Assemble staggered grid coordinate arrays
    [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi] = StaggeredMesh2dSquare(n); 
    dy = dx;
    dt = min((Re*dx^2)/4,1/(Ulid^2*Re)); % time step based on stability analysis
    
    A = NS2dLaplaceMatrix(n);
    
    reorder = 1;
    p=colamd(A); % find a permutation index map p
    [L,U]=lu(A(p,p)); % factor B=P’AP
    
    % Set initial flow field such that du/dx+dv/dy = 0
    u = zeros(n+2,n+1);             % initial horizontal velocity component
    v = zeros(n+1,n+2);             % initial vertical velocity component
    u(n+2,:) = Ulid;                % set lid velocity in u-array
    
    % Various other arrays
    P          = zeros(n);          % allocate pressure array s.t. P(:) = A\s(:);
    gmchist    = zeros(1,maxstep);  % allocate global mass conservation hist. vector
    cmchist    = zeros(1,maxstep);  % allocate max. cell mass conserva. hist. vector
    steadyhist = ones(1,maxstep);   % allocate global steady solution hist. vector
    
    %% Assemble and factorize Laplacian operator matrix
    % << Use your own function NS2dLaplaceMatrix to assemble Laplacian matrix >>
    % << Once model is working, do re-ordering + lu-factorization (use colamd+lu) >>
    
    %% Time integrate NS-equations to reach steady solution
    % << Make a for/while loop over the maximum number of allowed time steps. >>
    tic;  % Start timer
    for i=1:maxstep
    % << At each time step, do the following sequence of computations: >>
        step =  step + 1;
    % <<  - compute mass conservation on cell level, and from this on global level>>
        cmc = dy*(u(2:n+1,2:n+1)-u(2:n+1,1:n)) + dx*(v(2:n+1,2:n+1)-v(1:n,2:n+1));
        cmchist(1,i) = max(abs(cmc(:)));
        gmchist(1,i) = sum(cmc(:));
    
    % <<  - compute H1,H2 terms with your own function NS2dHfunctions >>
        [H1,H2] = NS2dHfunctions(n,dx,Re,u,v);
    
    % <<  - compute source array for Poisson equation incl. hom Neumann BC >>
        H1w = H1(:,1:n);
        H1e = H1(:,2:n+1);
        H2s = H2(1:n,:);
        H2n = H2(2:n+1,:);
    
        if correction == 1
            s = 1/dt * cmc + (H1e-H1w)*dy + (H2n - H2s)*dx;
        else
            s = (H1e-H1w)*dy + (H2n - H2s)*dx;
        end
        
    % <<  - solve Poisson equation for pressure - fixate it to avoid fluctuations.>>
        
        if reorder == 1
            b = s(:);
            x = zeros(size(b));
            x(p)=U\(L\b(p)); % solve for x 
            P_abs = zeros(size(s));
            P_abs(:) = x;
        else
            P_abs = zeros(n,n);
            P_abs(:) = A\s(:);
        end
    
        P_max = max(P_abs);
        P = P_abs - mean(P_abs(:)); 
    
    % <<  - compute velocity time derivatives, notes eq. (6.20)+(6.21) >>
        dP_dx = (P(:,2:n) - P(:,1:n-1)) / dx;
        dP_dy = (P(2:n,:) - P(1:n-1,:)) / dy;
        du_dt = H1(:,2:n) - dP_dx;
        dv_dt = H2(2:n,:) - dP_dy;
    
    
    % <<  - step velocity field forward using explicit Euler method >>
        u(2:n+1,2:n) = u(2:n+1,2:n) + dt * du_dt;
        v(2:n,2:n+1) = v(2:n,2:n+1) + dt * dv_dt;
        u(n+2,:) = Ulid;  % Boundary condition
    
    % <<  - compute the maximum rate of change (deviation from steady solution) >>
    % <<    from the velocity time derivatives, and store it in steadyhist >>
    
        dumax = max(abs(du_dt(:)));
        dvmax = max(abs(dv_dt(:)));
        dUmax = max(dumax,dvmax);
        steadyhist(i) = dUmax;
    
    % <<  - write some status message to screen if mod(step,statstride)==0 >>
        if mod(i,statstride) == 0
            fprintf('Step: %d, dUmax: %.4e, cmc_max: %.4e\n', i, dUmax, cmchist(i));
        end
        
    % <<  - interrupt time stepping if steadyhist(step) < steadytol (use break) >>
        if dUmax < steadytol
            fprintf('Steady state reached at step %d with dUmax = %.4e\n', i, dUmax);
            break
        end
    
    % << Extensions, once your model is working. Use if-statements to add >>
    % << extensions to the code - do not delete the original the code! >>
    % <<  - use re-ordered + lu-factorized matrices for solving the Poisson >>
    % <<    equation for the pressure. >>
    end
    
    if reorder == 1
        elapsed_time_LU(case_num) = toc;  % Stop timer and get the elapsed time
    else
        elapsed_time_no_LU(case_num) = toc;  % Stop timer and get the elapsed time
    end

    u_cell = 0.5 * (u(2:n+1,2:n+1) + u(2:n+1,1:n));   % average to vertical centers
    v_cell = 0.5 * (v(2:n+1,2:n+1) + v(1:n,2:n+1));   % average to horizontal centers
    
    % measurements in the table are made in the center of the domain
    center_idx = ceil(n/2);
    x = Xp(center_idx, :);
    y = Yp(:, center_idx);
    u_line = u_cell(:,center_idx); 
    v_line = v_cell(center_idx,:); 
    p_rel = P - P(center_idx,center_idx); % pressure relative to the center of the domain

    
% <<  - interpolate the results, if necessary, to the centerlines >>
    u_interp_all{case_num} = interp1([0; y; 1], [0; u_line; Ulid], vert(:,1), 'spline');
    v_interp_all{case_num} = interp1([0, x, 1], [0, v_line, 0], hor(:,1), 'spline');
    px_interp_all{case_num} = interp1(x, p_rel(center_idx,:), hor(:,1), 'spline');
    py_interp_all{case_num} = interp1(y, p_rel(:,center_idx), vert(:,1), 'spline');
    
% <<  - for each profile compute the max. trun. errors compared to benchmarks.>>
    ERROR(1,case_num) = max(abs(u_interp_all{case_num} - vert(:,2))); % u error
    ERROR(2,case_num) = max(abs(v_interp_all{case_num} - hor(:,2))); % v error
    ERROR(3,case_num) = max(abs(py_interp_all{case_num} - vert(:,3))); % py error
    ERROR(4,case_num) = max(abs(px_interp_all{case_num} - hor(:,3))); % px error

    x_all{case_num} = x;
    y_all{case_num} = y;
    u_all{case_num} = u_cell;
    v_all{case_num} = v_cell;
    p_all{case_num} = p_rel;

end

%% Plot benchmark results
% << Plot benchmark results for velocity and pressure profiles in centerlines,>>
% << e.g. using subplot(2,2,i) with i = 1,2,3,4 to show all four profiles. >>

% <<  - plot results on top of benchmark results in existing figure >>
%       (remember to include lables, title, legend etc. when the loop is done)>>

% Labels for each plot
labels = {'u', 'v', 'p_x', 'p_y'};  % For titles and filenames
ylabs = {'$u$', '$v$', '$p_x$', '$p_y$'};  % For y-axis labels in LaTeX
benchmarks = {vert(:,2), hor(:,2), hor(:,3), vert(:,3)};
positions = {vert(:,1), hor(:,1), hor(:,1), vert(:,1)};
interpolated = {u_interp_all, v_interp_all, px_interp_all, py_interp_all};

% Create and save 4 figures
for i = 1:4
    figure;
    hold on;
    grid on;
    
    for case_num = 1:length(Nvec)
        plot(positions{i}, interpolated{i}{case_num}, '-o', ...
            'DisplayName', ['N=' num2str(Nvec(case_num))]);
    end
    
    plot(positions{i}, benchmarks{i}, 'k--', 'DisplayName', 'Benchmark');
    
    %title([labels{i} ' profile']);
    xlabel('Position');
    ylabel(ylabs{i}, 'Interpreter', 'latex');
    
    % Add legend only for the first plot
    if i == 1
        legend show; % Display the legend for the first plot only
        legend('Location', 'northeast', 'FontSize', 12); % Position in the upper right, and enlarge the font
    end
    
    saveas(gcf, sprintf('plots/%s_profile.png', labels{i}));
end

% << Compute the asymptotic rates of convergence for all four profiles. >>
logh = log10(Nvec);  
for i = 1:4
    % Fit line in log-log space
    p = polyfit(logh, log10(ERROR(i,:)), 1);
    ROC(i) = p(1); 
    TC(i,:) = 10.^polyval(p, logh);  % Reconstruct trend line in linear space
end

%% Plot
% << Plot the truncation errors and tendency lines >>
for i = 1:4
    figure;
    loglog(Nvec, ERROR(i,:), 'o-', 'LineWidth', 2); hold on;
    loglog(Nvec, TC(i,:), '--', 'LineWidth', 1.5);
    %title([labels{i}, sprintf(' (ROC = %.2f)', ROC(i))]);
    xlabel('N'); ylabel('Truncation error');
    grid on;
    legend('Error','Trend line', 'Location', 'northeast', 'FontSize', 12); % Set location and font size
    saveas(gcf, sprintf('plots/error_vs_grid%d.png', i));  % Save as .png
end

% << Visualize streamlines and pressure field for one particular n >>
% << Remember: axis lables, title, grid lines, legends with rates of conv. etc>>