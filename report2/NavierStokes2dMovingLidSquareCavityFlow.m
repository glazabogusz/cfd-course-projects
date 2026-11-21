%% *****************************************************************************
% File: NavierStokes2dMovingLidSquareCavityFlow.m
%   Script for simulation of translating lid flow in square cavity described by 
%   incompressible unsteady 2D Navier-Stokes equations in nondimensional form.
%   Finite Volume Method on staggered cartesian mesh with uniform square cells. 
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator). 
%   NS equations integrated in time using the explicit Euler method.
%   SIMPLE scheme used to derive Poisson equation for pressure. Solved directly 
%   by LU-factorizing Laplacian matrix before time-stepping, then forward- and 
%   backward-substituting using Matlab's backslash operator in each time step.  
%
% Input        :
%   n          :  Number of cells along x,y-axis
%   Re         :  Global Reynolds number number
%   Ulid       :  Lid velocity (+-1, set =-1 to agree Bruneau & Saad 2006)
%   maxstep    :  Maximum number of time steps for explicit Euler method.
%   dt         :  Time step for explicit Euler method.
%   steadytol  :  Tolerance on maximum rate of change for termination at steady 
%                 solution, thus dUmax/dt < steadytol, dUmax = max(dumax,dvmax)   
%   statstride :  stride for writting status to screen
%
% Output       :
%   u          :  Horizontal velocity component array, u(1:n+2,1:n+1)
%   v          :  Vertical velocity component array, v(1:n+1,1:n+2)
%   p          :  Pressure field array, p(1:n,1:n)
%   cmc        :  Control volume mass conservation array, cmc(1:n,1:n)
%   gmchist    :  Global mass conservation history, gmchist(tstep)=sum(cmc(:))                   
%   cmchist    :  Max. cell mass conser. history, cmchist(step)=norm(cmc(:),inf)
%   A          :  Laplacian system matrix, A(1:n^2,1:n^2)
%   s          :  Source array for Poisson equation for pressure, s(1:n,1:n)
%   Xu,Yu      :  Coord. arrays for u-grid, Xu(1:n+2,1:n+1), Yu(1:n+2,1:n+1)
%   Xv,Yv      :  Coord. arrays for v-grid, Xv(1:n+1,1:n+2), Yv(1:n+1,1:n+2)
%   Xp,Yp      :  Coord. arrays for p-grid, Xp(1:n,1:n), Yp(1:n,1:n)
%   Xi,Yi,     :  Coord. arrays for fd-grid, Xi(1:n+1,1:n+1), Yi(1:n+1,1:n+1)
%   Plots of   :  i) velocity field streamlines, ii) pressure field contours. 
%                 iii) gmchist and steadyhist
%   Results written to file 'NS2d_Re<>_n<>.mat'
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
n          = 20;                 % number of cells along x,y-axis
Re         = 1000;                 % global Reynolds number
Ulid       = -1;                % lid velocity (+-1)
maxstep    = 1e6;               % maximum number of explicit Euler time steps 
steadytol  = 1e-8;              % tolerance on dUmax/dt to reach steady solution
statstride = 1e2;               % stride for writting status to screen

correction = 1;                 % correction not to include round-off error
step = 0;

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

elapsed_time = toc;  % Stop timer and get the elapsed time
fprintf('Computation time: %.4f seconds\n', elapsed_time);

%% Save steady solution
% << Remember: save your steady solution (u,v,p), and all other parameters and>>
% << arrays needed to redo/restart the simulation (very important for Re=1e3!)>>
% << You may want to truncate your gmchist and steadyhist to correspond to the >>
% << number of time steps taken, i.e. gmchist = gmchist(:,1:step+1); >>
% << Name datafile containing variables using sprintf, e.g. >>
% <<    save(sprintf('NS2d_Re%-3.2e_n%i',Re,n),n,Re,dt,step,u,v,p)  >>

step = i;
gmchist = gmchist(:,1:step);
cmchist = cmchist(:,1:step);
steadyhist = steadyhist(:,1:step);
filename = sprintf('NS2d_Re%-3.2e_n%i.mat', Re, n);
save(filename, 'n','Re','Ulid','dt','step','u','v','P','cmc','gmchist','cmchist','steadyhist','A','s','Xu','Yu','Xv','Yv','Xp','Yp','Xi','Yi');

%% Plots

% << produce various plots of the results, e.g streamlines, pressure field, >>
% << and gmchist, cmchist and steadyhist to show change of sol. with time >>
%

% Interpolate u and v onto pressure grid (centered)
u_interp = 0.5 * (u(2:n+1,2:n+1) + u(2:n+1,1:n));   % average to vertical centers
v_interp = 0.5 * (v(2:n+1,2:n+1) + v(1:n,2:n+1));   % average to horizontal centers

time = dt * (1:step);
gmchist = gmchist(1:step);
cmchist = cmchist(1:step);
steadyhist = steadyhist(1:step);

% Velocity streamlines
figure
streams = streamslice(Xp, Yp, u_interp, v_interp, 2); % 2 = density of streamlines
set(streams, 'Color', 'k')
axis equal tight
xlabel('x'); ylabel('y');
%title(sprintf('Streamlines, Re = %-3.2e, n = %i, step = %i, dt = %-4.3e', Re, n, step, dt))
saveas(gcf,sprintf('plots/stream_lines%.0f.png', n));

% Pressure field
figure
contourf(Xp, Yp, P, 20)
colorbar
axis equal tight
xlabel('x'); ylabel('y');
%title(sprintf('Pressure Contours, Re = %-3.2e, n = %i, step = %i, dt = %-4.3e', Re, n, step, dt))
saveas(gcf,sprintf('plots/pressure_field%.0f.png', n));


% Plot global mass conservation history
figure
plot(time, gmchist, 'b-', 'LineWidth', 1.5)
xlabel('Time')
ylabel('Global Mass Conservation')
% title('Global Mass Conservation History')
grid on
saveas(gcf,sprintf('plots/global_mass_conservation_Re%.0f.eps', Re));

% Plot max cell mass conservation history
figure
plot(time, cmchist, 'r-', 'LineWidth', 1.5)
xlabel('Time')
ylabel('Max Cell Mass Conservation')
% title('Max Cell Mass Conservation History')
grid on
saveas(gcf, sprintf('plots/cell_mass_conservation_Re%.0f.eps', Re));

% Plot max velocity change history (steady criterion)
figure
semilogy(time, steadyhist, 'k-', 'LineWidth', 1.5)
xlabel('Time')
ylabel('Max |du/dt| or |dv/dt|')
title('Velocity Convergence History')
grid on

% << hints: >>
% <<  - scale and trim plots properly, add axis-labels, title, legend, etc. >>
% <<  - title may contain variables 'Re,n,step,dt' etc, e.g. by
% <<       title(sprintf('Re = %-3.2e, n = %i, step = %i, dt = %-4.3e',Re,n,step,dt))  >>