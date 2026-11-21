%% ******************************************************************************
% File: FiniteVolumeConvectionDiffusion1D.m
%   Script for simulation of steady convection-diffusion heat problem on one-
%   dimensional using the Finite Volume Method Central difference fluxes applied
%   for the diffusive terms, and either central of upwinded difference fluxes
%   applied for the convective terms. Implementation does not include source
%   terms and is limited to a uniform convective velocity u = 1 and the 
%   following two sets of boundary conditions at x = 0 and x = L:
%     1) homogeneous Dirichlet BC at west wall (Ta = 0), 
%        inhomogeneous Dirichlet BC at east wall (Tb = 1)
%     2) homogeneous Dirichlet BC at west wall (Tw = 0), 
%        inhomogeneous Neuman BC at east wall (Db = 1)
%   Linear system of equations solved directly using Matlab's backslash operator
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
%   T           :   Temperature vector at cell nodes, T(1:n)
%   A           :   Convection-diffusion system matrix, A(1:n,1:n)
%   s           :   Source vector with BC contributions, s(1:n)
%   TT          :   Temperature field extrapolated to walls, TT(1:n+2)
%   CF,DF       :   Conv. and diff. fluxes through walls, CF=[CFw,CFe]
%   GHC         :   Global heat conservation, scalar (computed from wall fluxes)
%   solve_time  :   CPU-time for solving system of linear equations A*T(:)=s(:)
%   Plots of the temperature field etc.
%
%   Author      : <<author>>
%   Date        : <<date>>
%   Version     : <<version#>>
%*******************************************************************************

%% Clear
clear all
close all
clc

%% =================================== INPUT ===================================
%% Variables and parameters
n         = 5;                  % number of cells along x,y-axis
L         = 1;                  % size of square in x,y-direction
Pe        = 7;                  % global Peclet number
problem   = 1;                  % problem to solve (1 or 2 - see file header)
fvscheme  = 'uds';              % finite volume scheme ('cds' or 'uds')

%% Exact solutions, 
% Function Texact(x) wors on arrays
if problem == 1                 % BC: Tw=0, Te=1
    Texact = @(Pe,x) (exp(Pe*x)-1)/(exp(Pe)-1);   % notes (5.8)
elseif problem == 2             % BC: Tw=0, dTe/dx=1
    Texact = @(Pe,x) 1/Pe*exp(-Pe)*(exp(Pe*x)-1); % notes (5.9)
else
    disp('problem not implemented')
end
% ================================= INPUT DONE =================================



%% ============================== PRE-PROCESSING ===============================
% << assemble matrix A and source array s based on the arguments >>
% << (n,L,Pe,problem,fvscheme), following the 6-step approach below. >>

%% 1) Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xf cell face coordinate vector (of length n+1)
%          xc cell center coordinate vector (of length n)
dx = L/n;                       % cell size in x,y-direction
xf = 0:dx:L;                    % cell face coordinate vector along x,y-axis
xc = dx/2:dx:L-dx/2;            % cell center coordinates along x-axis

%% 2) Generate velocity field
% Input  : n, xf
% Output : UFw, UFe cell faces velocities (vectors of length n)
UFw = ones(1,n);                % western face velocity, u = 1
UFe = ones(1,n);                % eastern face velocity, u = 1
%% 3) Generate convective face flux matrices
% << Input  : n, dx, Pe, UFw, UFe >>
% << Output : Fw, Fe convetive fluxes (vectors of length n), e.g. >>
% <<          is Fe the convective fluxes at the east-faces of the cells. >>
Fe = Pe * UFe;
Fw = Pe * UFw;

%% 4) Generate coefficient matrices
% << Input  : n, dx, Pe, fvscheme, Fw, Fe >>
% << Output : aW, aE, aP coefficient vectors (of length n), e.g. is aE >>
% <<          the coefficients for the cells east of the one considered. >>

D = -1/dx;

% Dirichlet-type boundary condition
Ta = 0; 
Tb = 1;
Db = 1;

if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
    % << compute coefficients for cds scheme here >>
    % << check whether the cell Peclet number dx*Pe>=2 - then print warning >>
    if dx*Pe>=2
        disp('Warning P > 2')
    end
    aW = -Fw/2 + D * ones(n,1);
    aE = Fe/2 + D * ones(n,1);
    aP = -(aE + aW) + Fe - Fw; 
 
elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
    aW = (D + min(0,-Fe(1)))* ones(n,1);
    aE = (D + min(0,Fe(1))) * ones(n,1);
    aP = -(aE + aW) + Fe - Fw;
   

else                              % fvscheme unknown
    disp('fvscheme not implemented')
end

%% 5) Impose boundary conditions and compute source vector from BC
% << Input  : n, dx, problem, aW, aE, aP >>
% << Output : s source term vector (of length n) >>
% <<          aW, aE, aP modified coefficient vectors (of lengthn) >>
s = zeros(n,1);

% We use first-order boundary conditions for both schemes. For CDS they could
% be 2nd order, but this only affects the boundary points.
aP(1) = aP(1)-aW(1); s(1)=s(1)-2*Ta*aW(1); aW(1)=0;
if problem == 1
    aP(n) = aP(n)-aE(n); s(n)=s(n)-2*Tb*aE(n); aE(n)=0;
elseif problem == 2
    aP(n) = aP(n)+aE(n); s(n)=s(n)-aE(n)*dx*Db; aE(n)=0;
end


%% 6) Assemble system matrix
% << Input  : n, aW, aE, aP  >>
% << Output : A 3-diagonal system matrix (array of size n*n) >>
% 
% << hints  : >>
% <<  - use function 'spdiags' to assemble A. The coefficient arrays >>
% <<    aW, aE, aP must in advance have been turned into column >>
% <<    vectors by setting aW=aW(:); etc. >>
% <<  - since the spdiags function places the first element of each >>
% <<    coefficient vector in the first column of the matrix, the coefficient >>
% <<    vectors aW,aE will need to be appropriately shifted and padded, >>
% <<    that is removing existing elements and adding zeros to either the start >>
% <<    or the end of the vectors corresponding to the diagonals in A, where >>
% <<    the vectors will be located.
% <<  - e.g. 1D-laplacian of report 1, ex. 2.4.1 may be assembled as follows >>
% <<     n      = 5;                                         >>
% <<     ovec   = ones(n-1,1);                               >>
% <<     Adiags = [[ovec;0;0],[1;-2*ovec;1],[0;0;ovec]];     >>
% <<     A      = spdiags(Adiags,-1:1,n+1,n+1);              >>
% ============================ PRE-PROCESSING DONE =============================

for i=2:(n-1)
        A(i,i) = aP(i,1);
        A(i,i-1) = aW(i,1);
        A(i,i+1) = aE(i,1);
end

A(1,1:2) = [aP(1), aE(1)];
A(n,n-1:n) = [aW(n), aP(n)];


%% ===================== SOLVE SYSTEM OF LINEAR EQUATIONS ======================
% Solve equations using direct backslash solve
if exist('A','var') && exist('s','var')
    solve_time = cputime;       % start timing using Matlab's cputime func.
    T = A\s(:);                 % temperature field
    solve_time = cputime-solve_time;% stop timing for linear solver
end
% ==================== SYSTEM OF LINEAR EQUATIONS SOLVED ======================



%% ============================== POST-PROCESSING ==============================
% << compute extended temperature field TT, net convective and diffusive  >>
% << fluxes through the walls CF,DF and global heat conservation GHC.  >>

%% 1) Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xf cell face coordinate vector (of length n+1)
%          xc cell center coordinate vector (of length n)
dx = L/n;                       % cell size in x,y-direction
xf = 0:dx:L;                    % cell face coordinate vector along x,y-axis
xc = dx/2:dx:L-dx/2;            % cell center coordinates along x-axis

%% 2) Generate velocity field
% Input  : n, xf
% Output : UFw, UFe cell faces velocities (vectors of length n)
UFw = ones(1,n);                % western face velocity, u = 1
UFe = ones(1,n);                % eastern face velocity, u = 1

%% 3) Extrapolate temperature field to walls
% << Input  : n, dx, fvscheme, problem, T >>
% << Output : TT extended temperature field vector (of length n+2) >>
%
% << hints: >>
% <<  - see notes section 5.3 >>
TT = zeros(1,n+2);              % initialize extended temperature array
TT(2:n+1) = T;

if strcmp(fvscheme,'cds')     % CDS FV-scheme applied
    TT(1) = Ta;
    if problem == 1
        TT(n+2) = Tb;
    elseif problem == 2
        TT(n+2) = T(n) + 0.5*dx*Db;
    end
elseif strcmp(fvscheme,'uds') % UDS FV-scheme applied
    TT(1) = 2*Ta - T(1);
    TT(n+2) = T(n);
else                              % fvscheme unknown
    disp('fvscheme not implemented')
end
   
%% 4) Compute temperature gradients at walls
% << Input  : n, dx, problem, T
% << Output : DTw, DTe, temperature wall gradients (scalars) >>
%
% << hints: >>
% <<  - see notes section 5.3 >>

DTw = 2/dx * (T(1)-Ta);
if problem == 1
    DTe = 2/dx * (Tb - T(n));
elseif problem == 2
    DTe = Db;
end
%% 5) Convective and diffusive wall fluxes and global conservation of heat
% << Input  : n, dx, TT, DTw, DTe, Pe, UFw, UFe >>
% << Output : CF=[CFw,CFe] convective fluxes at walls >>
% <<          DF=[DFw,DFe] diffusive fluxes at walls >>
% <<          GHC global heat conservation (scalar) >>
%
% << hints: >>
% <<  - see notes section 5.4 >>

CFw = UFw(1) * TT(1);
CFe = UFe(n) * TT(n+2);
CF = [CFw, CFe];

DFw = DTw / Pe;
DFe = DTe / Pe;
DF = [DFw, DFe];

GHC = CFw - DFw - (CFe - DFe);


% ============================ POST-PROCESSING DONE ============================



%% =================================== PLOTS ===================================
% << Plot both the temperature field computed by the finite volume method >>
% << and the exact temperature field Tex=Texact(Pe,xt). >>
% << Further plot the truncation error. >>

%% Coordinate arrays
% Input  : n, L
% Output : dx cell size
%          xc cell center coordinate vector (of length n)
%          xt cell center coor. vector extended with walls (of length n+2)
dx = L/n;                       % cell size in x,y-direction
xc = dx/2:dx:L-dx/2;            % cell center coordinate vector along x,y-axis
xt = [0 xc 1];                  % extended cell center coor. vector, incl. walls

%% Compute exact solution
Tex = Texact(Pe,xt);            % exact solution for given problem,Pe,n

%% Temperature field
figure(1)
subplot(211)
plot(xt,TT,'b*-') 
hold('on'), box('on'), grid('on')
plot(xt,Tex,'r--')  
xlabel('x')
ylabel('T(x),')
title(sprintf('Temperature field (problem %i, n = %i, Pe = %-3.2f, FV %s)',...
    problem,n,Pe,fvscheme))
legend('FV solution','Exact solution','location','northwest')
subplot(212)
plot(xt,TT-Tex,'b*-')
xlabel('x')
ylabel('T_{FV}(x)-T_{ex}(x)')
title('Truncation error on T(x)')
box('on'), grid('on')
print('-f1','-depsc2',...
    sprintf('TemperatureField1D_problem%i_n%i_Pe%-3.2f_FV%s.eps',...
    problem,n,Pe,fvscheme))
% ================================= PLOTS DONE =================================