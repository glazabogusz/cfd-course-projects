%% << Function for computation of H1,H2 terms >>
% << Implement as script -- then turn into function >>
% << Benchmark results for development, given input: >>
% <<   n=7;                                                         >>
% <<   Re=3;                                                        >>
% <<   [dx,Xu,Yu,Xv,Yv,Xp,Yp,Xi,Yi]=StaggeredMesh2dSquare(n);       >>
% <<   u=[zeros(1,n+1);Xu(2:n+1,:).*(1-Xu(2:n+1,:));ones(1,n+1)];   >>
% <<   v=[zeros(n+1,1),Yv(:,2:n+1).*(1-Yv(:,2:n+1)),zeros(n+1,1)];  >>
% << returns: >>
% <<   H1 = [ 0       -4.932      -7.6744      -8.9436      -8.8095      -7.3421      -4.6113            0
%             0     -0.89699     -0.94947     -0.87366     -0.73955      -0.6171     -0.57629            0
%             0       -0.862     -0.89116     -0.80369     -0.66958     -0.55879      -0.5413            0
%             0     -0.82702     -0.83285     -0.73372     -0.59961     -0.50049     -0.50632            0
%             0     -0.79203     -0.77454     -0.66375     -0.52964     -0.44218     -0.47133            0
%             0     -0.75705     -0.71623     -0.59378     -0.45967     -0.38387     -0.43635            0
%             0       27.945       25.342       24.143       24.277       25.674       28.265            0 ];  >>
% << H2 = [      0            0            0            0            0            0            0
%           -4.932     -0.89699       -0.862     -0.82702     -0.79203     -0.75705      -4.7221
%          -7.6744     -0.94947     -0.89116     -0.83285     -0.77454     -0.71623      -7.3246
%          -8.9436     -0.87366     -0.80369     -0.73372     -0.66375     -0.59378      -8.5238
%          -8.8095     -0.73955     -0.66958     -0.59961     -0.52964     -0.45967      -8.3897
%          -7.3421      -0.6171     -0.55879     -0.50049     -0.44218     -0.38387      -6.9922
%          -4.6113     -0.57629      -0.5413     -0.50632     -0.47133     -0.43635      -4.4014
%                0            0            0            0            0            0            0 ]; >>


function [H1,H2] = NS2dHfunctions(n,dx,Re,u,v)
dy = dx;

%  east
ue = (u(2:n+1, 3:n+1) + u(2:n+1, 2:n))/2;
duedx = (u(2:n+1, 3:n+1) - u(2:n+1, 2:n))/dx; % uxe in the slides

%  north 
un = (u(3:n+1, 2:n) + u(2:n, 2:n))/2;
un(n, :) = u(n+2, 2:n); 

dundy = (u(3:n+1, 2:n) - u(2:n, 2:n))/dy;
dundy(n, :) = 2/dy * (u(n+2, 2:n) - u(n+1, 2:n));

%  south 
us(2:n, :) = (u(3:n+1, 2:n) + u(2:n, 2:n))/2; 
us(1, :) = u(1, 2:n);

dusdy(2:n, :) = (u(3:n+1, 2:n) - u(2:n, 2:n))/dy;
dusdy(1, :) = 2/dy * (u(2, 2:n) - u(1, 2:n));

%  west
uw = (u(2:n+1, 2:n) + u(2:n+1, 1:n-1))/2;
duwdx = (u(2:n+1, 2:n) - u(2:n+1, 1:n-1))/dx;

%  V, north and south
vn = (v(2:n, 3:n+1) + v(2:n, 2:n))/2; 
vn(n, :) = v(n+1, 2:n); 

vs(2:n, :) = (v(2:n, 3:n+1) + v(2:n, 2:n))/2;
vs(1, :) = v(1, 2:n);

% square of velocity
u2 = zeros(n+2, n+1);
u2(2:n+1, 2:n) = (ue.*ue - uw.*uw);

% derivatives of dudx and duxy
dudx_dx = zeros(n+2, n+1);
dudx_dx(2:n+1, 2:n) = (duedx-duwdx)/dx;
dudy_dy = zeros(n+2, n+1);
dudy_dy(2:n+1, 2:n) = (dundy-dusdy)/dy; 

% u*v on u grid (north, south)
uv = zeros(n+2, n+1);
uv(2:n+1, 2:n) = un.*vn - us.*vs;

% H1
H1 = 1/Re * (dudx_dx + dudy_dy) - u2/dx - uv/dy;
H1 = H1(2:n+1, :);

%%%%%



% north, v grid

vn_v = (v(3:n+1, 2:n+1) + v(2:n, 2:n+1))/2;
dvndy_v = (v(3:n+1, 2:n+1) - v(2:n, 2:n+1))/dy;

% south, v grid

vs_v = (v(2:n, 2:n+1) + v(1:n-1, 2:n+1))/2;
dvsdy_v = (v(2:n, 2:n+1) - v(1:n-1, 2:n+1))/dy;


%V east and derivative on v grid
ve_v = (v(2:n, 3:n+1) + v(2:n, 2:n))/2;
ve_v(:, n) = v(2:n, n+2);

dvedx_v = (v(2:n, 3:n+1)-v(2:n, 2:n))/dx;
dvedx_v(:, n) = 2/dx * (v(2:n, n+2) - v(2:n, n+1));

%V west and derivative on v grid
vw_v(:, 2:n) = (v(2:n, 3:n+1) + v(2:n, 2:n))/2;
vw_v(:, 1) = v(2:n, 1);

dvwdx_v(:, 2:n) = (v(2:n, 3:n+1)-v(2:n, 2:n))/dx;
dvwdx_v(:, 1) = 2/dx * (v(2:n, 2)-v(2:n, 1));

%  U east and west on v grid
ue_v = (u(3:n+1, 2:n)+u(2:n, 2:n))/2;
ue_v(:, n) = u(2:n, n+1);

uw_v(:, 2:n) = (u(3:n+1, 2:n)+u(2:n, 2:n))/2;
uw_v(:, 1) = u(2:n, 1);

% squares
v2_v = zeros(n+1, n+2);
v2_v(2:n, 2:n+1) = (vn_v.*vn_v - vs_v.*vs_v);

dvdx_dx = zeros(n+1, n+2);
dvdx_dx(2:n, 2:n+1) = (dvedx_v - dvwdx_v)/dx; % 
dvdy_dy = zeros(n+1, n+2);
dvdy_dy(2:n, 2:n+1) = (dvndy_v-dvsdy_v)/dy; % 


uv_v = zeros(n+1, n+2);
uv_v(2:n, 2:n+1) = (ue_v.*ve_v - uw_v.*vw_v);

H2 = 1/Re * (dvdx_dx + dvdy_dy) - v2_v/dy - uv_v/dx;
H2 = H2(:, 2:n+1);






%% *****************************************************************************
% NS2dHfunctions:
%   Computes H1- and H2-functions for the 2D incompressible Navier-Stokes 
%   equations in nondimensional form.
%   Finite Volume Method on staggered Cartesian mesh with uniform square cells. 
%   Central difference scheme (CDS) applied throughout the model (convective
%   and diffusive fluxes, gradient, divergence and Laplace operator). 
%   H1- and H2-expressions:
%        H1 = 1/dx*(1/Re*dudx-u*u)|_w^e + 1/dy*(1/Re*dudy-u*v)|_s^n
%        H2 = 1/dx*(1/Re*dvdx-u*v)|_w^e + 1/dy*(1/Re*dvdy-v*v)|_s^n
%
% Syntax:
%   [H1,H2] = NS2dHfunctions(n,dx,Re,u,v)
%
% Input  :
%   n    :  Number of cells along x,y-axis
%   dx   :  Cell size in x,y
%   Re   :  Global Reynolds number number
%   u    :  Horizontal velocity component array, u(1:n+2,1:n+1)
%   v    :  Vertical velocity component array, v(1:n+1,1:n+2)
%
% Output :
%   H1   :  H1-function for horizontal momentum equation, H1(1:n,1:n+1) 
%   H2   :  H2-function for vertical momentum equation,  H2(1:n+1,1:n)
%
%   Author      : <<author>>
%   Date        : <<date>>
%   Version     : <<version#>>
%*******************************************************************************

%% Compute H1,H2 functions
% << compute H1,H2 terms from notes p. 135 eq. (6.20)+(6.21) >>
%
% << hints: 
% <<  - H1 is not needed at the north+south walls. 
% <<  - H2 is not needed at the east+west walls. 
% <<  - Hence there is no reason to compute those there. Thus notice the >>
% <<    suggested sizes for H1,H2 in the function documentation. >>

end