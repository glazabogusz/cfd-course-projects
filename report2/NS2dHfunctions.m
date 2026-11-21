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

%% Compute H1,H2 functions
% << compute H1,H2 terms from notes p. 135 eq. (6.20)+(6.21) >>
%
% << hints: 
% <<  - H1 is not needed at the north+south walls. 
% <<  - H2 is not needed at the east+west walls. 
% <<  - Hence there is no reason to compute those there. Thus notice the >>
% <<    suggested sizes for H1,H2 in the function documentation. >>

end