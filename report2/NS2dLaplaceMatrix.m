% Syntax:
%   [A] = NS2dLaplaceMatrix(n)
%
% Input  :
%   n    :  Number of cells along x,y-axis
%
% Output :
%   A    :  Laplacian operator matrix, A(1:n^2,1:n^2)

%% Assemble Laplacian operator matrix
% << implement the Laplacian matrix with homogeneous Neumann boundary >>
% << conditions on all walls >>
%
% << hints: >>
% <<  - the coefficients from the boundary cells, in general required for >>
% <<    modifcation of the source array, may be neglected since the boundary >>
% <<    conditions are homogeneous anyway, and thus no modification will occur >>
% <<    to the source array. Hence only A is needed as output from this function.>>


function [L] = NS2dLaplaceMatrix(n)
dx = 1/n; dy = 1/n;

%derivatives:
Dwest = dy/dx;
Deast = dy/dx;
Dsouth = dx/dy;
Dnorth = dx/dy;

% coefficients

a_west(1:n, 1:n) = Dwest;
a_south(1:n, 1:n) = Dsouth;
a_north(1:n, 1:n) = Dnorth;
a_east(1:n, 1:n) = Deast;
a_P(1:n, 1:n) = -(a_west + a_east + a_south + a_north);

% West Boundary
a_P(:, 1) = a_P(:, 1) + a_west(:, 1);
a_west(:, 1) = 0;

% East Boundary
a_P(:,n) = a_P(:,n) + a_east(:,n);
a_east(:,n) = 0;

% North Boundary
a_P(n, :) = a_P(n, :) + a_north(n, :);
a_north(n, :) = 0;

% South Boundary
a_P(1, :) = a_P(1, :) + a_south(1, :);
a_south(1, :) = 0;

% vector b:
b1 = [a_west(n+1:end), zeros(n, 1)'];
b2 = [a_south(2:end), 0];
b3 = a_P(:)';
b4 = [0, a_north(1:end-1)];
b5 = [zeros(n, 1)', a_east(1:end-n)];

L = spdiags([b1; b2; b3; b4; b5]', [-n, -1:1, n], n^2, n^2);

end
