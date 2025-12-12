%% Function that fits an ellipse using four known coordinates
% Name: Mees Bouwman
% Date: 19-10-2025
% Study program: Applied Physics
% Year: 3
% Assignment/Project: Smart Solutions Semester (3S) - Blood & Heat

% The mathematical derivation used for fitting an ellipse
% through four given points is described in the following math forum:
% https://math.stackexchange.com/questions/4471307/approximating-an-ellipse-given-4-points

% This function uses another function that calculates the area
% of the drawn ellipse, called '[Opp] = ellipsoppvervlak(a)'.
% This is the second function (sub-function) within this script.

function [r0, rx, ry, semi_major, semi_minor] = ellipseFourPoints(P1, P2, P3, P4)
% r0 = 'center of ellipse' a 2x1 vector
% rx and ry are the coordinates for plotting the ellipse, an array
% semi_major = 'half the of longest length', a scalar (float)
% semi_minor = 'half the of shortest length', a scalar (float)
% P1 to P4 are the four points/coordinates where a ellipse is fitted to,
% these are an 1x2 vector. P(1) is the x-component and P(2) is the
% y-component.

%% More info on what the matricis mean, please refer to Appendix B in the report
Q = [P1(1)*P1(2), P1(2)^2, P1(1), P1(2), 1;
    P2(1)*P2(2), P2(2)^2, P2(1), P2(2), 1;
    P3(1)*P3(2), P3(2)^2, P3(1), P3(2), 1;
    P4(1)*P4(2), P4(2)^2, P4(1), P4(2), 1;];

b = transpose([-P1(1)^2, -P2(1)^2, -P3(1)^2, -P4(1)^2]);

% \ is an inverse operator for matrices, the matrice before the \ sign is taken an inverse of, and ' is equal to transpose()
u = (Q * Q') \ b;
V = Q' * u;
w = null(Q);

% norm() returns the magnitude of a vector, a scalar
W = w / norm(w);

% Solving for the roots
W1 = W(1); W2 = W(2); V1 = V(1); V2 = V(2);
poly_coeffi = [W1^2, (2 * V1 * W1 - 4 * W2), (V1^2 - 4 * V2)];
roots_poly_coeffi = roots(poly_coeffi);

% t1 < t2
t1 = min(roots_poly_coeffi);
t2 = max(roots_poly_coeffi);

% t ∈ [t1 t2] represents the range of values where a valid ellipse exists
% t_optimal is the value that corresponds to the ellipse with the minimum area
% This is determined using the function ellipsoppvervlak. The whole function is given below.
Opp_function = @(t) ellipsoppvervlak(V + t * W);
t_optimaal = fminbnd(Opp_function, t1, t2);
a = transpose(V + t_optimaal * W);

A = 1; B = a(1); C = a(2); D = a(3); E = a(4); F = a(5);

q = [A, B/2;
    B/2, C];
G = [D;
    E];

r0 = -0.5 * (q \ G);   % center ellips

p = q / (transpose(r0) * q * r0 - F);

% The rotation matrix R and the diagonal D can be calculated
% analyticaly using equation 23 to 36, or just by calulating
% the eigenvalues of p using the command eig(). This retruns
% eigenvectors and can be split into [R, D].
% The second method is more efficient and compact. 
[R, D] = eig(p);

D11 = D(1,1);
D22 = D(2,2);

theta = linspace(0, 2 * pi, 1000);
z = [cos(theta) / sqrt(D11);
    sin(theta) / sqrt(D22)];

r = r0 + R * z;
rx = r(1,:);
ry = r(2,:);

a1 = 1 / sqrt(D11);
a2 = 1 / sqrt(D22);

semi_major = max(a1, a2);
semi_minor = min(a1, a2);
end

%% ---------------------------------- Sub-function ---------------------------------- %%
% Earlier in the code this was used to find the minimum single variable t
% to calulated the ellipse with the minimum area.
    
% Opp_function = @(t) ellipsoppvervlak(V + t * W);
% t_optimaal = fminbnd(Opp_function, t1, t2);
% a = transpose(V + t_optimaal * W);

function Opp = ellipsoppvervlak(a)
A = 1; B = a(1); C = a(2); D = a(3); E = a(4); F = a(5);
q = [A, B/2;
    B/2, C];
G = [D;
    E];
r0 = (-0.5) * (q \ G);    % center ellips
p = q / (transpose(r0) * q * r0 - F);
[~, D] = eig(p);
D11 = D(1,1);
D22 = D(2,2);
a = 1 / (sqrt(D11));
b = 1 / (sqrt(D22));
Opp = pi * a * b;
end

