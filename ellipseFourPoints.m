%% Functie dat een ellipse kan fitten met vier bekende coördinaten
% Naam: Mees Bouwman
% Datum: 19-10-2025
% Opleiding: Technische Natuurkunde
% Leerjaar: 3
% Opdracht/project: Smart Solutions Semester (3S) - Blood & Heat

% De wiskundige afleiding dat is gebruikt voor het fitten van een ellips
% met vier gegeven punten staat beschreven in deze wiskunde forum.
% https://math.stackexchange.com/questions/4471307/approximating-an-ellipse-given-4-points

% Deze functie maakt gebruik van een andere functie dat de oppervlakte
% berekend van de getekende ellips genaamd '[Opp] = ellipsoppvervlak(a)'.
% Dit is de tweede functie (sub-function) in dit script.

function [r0, rx, ry, semi_major, semi_minor] = ellipseFourPoints(P1, P2, P3, P4)
%
Q = [P1(1)*P1(2), P1(2)^2, P1(1), P1(2), 1;
    P2(1)*P2(2), P2(2)^2, P2(1), P2(2), 1;
    P3(1)*P3(2), P3(2)^2, P3(1), P3(2), 1;
    P4(1)*P4(2), P4(2)^2, P4(1), P4(2), 1;];
%
b = transpose([-P1(1)^2, -P2(1)^2, -P3(1)^2, -P4(1)^2]);
%
u = (Q * Q') \ b;
V = Q' * u;
w = null(Q);
W = w / norm(w);
%
W1 = W(1); W2 = W(2); V1 = V(1); V2 = V(2);
poly_coeffi = [W1^2, (2 * V1 * W1 - 4 * W2), (V1^2 - 4 * V2)];
roots_poly_coeffi = roots(poly_coeffi);
% t1 < t2
t1 = min(roots_poly_coeffi);
t2 = max(roots_poly_coeffi);
% t ∈ [t1 t2] geeft waardes waarin een geldige ellips bestaat
% t_optimaal is de waarde die een ellips geeft met de minimale oppervlakte
% Dit is wordt met behulp van de function beneden gedaan.
Opp_function = @(t) ellipsoppvervlak(V + t * W);
t_optimaal = fminbnd(Opp_function, t1, t2);
a = transpose(V + t_optimaal * W);
%
A = 1; B = a(1); C = a(2); D = a(3); E = a(4); F = a(5);
%
q = [A, B/2;
    B/2, C];
G = [D;
    E];
%
r0 = -0.5 * (q \ G);   % center ellips
%
p = q / (transpose(r0) * q * r0 - F);
%
[R, D] = eig(p);
%
D11 = D(1,1);
D22 = D(2,2);
%
theta = linspace(0, 2 * pi, 1000);
z = [cos(theta) / sqrt(D11);
    sin(theta) / sqrt(D22)];
%
r = r0 + R * z;
rx = r(1,:);
ry = r(2,:);
%
a1 = 1 / sqrt(D11);
a2 = 1 / sqrt(D22);
%
semi_major = max(a1, a2);
semi_minor = min(a1, a2);
end

%% ----------------------------------------------------------------------------------------- %%
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

