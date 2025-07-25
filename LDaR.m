function z = LDaR(M, m, a)
% Rigorous computation of a lower bound for D_a(M, m) using INTLAB
tic
longprecision(32);

% === Input validation ===
if ~isa(M, 'intval') || ~isa(m, 'intval') || ~isa(a, 'intval')
    error('Inputs M, m, and a must be INTLAB intervals.');
end

% === Definitions ===
rM = a .* M ./ (1 + M);
rm = a .* m ./ (1 + m);

bp = (-0.5 + (1 + M) .* m ./ M) ./ a;
ap = (0.5 + (1 + M) .* m ./ M) ./ a;
am = 0.5 .* rM ./ (a .* rm) + M ./ rM;
bm = -0.5 .* rM ./ (a .* rm) + M ./ rM;
bmm = max(-1, bm);

% === Small epsilon for interval padding ===
delta = intval('1e-15');

% === Define helper functions ===
    function z = zeta(s, M, m, a)
        z = M .* (s <= bm) + ...
            (0.5 .* rm .* a .* (s - bm).^2 + M) .* (s > bm & s <= am) + ...
            rM .* s .* (s > am & s <= ap) + ...
            (0.5 .* rM .* a .* (s - bp).^2 + m) .* (s > ap & s <= bp) + ...
            m .* (s > bp);
    end

    function z = frazeta(s, M, m, a)
        z = a.^2 .* zeta(s - 1, M, m, a) ./ ((1 + zeta(s, M, m, a)).^2 .* (1 + zeta(s - 1, M, m, a)));
    end

% === Initial evaluations at endpoints ===
X = intval(bmm + 1);
Y = intval(ap);
x1 = inf(frazeta(X, M, m, a));
x2 = inf(frazeta(Y, M, m, a));

% === Polynomial setup with scalar midpoint extraction ===
a0   = mid(a);
m0   = mid(m);
M0   = mid(M);
rM0  = a0 * M0 / (1 + M0);
rm0  = a0 * m0 / (1 + m0);
bm0  = -0.5 * rM0 / (a0 * rm0) + M0 / rM0;

c1 = -0.5 * a0^2 * rM0 * rm0^2;
c2 = 0;
c3 = -2 * a0 * M0 * rm0 * rM0;
c4 = a0 * rm0 * (1 + rM0 * (1 + bm0));
c5 = -2 * a0 * M0^2;

Coefs = [c1 c2 c3 c4 c5];  % now pure double
p = polynom(Coefs);
r = roots(Coefs);  % 

% === Evaluate frazeta at valid roots ===
x3 = zeros(1, 4);  % Preallocate
for j = 1:4
    rj = r(j);
    if isreal(rj) && ...
       rj >= max(0, -1 - bm) && ...
       rj < ap - 1 - bm && ...
       rj < am + 1 - 1 - bm

        Sj = verifypoly(p, rj) + 1 + bm;
        x3(j) = inf(frazeta(Sj, M, m, a));
    else
        x3(j) = 0;
    end
end

% === Case distinction for remaining domain ===
if (ap - am >= 1)
    T1 = hull(am + 1 + delta, ap - delta);
    tak1 = inf(frazeta(T1, M, m, a));
    
    T2 = hull(ap + delta, min(bp, intval(1)));
    tak2 = inf(frazeta(T2, M, m, a));
    
    x4 = max(tak1, tak2);

else
    s_range = ap + (am - ap + 1);

    T3a = hull(ap + delta, ap + s_range / 8);
    T3b = hull(ap + s_range / 8, ap + s_range / 4);
    T3c = hull(ap + s_range / 4, ap + s_range / 2);
    T3d = hull(ap + s_range / 2, am + 1 - delta);

    tak3a = inf(frazeta(T3a, M, m, a));
    tak3b = inf(frazeta(T3b, M, m, a));
    tak3c = inf(frazeta(T3c, M, m, a));
    tak3d = inf(frazeta(T3d, M, m, a));

    bp1 = min(bp, intval(1));
    T4_width = bp1 - (am + 1);

    T4a = hull(am + 1 + delta, am + 1 + T4_width / 8);
    T4b = hull(am + 1 + T4_width / 8, am + 1 + T4_width / 4);
    T4c = hull(am + 1 + T4_width / 4, am + 1 + T4_width / 2);
    T4d = hull(am + 1 + T4_width / 2, bp1);

    tak4a = inf(frazeta(T4a, M, m, a));
    tak4b = inf(frazeta(T4b, M, m, a));
    tak4c = inf(frazeta(T4c, M, m, a));
    tak4d = inf(frazeta(T4d, M, m, a));

    x4 = max([tak3a, tak3b, tak3c, tak3d, tak4a, tak4b, tak4c, tak4d]);
end

% === Final aggregation ===
v = [x1, x2, x3, x4];
z = max(v);

toc
end
