function z = UDaR(M,m,a) % rigorous computation of an upper bound for D_a(M,m)
tic
longprecision(32);
% Ensure m is an interval if passed as a double
if ~isa(M, 'intval'); M = intval(M); end
if ~isa(m, 'intval'); m = intval(m); end
if ~isa(a, 'intval'); a = intval(a); end

rM = a.*M./(1+M);
rm = a.*m./(1+m);

%key points in the definition of \tilde z=:zet, bm < am <ap<bp

bp = (-0.5+(1+M).*m./M)./a;
ap = (0.5+(1+M).*m./M)./a;
am = 0.5.*rM./(a.*rm)+M./rM;
bm = -0.5.*rM./(a.*rm)+M./rM;
bmm = max(-1,bm);

%computation of \tilde z

    function z = zet(s,M,m,a)
        z = M.*(s <=bm)+ (0.5.*rm.*a.*(s-bm).^2+M).*(s>bm & s <=am)+ rM.*s.*(s >am & s <=ap) + (0.5.*rM.*a.*(s-bp).^2+m).*(s>ap & s <=bp)+m.*(s>bp);
    end

%computation of the fractional function to maximize

    function z = frazet(s,M,m,a)
        z = a.^2.*zet(s-1,M,m,a)./((1+zet(s,M,m,a)).^2.*(1+zet(s-1,M,m,a)));
    end

X = intval(bmm+1); %left terminal point
Y1 = intval(ap); %one of connection points
Y2 = intval(am+1); %other connection point
x1 = sup(frazet(X,M,m,a)); % upper bounds for evaluations in the terminal points
x2 = sup(frazet(Y1,M,m,a));
x2a = sup(frazet(Y2,M,m,a));

%the next polynomial determines the critical points (if any) of frazet on the
%interval [0, \min(ap, am+1)]. If these points do not exists, we evaluate at Y

% Use only double-compatible coefficients for root-finding

c1 = -0.5 * mid(a)^2 * mid(rM) * mid(rm)^2;
c2 = 0;
c3 = -2 * mid(a) * mid(M) * mid(rm) * mid(rM);
c4 = mid(a) * mid(rm) * (1 + mid(rM) * (1 + mid(bm)));
c5 = -2 * mid(a) * mid(M)^2;

Co = [c1 c2 c3 c4 c5];
p = polynom(Co);
r = roots(Co);

x3 = zeros(1,4);
for j = 1:4
    if (isreal(r(j)) > 0) && (r(j) >= max(0, -1 - mid(bm))) && (r(j) < mid(ap) - 1 - mid(bm)) && (r(j) < mid(am) + 1 - 1 - mid(bm))
        xj = verifypoly(p, r(j)); 
        %if there exists some approx. to real root r(j) in the interval, xj gives a verified interval around it
        Xint = xj + 1 + bm;
        x3(j) = sup(frazet(Xint, M, m, a)); %upper bound for evaluations in these critical points
    else
        x3(j) = 0; %if the critical points do not exist
    end
end

% hence, max(x1,x2,x3) compute an upper bound for Da on [0, \min(ap,am+1)].
%Still, we need consider the case when ap ≥ am+1 and when ap < am+1

x4 = 0;
if sup(ap-am) >= 1
    T1 = hull(am+1, ap);
    tak1 = sup(frazet(T1, M, m, a));
    T2 = hull(ap, min(bp,1));
    tak2 = sup(frazet(T2, M, m, a));
    x4 = max(tak1, tak2);
else
    L = 32; %first, we look for the upper bound on [ap,am+1] by using a direct interval evaluation
    tak3 = zeros(1,L);
    for k = 1:L
        T3 = hull(ap + (k-1)*(am-ap+1)/L, ap + k*(am-ap+1)/L);
        tak3(k) = sup(frazet(T3, M, m, a));
    end
    H = 10; %then we look for the upper bound on [am+1, min(bp,1)]
    tak4 = zeros(1,H);
    for i = 1:H
        T4 = hull(am+1 + (i-1)*(min(bp,1)-am-1)/H, am+1 + i*(min(bp,1)-am-1)/H);
        tak4(i) = sup(frazet(T4, M, m, a));
    end
    x4 = max([tak3 tak4]);
end

v = [x1 x2 x2a x3 x4];
z = max(v);
toc
end
 