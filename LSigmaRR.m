function q = LSigmaRR(M, m)
% LSigma computes a lower bound for Sigma(m, M) using rigorous INTLAB interval arithmetic
% Inputs: M, m - must be of type 'intval'
% Output: q - lower bound for Sigma(m, M)

longprecision(32);        % High precision for INTLAB
a = -intval(37)/intval(24);
tic

% Ensure m is an interval if passed as a double
if ~isa(M, 'intval'); M = intval(M); end
if ~isa(m, 'intval'); m = intval(m); end

rM = M ./ (1 + M);
rm = m ./ (1 + m);
Za = LDaR(M, m, a) ./ a^2;
B = sqrt(2 * (1 + m) .* Za);
ami = (-intval(24)/intval(37)) * (1 + m - 0.5 * rm ./ Za);
bmi = (-intval(24)/intval(37)) * (1 + m + 0.5 * rm ./ Za);
C = a .* sqrt(0.5 * Za ./ (1 + m)) .* (bmi + 1);

% Branch selection using interval bounds
if inf(ami) >= -1
    I13 = a * rm - 1 - 0.5 * rm .* (1 + rm) ./ Za + log(1 + m + 0.5 * rm.^2 ./ Za) ./ rm;
    I2 = rm ./ Za - 2 * atan(rm ./ B) ./ B;
    q = I2 + I13;
elseif (inf(ami) < -1) && (sup(ami) > -1) 
I13 = a * rm - 1 - 0.5 * rm .* (1 + rm) ./ Za + log(1 + m + 0.5 * rm.^2 ./ Za) ./ rm;
    I2 = rm ./ Za - 2 * atan(rm ./ B) ./ B;
    qT = I2 + I13;
    qR = -1 - m - 0.5 * rm ./ Za + ...
        log(1 + m + 0.5 * rm.^2 ./ Za) ./ rm + ...
        a .* (bmi + 1) - ...
        (2 ./ Za) .* atan(C);
    q=hull(qR,qT); 
elseif inf(bmi) >= -1
    q = -1 - m - 0.5 * rm ./ Za + ...
        log(1 + m + 0.5 * rm.^2 ./ Za) ./ rm + ...
        a .* (bmi + 1) - ...
        (2 ./ Za) .* atan(C);
elseif (inf(bmi) <-1) &&   (sup(bmi) >-1)
    qD = -1 - m - 0.5 * rm ./ Za + ...
        log(1 + m + 0.5 * rm.^2 ./ Za) ./ rm + ...
        a .* (bmi + 1) - ...
        (2 ./ Za) .* atan(C);
    qE=a + log(1 - a .* rm) ./ rm;
    q=hull(qD,qE); 
else
    q = a + log(1 - a .* rm) ./ rm;
end

toc
end
