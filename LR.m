function q = LR(M, m)
%tic;
% This function computes L_a(M, m) using INTLAB interval arithmetic
longprecision(32);
if ~isa(M, 'intval') || ~isa(m, 'intval')
        error('Inputs M and m must be INTLAB intervals (intval)');
end
% Interval-safe constant definition
a = -intval(37) / intval(24);

% Interval expressions
rM = M ./ (1 + M);
rm = m ./ (1 + m);
B = sqrt(-2 * (1 + M) .* rm);
C = B ./ (a .* rm);

apl = (-intval(24)/intval(37)) * (1 + M - 0.5 * rM ./ rm);
bpl = (-intval(24)/intval(37)) * (1 + M + 0.5 * rM ./ rm);

% Logical conditions on intervals (using infimum/supremum)
if inf(apl) >= -1
    I13 = a .* rM - 1 - 0.5 * rM .* (1 + rM) ./ rm + ...
          log(1 + M + 0.5 * (rM).^2 ./ rm) ./ rM;
    I2 = rM ./ rm + log((rM + B) ./ (-rM + B)) ./ B;
    q = I2 + I13;
elseif inf(bpl) >= -1
    q = -1 - M - 0.5 * rM ./ rm + ...
        log(1 + M + 0.5 * (rM).^2 ./ rm) ./ rM + ...
        a .* (bpl + 1) - ...
        log( abs((apl + 1 - C) ./ (apl + 1 + C)) .* ...
             abs((bpl - apl - C) ./ (bpl - apl + C)) ) ./ B;
elseif sup(bpl) <= -1
    q = a + log(1 - a .* rM) ./ rM; 
elseif inf(apl) < -1 && sup(apl) > -1
    I13 = a .* rM - 1 - 0.5 * rM .* (1 + rM) ./ rm + ...
          log(1 + M + 0.5 * (rM).^2 ./ rm) ./ rM;
    I2 = rM ./ rm + log((rM + B) ./ (-rM + B)) ./ B;
    q1 = I2 + I13;
    q2 = -1 - M - 0.5 * rM ./ rm + ...
        log(1 + M + 0.5 * (rM).^2 ./ rm) ./ rM + ...
        a .* (bpl + 1) - ...
        log( abs((apl + 1 - C) ./ (apl + 1 + C)) .* ...
             abs((bpl - apl - C) ./ (bpl - apl + C)) ) ./ B;
     q= hull(q1,q2);
elseif inf(bpl) < -1 && sup(bpl) > -1
    q3 = -1 - M - 0.5 * rM ./ rm + ...
        log(1 + M + 0.5 * (rM).^2 ./ rm) ./ rM + ...
        a .* (bpl + 1) - ...
        log( abs((apl + 1 - C) ./ (apl + 1 + C)) .* ...
             abs((bpl - apl - C) ./ (bpl - apl + C)) ) ./ B;
    q4 = a + log(1 - a .* rM) ./ rM; 
    q=hull(q3,q4); 
else
    warning('No condition met. Returning NaN.');
    q = NaN;
end
%elapsed = toc;  % End timer and store elapsed time
%disp(['Elapsed time: ', num2str(elapsed), ' seconds']);
end
