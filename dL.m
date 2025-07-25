function q = dL(M, m)
%tic;
% This function computes partial derivative dL_a(M, m)/dm using INTLAB interval arithmetic
longprecision(32);
if ~isa(M, 'intval') || ~isa(m, 'intval')
        error('Inputs M and m must be INTLAB intervals (intval)');
end
a = -intval(37)/intval(24);
rM = M ./ (1 + M);
rm = m ./ (1 + m);
drm = 1 ./ (1 + m).^2;

B = sqrt(-2 * (1 + M) .* rm);
C = B ./ (a .* rm);
DB = -(1 + M) ./ B;

apl = (-intval(24)/intval(37)) .* (1 + M - 0.5 .* rM ./ rm);
bpl = (-intval(24)/intval(37)) .* (1 + M + 0.5 .* rM ./ rm);
dbpl = (intval(12)/intval(37)) .* rM .* drm ./ rm.^2;
dapl = -dbpl;

DC = (-intval(24)/intval(37)) .* (DB .* rm - B) .* drm ./ rm.^2;

% Check if M satisfies the domain condition
threshold = -(intval(37)/intval(24)) .* m ./ (1 + m) - 1 + (1 + m) .* log(1 + m) ./ m;
if inf(M) >= sup(threshold)
    disp('out')
    q = NaN;
    return
end

if inf(apl) >= -1
    I13 = (0.5 .* rM .* (1 + rM) ./ rm.^2 - ...
           (0.5 .* rM.^2 ./ rm.^2) .* (1 + M + 0.5 .* rM.^2 ./ rm).^(-1) ./ rM) .* drm;
    I2 = (-rM ./ rm.^2 - DB .* log((rM + B) ./ (-rM + B)) ./ B.^2 + ...
          DB .* ((rM + B).^(-1) - (-rM + B).^(-1)) ./ B) .* drm;
    q = I2 + I13;
elseif sup(bpl) > -1
    q = 0.5 .* rM .* drm ./ rm.^2 - ...
        (0.5 .* rM.^2 ./ rm.^2) .* drm .* (1 + M + 0.5 .* rM.^2 ./ rm).^(-1) ./ rM + ...
        a .* dbpl + ...
        log(abs((apl + 1 - C) ./ (apl + 1 + C)) .* abs((bpl - apl - C) ./ (bpl - apl + C))) .* ...
        DB .* drm ./ B.^2 - ...
        B.^(-1) .* ((dapl - DC) ./ (apl + 1 - C) - ...
                   (dapl + DC) ./ (apl + 1 + C) + ...
                   (dbpl - dapl - DC) ./ (bpl - apl - C) - ...
                   (dbpl - dapl + DC) ./ (bpl - apl + C));
elseif (sup(bpl) <= -1) 
    q = 0;
elseif inf(apl) < -1 && sup(apl) >= -1
I13 = (0.5 .* rM .* (1 + rM) ./ rm.^2 - ...
           (0.5 .* rM.^2 ./ rm.^2) .* (1 + M + 0.5 .* rM.^2 ./ rm).^(-1) ./ rM) .* drm;
    I2 = (-rM ./ rm.^2 - DB .* log((rM + B) ./ (-rM + B)) ./ B.^2 + ...
          DB .* ((rM + B).^(-1) - (-rM + B).^(-1)) ./ B) .* drm;
    q1 = I2 + I13;
    q2 = 0.5 .* rM .* drm ./ rm.^2 - ...
        (0.5 .* rM.^2 ./ rm.^2) .* drm .* (1 + M + 0.5 .* rM.^2 ./ rm).^(-1) ./ rM + ...
        a .* dbpl + ...
        log(abs((apl + 1 - C) ./ (apl + 1 + C)) .* abs((bpl - apl - C) ./ (bpl - apl + C))) .* ...
        DB .* drm ./ B.^2 - ...
        B.^(-1) .* ((dapl - DC) ./ (apl + 1 - C) - ...
                   (dapl + DC) ./ (apl + 1 + C) + ...
                   (dbpl - dapl - DC) ./ (bpl - apl - C) - ...
                   (dbpl - dapl + DC) ./ (bpl - apl + C));
    q=hull(q1,q2); 
else
    warning('No condition met. Returning NaN.');
    q = NaN;
end
%elapsed = toc;  % End timer and store elapsed time
%disp(['Elapsed time: ', num2str(elapsed), ' seconds']);
end
