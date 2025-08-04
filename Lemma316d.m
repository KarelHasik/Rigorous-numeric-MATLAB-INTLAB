X = intval('-0.009');
j = 1;
rndold = getround;  % Save rounding mode
while (sup(X) >= -0.24836) && (sup(C(X)) < inf(X))
    setround(1);                         % upward rounding
    X = intval(sup(C(X)));
    setround(rndold);                    % restore rounding mode
    j = j + 1; 
end

disp('Final result:')
disp(['sup(X) = ', num2str(sup(X),10)])
disp(['j = ', num2str(j)])
disp(['sup(C(X)) = ', num2str(sup(C(X)),10)])


function y = C(m)
  a = -intval(37)/intval(24);
  one = intval(1);
  Am = a*m/(one + m) - one + (one + m)*log(one + m)/m;
  y = LR(Am, m);  
end

%Final result:
%sup(X) = -0.2483603907
%j = 176
%sup(C(X)) = -0.2483608815