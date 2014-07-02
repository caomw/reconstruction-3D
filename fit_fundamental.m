function [ ans ] = fit_fundamental( matches )
%Finds the fundamental matrix based on provided matches(with normalizing).

    syms a11 a12 a13 a21 a22 a23 a31 a32 a33;    
    FundamentalMatrix = [ a11 a12 a13 ; a21 a22 a23 ; a31 a32 a33 ];
    shouldBeMin = 0;
    for i = 1:size(matches,1)
        shouldBeMin = shouldBeMin + ( [matches(i, 1:2) 1] * FundamentalMatrix * ([matches(i, 3:4) 1])' )^2;
    end
    
    component = [a11 a12 a13 a21 a22 a23 a31 a32 a33];
    for i = 1:9
        differentiated(i) = diff(shouldBeMin, component(i));
    end
    
    for i = 1:8
        constants(i) = (-1) * subs(differentiated(i), component(1:9), [zeros(1,8) 1]);
    end
    
    coefficients = zeros(8);
    for i = 1:8
        for j = 1:8
            coefficients(i, j) = diff(differentiated(i), component(j));
        end
    end
   
 
    A = coefficients\constants';
    A = [A' 1];
    FM1 = [A(1:3); A(4:6); A(7:9)] ; 
    
    % We desire no change in z direction so...
    [U, S, V] = svd(FM1);
    S(3, 3) = 0;
    ans = U * S * V';
    
end







