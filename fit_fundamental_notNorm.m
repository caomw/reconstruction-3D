function [ ans ] = fit_fundamental_notNorm( matches )
% Finds the fundamental matrix based on provided matches (without
% normalizing).

    syms a11 a12 a13 a21 a22 a23 a31 a32 a33;    
    fundamentalMatrix = [ a11 a12 a13 ; a21 a22 a23 ; a31 a32 a33 ];
    shouldBeMin = 0;
    for i = 1:size(matches,1)
        shouldBeMin = shouldBeMin + ( [matches(i, 1:2) 1] * fundamentalMatrix * ([matches(i, 3:4) 1])' )^2;
    end
    
    
    
    component = [a11 a12 a13 a21 a22 a23 a31 a32 a33];
    for i = 1:9
        differentiated(i) = diff(shouldBeMin, component(i));
    end
    
    A = zeros(9);
    for i = 1:9
        for j = 1:9
            A(i, j) = diff(differentiated(i), component(j));
        end
    end
     
    [U, S, V] = svd(A);
    lastColumnOfSVD = V(:,end);
    FM1 = ([lastColumnOfSVD(1:3), lastColumnOfSVD(4:6), lastColumnOfSVD(7:9)])';     
    
    [U, S, V] = svd(FM1);
    S(3, 3) = 0;
    ans = U * S * V';
end



