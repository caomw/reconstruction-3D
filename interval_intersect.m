function answer = intervalIntersect(ints,n)
    
    pcount = size(ints,2);
    intsT = sortrows((ints)', 1);

    opens = zeros(1,19);
    inters = zeros(2000,2);
    intersNum = 0;
    all_open = 0;
    for z=1:pcount
        opens(intsT(z,3)) = opens(intsT(z,3)) + intsT(z,2);
         if sum(opens) >= n && intsT(z,2) == 1
            all_open = 1;
         end
        
        if intsT(z,2) == 1
            if all_open == 1
                inters(intersNum+1,1) = z;
            end
        else
            if all_open == 1
                all_open = 0;
                inters(intersNum+1,2) = z;
                intersNum = intersNum+1;
            end
        end
    end
    
    answer = intsT(reshape(inters(1:intersNum,:), 1, 2*intersNum),:);
    answer = (answer)' ;

end
