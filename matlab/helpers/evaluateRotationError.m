function rotation_error = evaluateRotationError(R_gt,R)
    
    temp = size(size(R));
    numberSolutions = 1;
    if temp(1,2) == 3
        temp2 = size(R);
        numberSolutions = temp2(1,3);
    end
    
    if numberSolutions == 1
        
        %rotation_error = norm(rodrigues(R_gt'*R));
        rotation_error = norm( rodrigues(R_gt) - rodrigues(R) );
        
    else
        
        rotation_errors = zeros(1,numberSolutions);
        
        index = 0;
        
        for i=1:numberSolutions
            
            %Check if there is any Nan
            if ~isnan(R(1,1,i))
                index = index + 1;
                %rotation_errors(1,index) = norm(rodrigues(R_gt'*R(:,:,i)));
                rotation_errors(1,index) = norm( rodrigues(R_gt) - rodrigues(R(:,:,i)) );
            end
        end
        
        %find the smallest error (we are the most "nice" to algorithms returning multiple solutions,
        %and do the disambiguation by hand)        
        [~,minIndex] = min(rotation_errors(1,1:index));
        rotation_error = rotation_errors(1,minIndex);
        
    end

end
