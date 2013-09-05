function [position_error,rotation_error] = evaluateTransformationError(T_gt,T)
    
    temp = size(size(T));
    numberSolutions = 1;
    if temp(1,2) == 3
        temp2 = size(T);
        numberSolutions = temp2(1,3);
    end
    
    if numberSolutions == 1
        
        position_error = norm(T_gt(:,4)-T(:,4));
        rotation_error = norm(rodrigues(T_gt(:,1:3)'*T(:,1:3)));
        
    else
        
        position_errors = zeros(1,numberSolutions);
        rotation_errors = zeros(1,numberSolutions);
        
        index = 0;
        
        for i=1:numberSolutions
            
            %Check if there is any Nan
            if ~isnan(T(1,1,i))
                index = index + 1;
                position_errors(1,index) = norm(T_gt(:,4)-T(:,4,i));
                rotation_errors(1,index) = norm(rodrigues(T_gt(:,1:3)'*T(:,1:3,i)));
            end
        end
        
        %find the smallest error (we are the most "nice" to algorithms returning multiple solutions,
        %and do the disambiguation by hand)
        [~,minIndex] = min(position_errors(1,1:index));
        position_error = position_errors(1,minIndex);
        rotation_error = rotation_errors(1,minIndex);
        
    end

end