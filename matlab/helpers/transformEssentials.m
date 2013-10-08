function Rs = transformEssentials(Es)

    temp = size(size(Es));
    numberSolutions = 1;
    if temp(1,2) == 3
        temp2 = size(Es);
        numberSolutions = temp2(1,3);
    end
    
    if numberSolutions == 1
        
        Rs = zeros(3,3,2);
        [U,~,V] = svd(Es);
        W = [0 -1 0; 1 0 0; 0 0 1];
        Rs(1:3,1:3) = U * W * V';
        Rs(1:3,4:6) = U * W' * V';
        
        if( det(Rs(1:3,1:3)) < 0 )
            Rs(1:3,1:3) = -Rs(1:3,1:3);
        end
        if( det(Rs(1:3,4:6)) < 0 )
            Rs(1:3,4:6) = -Rs(1:3,4:6);
        end
        
    else
        
        Rs_temp = zeros(3,3,2*numberSolutions);
        index = 0;
        
        for i=1:numberSolutions
            
            %Check if there is any Nan
            if ~isnan(Es(1,1,i))
                [U,~,V] = svd(Es(:,:,i));
                W = [0 -1 0; 1 0 0; 0 0 1];
                index = index + 1;
                Rs_temp( :,:, index ) = U * W * V';
                if(det(Rs_temp( :,:, index )) < 0)
                    Rs_temp( :,:, index ) = -Rs_temp( :,:, index );
                end
                index = index + 1;
                Rs_temp( :,:, index ) = U * W' * V';
                if(det(Rs_temp( :,:, index )) < 0)
                    Rs_temp( :,:, index ) = -Rs_temp( :,:, index );
                end
            end
        end
        
        Rs = Rs_temp(:,:,1:index);
        
    end
                
end
