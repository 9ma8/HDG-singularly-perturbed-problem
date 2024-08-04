function result=FE_local_basis_1D_reference(x,basis_type,basis_index,derivative_degree)

%We will use "FE" to replace "finite element" in the comments.
%This is for the local basis functions of 1D FE.
%x: the coordinate of the point where we want to evaluate the local FE basis function.
%basis_type: the type of the FE.
%basis_type=101:1D linear FE.
%basis_index: the index of basis function to specify which basis function we want to use.
%derivative_degree:the derivative degree of the FE basis function.

if basis_type==101
%     vertices = sort(vertices, 'ascend');
    if derivative_degree==0
        
        if basis_index==1
            result=(1-x)/2;
        elseif basis_index==2
            result=(x+1)/2;
        end

    elseif derivative_degree==1
        
        if basis_index==1
            result=-1/2;
        elseif basis_index==2
            result=1/2;
        end
        
    end
    
elseif basis_type==102
  
    
    if derivative_degree==0
        
        if basis_index==1
            result=(x.^2-x)./2;
        elseif basis_index==3
            result=(x.^2+x)./2;
        elseif basis_index==2
            result=-x.^2+1;
        end

    elseif derivative_degree==1
        
        if basis_index==1
            result=(2*x-1)./2;
        elseif basis_index==3
            result=(2*x+1)./2;
        elseif basis_index==2
            result=-2*x;
        end

    elseif derivative_degree==2
        
        if basis_index==1
            result=1;
        elseif basis_index==3
            result=1;
        elseif basis_index==2
            result=-2;
        end
    end
        

elseif basis_type==103
    
    %bottom=(vertices(2)-vertices(1))^2;
    
    if derivative_degree==0
        
        if basis_index==1
            result=-9/2*((x-vertices(1))/(vertices(2)-vertices(1))).^3+9*((x-vertices(1))/(vertices(2)-vertices(1))).^2-11/2*((x-vertices(1))/(vertices(2)-vertices(1)))+1;
        elseif basis_index==4
            result=9/2*((x-vertices(1))/(vertices(2)-vertices(1))).^3-9/2*((x-vertices(1))/(vertices(2)-vertices(1))).^2+((x-vertices(1))/(vertices(2)-vertices(1)));
        elseif basis_index==2
            result=27/2*((x-vertices(1))/(vertices(2)-vertices(1))).^3-45/2*((x-vertices(1))/(vertices(2)-vertices(1))).^2+9*((x-vertices(1))/(vertices(2)-vertices(1)));
        elseif basis_index==3
            result=-27/2*((x-vertices(1))/(vertices(2)-vertices(1))).^3+18*((x-vertices(1))/(vertices(2)-vertices(1))).^2-9/2*((x-vertices(1))/(vertices(2)-vertices(1)));
        end

    elseif derivative_degree==1
        
        if basis_index==1
            result=-27/(2*(vertices(2)-vertices(1)).^3)*(x-vertices(1)).^2+18/((vertices(2)-vertices(1)).^2)*(x-vertices(1))-11/(2*(vertices(2)-vertices(1)));
        elseif basis_index==4
            result=27/(2*(vertices(2)-vertices(1)).^3)*(x-vertices(1)).^2-9/((vertices(2)-vertices(1)).^2)*(x-vertices(1))+1/(vertices(2)-vertices(1));
        elseif basis_index==2
            result=81/(2*(vertices(2)-vertices(1)).^3)*(x-vertices(1)).^2-45/((vertices(2)-vertices(1)).^2)*(x-vertices(1))+9/(vertices(2)-vertices(1));
        elseif basis_index==3
            result=-81/(2*(vertices(2)-vertices(1)).^3)*(x-vertices(1)).^2+36/((vertices(2)-vertices(1)).^2)*(x-vertices(1))-9/(2*(vertices(2)-vertices(1)));
        end

    elseif derivative_degree==2
        
        if basis_index==1
            result=-27/((vertices(2)-vertices(1)).^3)*(x-vertices(1))+18/((vertices(2)-vertices(1)).^2);
        elseif basis_index==4
            result=27/((vertices(2)-vertices(1)).^3)*(x-vertices(1))-9/((vertices(2)-vertices(1)).^2);
        elseif basis_index==2
            result=81/((vertices(2)-vertices(1)).^3)*(x-vertices(1))-45/((vertices(2)-vertices(1)).^2);
        elseif basis_index==3
            result=-81/((vertices(2)-vertices(1)).^3)*(x-vertices(1))+36/((vertices(2)-vertices(1)).^2);
        end
    elseif derivative_degree==3
        
        if basis_index==1
            result=-27/((vertices(2)-vertices(1)).^3);
        elseif basis_index==4
            result=27/((vertices(2)-vertices(1)).^3);
        elseif basis_index==2
            result=81/((vertices(2)-vertices(1)).^3);
        elseif basis_index==3
            result=-81/((vertices(2)-vertices(1)).^3);
        end
    end
elseif basis_type==104
    
    %bottom=(vertices(2)-vertices(1))^2;
    
    if derivative_degree==0
        
        if basis_index==1
            result=32/(3*(vertices(2)-vertices(1)).^4)*(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(2));
        elseif basis_index==5
            result=32/(3*(vertices(2)-vertices(1)).^4)*(x-vertices(1))*(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-(vertices(1)+3*vertices(2))/4);
        elseif basis_index==2
            result=-128/(3*(vertices(2)-vertices(1)).^4)*(x-vertices(1))*(x-(vertices(1)+vertices(2))/2)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(2));
        elseif basis_index==3
            result=64/((vertices(2)-vertices(1)).^4)*(x-vertices(1))*(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(2));
        elseif basis_index==4
            result=-128/(3*(vertices(2)-vertices(1)).^4)*(x-vertices(1))*(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-vertices(2));
        end

    elseif derivative_degree==1
        
        if basis_index==1
            result=32/(3*(vertices(2)-vertices(1)).^4)*((x-(vertices(2)+vertices(1))/2)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(2))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(2))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-vertices(2))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-(vertices(1)+3*vertices(2))/4));
        elseif basis_index==5
            result=32/(3*(vertices(2)-vertices(1)).^4)*((x-(vertices(2)+vertices(1))/2)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(1))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+3*vertices(2))/4)*(x-vertices(1))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-vertices(1))+(x-(vertices(2)+3*vertices(1))/4)*(x-(vertices(1)+vertices(2))/2)*(x-(vertices(1)+3*vertices(2))/4));
        elseif basis_index==2
            result=-128/(3*(vertices(2)-vertices(1)).^4)*((x-(vertices(2)+vertices(1))/2)*(x-(3*vertices(2)+vertices(1))/4)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(2)+vertices(1))/4)*(x-vertices(2))+(x-vertices(1))*(x-(vertices(2)+vertices(1))/2)*(x-vertices(2))+(x-vertices(1))*(x-(vertices(2)+vertices(1))/2)*(x-(3*vertices(2)+vertices(1))/4));
        elseif basis_index==3
            result=64/((vertices(2)-vertices(1)).^4)*((x-(3*vertices(1)+vertices(2))/4)*(x-(3*vertices(2)+vertices(1))/4)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(2)+vertices(1))/4)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(1)+vertices(2))/4)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(1)+vertices(2))/4)*(x-(3*vertices(2)+vertices(1))/4));
        elseif basis_index==4
            result=-128/(3*(vertices(2)-vertices(1)).^4)*((x-(3*vertices(1)+vertices(2))/4)*(x-(vertices(2)+vertices(1))/2)*(x-vertices(2))+(x-vertices(1))*(x-(vertices(2)+vertices(1))/2)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(1)+vertices(2))/4)*(x-vertices(2))+(x-vertices(1))*(x-(3*vertices(1)+vertices(2))/4)*(x-(vertices(2)+vertices(1))/2));
        end
    end
% 
%     elseif derivative_degree==2
%         
%         if basis_index==1
%             result=-27/((vertices(2)-vertices(1)).^3)*(x-vertices(1))+18/((vertices(2)-vertices(1)).^2);
%         elseif basis_index==2
%             result=27/((vertices(2)-vertices(1)).^3)*(x-vertices(1))-9/((vertices(2)-vertices(1)).^2);
%         elseif basis_index==3
%             result=81/((vertices(2)-vertices(1)).^3)*(x-vertices(1))-45/((vertices(2)-vertices(1)).^2);
%         elseif basis_index==4
%             result=-81/((vertices(2)-vertices(1)).^3)*(x-vertices(1))+36/((vertices(2)-vertices(1)).^2);
%         end
%     end
%     elseif derivative_degree==3
%         
%         if basis_index==1
%             result=-27/((vertices(2)-vertices(1)).^3);
%         elseif basis_index==2
%             result=27/((vertices(2)-vertices(1)).^3);
%         elseif basis_index==3
%             result=81/((vertices(2)-vertices(1)).^3);
%         elseif basis_index==3
%             result=-81/((vertices(2)-vertices(1)).^3);
%         end
%     end
end


