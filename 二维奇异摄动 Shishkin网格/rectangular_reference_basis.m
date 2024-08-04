function r=rectangular_reference_basis(x,y,basis_type,basis_index,derivative_degree_x,derivative_degree_y)

%We will use "FE" to replace "finite element" in the comments.
%This is the reference FE basis function on triangle ABC where A=(0,0), B=(1,0) and C=(0,1).
%x,y: the coordinates of the point where we want to evaluate the reference FE basis function.
%basis_type: the type of the FE.
%basis_type=1:2D linear FE.
%basis_type=2:2D Lagrange quadratic FE.
%basis_index: the index of FE basis function to specify which FE basis function we want to use.
%derivative_degree_x:the derivative degree of the FE basis function with respect to x.
%derivative_degree_y:the derivative degree of the FE basis function with respect to y.
%More explanation is in my "Notes for tool box of standard triangular FE" section 1-4.


if basis_type==2
    
    if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            r=(x.*y-y.*(x.^2)-x.*(y.^2)+(x.^2).*(y.^2))/4;
        elseif basis_index==2
            r=(-x.*y-y.*(x.^2)+x.*(y.^2)+(x.^2).*(y.^2))/4;
        elseif basis_index==3
            r=(x.*y+y.*(x.^2)+x.*(y.^2)+(x.^2).*(y.^2))/4;
        elseif basis_index==4
            r=(-x.*y+y.*(x.^2)-x.*(y.^2)+(x.^2).*(y.^2))/4;
        elseif basis_index==5
            r=(-y+y.^2+y.*(x.^2)-(x.^2).*(y.^2))/2;
        elseif basis_index==6
            r=(x+x.^2-x.*(y.^2)-(x.^2).*(y.^2))/2;
         elseif basis_index==7
             r=(y+y.^2-y.*(x.^2)-(x.^2).*(y.^2))/2;
         elseif basis_index==8
              r=(-x+x.^2+x.*(y.^2)-(x.^2).*(y.^2))/2;
         elseif basis_index==9
               r=1-x.^2-y.^2+(x.^2).*(y.^2);       
        end
             
    elseif derivative_degree_x==1&&derivative_degree_y==0
 
        if basis_index==1
            r=(y-2*x.*y-y.^2+2*x.*(y.^2))/4;
        elseif basis_index==2
            r=(-y-2*x.*y+y.^2+2*x.*(y.^2))/4;
        elseif basis_index==3
            r=(y+2*x.*y+y.^2+2*x.*(y.^2))/4;
        elseif basis_index==4
            r=(-y+2*x.*y-y.^2+2*x.*(y.^2))/4;
        elseif basis_index==5
            r=x.*y-x.*(y.^2);
        elseif basis_index==6
            r=(1+2*x-y.^2-2*x.*(y.^2))/2;
         elseif basis_index==7
             r=-x.*y-x.*(y.^2);
         elseif basis_index==8
             r=(-1+2*x+y.^2-2*x.*(y.^2))/2;
         elseif basis_index==9
             r=-2*x+2*x.*(y.^2);
        end           

                      
    elseif derivative_degree_x==0&&derivative_degree_y==1
            
        if basis_index==1
            r=(x-x.^2-2*x.*y+2*y.*(x.^2))/4;
        elseif basis_index==2
            r=(-x-x.^2+2*x.*y+2*y.*(x.^2))/4;
        elseif basis_index==3
            r=(x+x.^2+2*x.*y+2*y.*(x.^2))/4;
        elseif basis_index==4
            r=(-x+x.^2-2*x.*y+2*y.*(x.^2))/4;
        elseif basis_index==5
            r=(-1+2*y+x.^2-2*y.*(x.^2))/2;
        elseif basis_index==6
            r=-x.*y-y.*(x.^2);
        elseif basis_index==7
            r=(1+2*y-x.^2-2*y.*(x.^2))/2;
        elseif basis_index==8
            r=x.*y-y.*(x.^2);
        elseif basis_index==9
            r=-2*y+2*y.*(x.^2);
        end
      
    elseif derivative_degree_x==2&&derivative_degree_y==0  
        
        if basis_index==1
            r=(-y+y.^2)/2;
        elseif basis_index==2
            r=(-y+y.^2)/2;
        elseif basis_index==3
            r=(y+y.^2)/2;
        elseif basis_index==4
            r=(y+y.^2)/2;
        elseif basis_index==5
            r=y-y.^2;
        elseif basis_index==6
            r=1-y.^2;
        elseif basis_index==7
            r=-y-y.^2;
        elseif basis_index==8
            r=1-y.^2;
        elseif basis_index==9
            r=-2+2*y.^2;
        end

    elseif derivative_degree_x==0&&derivative_degree_y==2 

        if basis_index==1
            r=(-x+x.^2)/2;
        elseif basis_index==2
            r=(x+x.^2)/2;
        elseif basis_index==3
            r=(x+x.^2)/2;
        elseif basis_index==4
            r=(-x+x.^2)/2;
        elseif basis_index==5
            r=1-x.^2;
        elseif basis_index==6
            r=-x-x.^2;
        elseif basis_index==7
            r=1-x.^2;
        elseif basis_index==8
            r=x-x.^2;
        elseif basis_index==9
            r=-2+2*x.^2;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==1 

        if basis_index==1
            r=(1-2*x-2*y+4*x.*y)/4;
        elseif basis_index==2
            r=(-1-2*x+2*y+4*x.*y)/4;
        elseif basis_index==3
            r=(1+2*x+2*y+4*x.*y)/4;
        elseif basis_index==4
            r=(-1+2*x-2*y+4*x.*y)/4;
        elseif basis_index==5
            r=x-2*x.*y;
        elseif basis_index==6
            r=-y-2*x.*y;
        elseif basis_index==7
            r=-x-2*x.*y;
        elseif basis_index==8
            r=y-2*x.*y;
        elseif basis_index==9
            r=4*x.*y;
        end 
        
     elseif derivative_degree_x==2&&derivative_degree_y==1 

        if basis_index==1
            r=(-1+2*y)/2;
        elseif basis_index==2
            r=(-1+2*y)/2;
        elseif basis_index==3
            r=(1+2*y)/2;
        elseif basis_index==4
            r=(1+2*y)/2;
        elseif basis_index==5
            r=1-2*y;
        elseif basis_index==6
            r=-2*y;
        elseif basis_index==7
            r=-1-2*y;
        elseif basis_index==8
            r=-2*y;
        elseif basis_index==9
            r=4*y;
        end 
        
      elseif derivative_degree_x==1&&derivative_degree_y==2

        if basis_index==1
            r=(-1+2*x)/2;
        elseif basis_index==2
            r=(1+2*x)/2;
        elseif basis_index==3
            r=(1+2*x)/2;
        elseif basis_index==4
            r=(-1+2*x)/2;
        elseif basis_index==5
            r=-2*x;
        elseif basis_index==6
            r=-1-2*x;
        elseif basis_index==7
            r=-2*x;
        elseif basis_index==8
            r=1-2*x;
        elseif basis_index==9
            r=4*x;
        end 
        
      elseif derivative_degree_x==2&&derivative_degree_y==2 

        if basis_index==1
            r=1;
        elseif basis_index==2
            r=1;
        elseif basis_index==3
            r=1;
        elseif basis_index==4
            r=1;
        elseif basis_index==5
            r=-2;
        elseif basis_index==6
            r=-2;
        elseif basis_index==7
            r=-2;
        elseif basis_index==8
            r=-2;
        elseif basis_index==9
            r=4;
        end 
      
    end


elseif basis_type==1

    if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            r=(1-x-y+x.*y)/4;
        elseif basis_index==2
            r=(1+x-y-x.*y)/4;
        elseif basis_index==3
            r=(1+x+y+x.*y)/4;
        elseif basis_index==4
            r=(1-x+y-x.*y)/4;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==0
        
        if basis_index==1
            r=(y-1)/4;
        elseif basis_index==2
            r=(1-y)/4;
        elseif basis_index==3
            r=(1+y)/4;
        elseif basis_index==4
            r=(-1-y)/4;    
        end

    elseif derivative_degree_x==0&&derivative_degree_y==1
        
        if basis_index==1
            r=(-1+x)/4;
        elseif basis_index==2
            r=(-1-x)/4;
        elseif basis_index==3
            r=(1+x)/4;
        elseif basis_index==4
            r=(1-x)/4;
        end
     elseif derivative_degree_x==1&&derivative_degree_y==1
         if basis_index==1
            r=1/4;
        elseif basis_index==2
            r=-1/4;
        elseif basis_index==3
            r=1/4;
        elseif basis_index==4
            r=-1/4;
         end
        
    end
       
end