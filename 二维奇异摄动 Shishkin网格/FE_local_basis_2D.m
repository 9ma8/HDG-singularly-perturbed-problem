function result=FE_local_basis_2D(x,y,vertices,basis_type,basis_index,der_x,der_y)


J_11=vertices(1,3)-vertices(1,1);
J_21=vertices(2,3)-vertices(2,1);
% J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));
x_hat=(2*(x-vertices(1,1))-J_11)/J_11;
y_hat=(2*(y-vertices(2,1))-J_21)/J_21;

result=zeros(length(x_hat),1);
if der_x==0 && der_y==0
    result=rectangular_reference_basis(x_hat,y_hat,basis_type,basis_index,0,0);
elseif der_x==1 && der_y==0
    result=2*(rectangular_reference_basis(x_hat,y_hat,basis_type,basis_index,1,0))./(J_11); 
elseif der_x==0 && der_y==1
    result= 2*(rectangular_reference_basis(x_hat,y_hat,basis_type,basis_index,0,1))./(J_21);
end


end