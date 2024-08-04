function result=FE_local_basis_1D_local(x,vertices,basis_type,basis_index,der_x)


J_11=vertices(2)-vertices(1);
% J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));
x_hat=(2*x - (vertices(2)+vertices(1)))./(J_11);

result=zeros(length(x_hat),1);
if der_x==0
    result=FE_local_basis_1D_reference(x_hat,basis_type,basis_index,0);
elseif der_x==1 
    result=2*(FE_local_basis_1D_reference(x_hat,basis_type,basis_index,1))./(J_11); 
end
end