function result=Gauss_quadrature_2D_volume_test_uD_v(g_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_u_test,basis_index_test,der_x_test,der_y_test)
Gpn=size(Gauss_nodes,2);    % number of Gauss points.
r=0;
for k=1:Gpn  %loop on the Gaussian nodes
         %将参考基函数映射到小区间上
%          c(k)=g_fun(Gauss_nodes(1,k),Gauss_nodes(2,k));
        r=r+Gauss_weights(1,k)*g_fun(Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_u_test,basis_index_test,der_x_test,der_y_test);
end
result=r;
end