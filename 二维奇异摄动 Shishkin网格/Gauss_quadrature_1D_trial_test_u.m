function int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_u_trial,basis_index_trial,der_x_trial,der_y_trial,basis_type_u_test,basis_index_test,der_x_test,der_y_test)
%Gauss_quadrature_1D_trial_test使用高斯正交计算积分。
%它接收高斯节点和权值、单元的顶点、基本类型、基索引和有限元空间的导数阶数，用于试验和测试有限元空间。
%顶点1表示trial基的顶点，顶点2表示test三角形的顶点。

% int_value=0;   %初始化
Gpn=size(Gauss_nodes,2);   %Gauss节点的个数
r=0;
for k=1:Gpn     %遍历所有高斯节点
        %将参考基函数映射到参考单元上
        r=r+Gauss_weights(1,k)*feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices1,basis_type_u_trial,basis_index_trial,der_x_trial,der_y_trial)*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices2,basis_type_u_test,basis_index_test,der_x_test,der_y_test);
end
int_value=r;
end