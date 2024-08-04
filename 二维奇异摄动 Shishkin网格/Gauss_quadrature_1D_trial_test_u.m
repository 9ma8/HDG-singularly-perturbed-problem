function int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_u_trial,basis_index_trial,der_x_trial,der_y_trial,basis_type_u_test,basis_index_test,der_x_test,der_y_test)
%Gauss_quadrature_1D_trial_testʹ�ø�˹����������֡�
%�����ո�˹�ڵ��Ȩֵ����Ԫ�Ķ��㡢�������͡�������������Ԫ�ռ�ĵ�����������������Ͳ�������Ԫ�ռ䡣
%����1��ʾtrial���Ķ��㣬����2��ʾtest�����εĶ��㡣

% int_value=0;   %��ʼ��
Gpn=size(Gauss_nodes,2);   %Gauss�ڵ�ĸ���
r=0;
for k=1:Gpn     %�������и�˹�ڵ�
        %���ο�������ӳ�䵽�ο���Ԫ��
        r=r+Gauss_weights(1,k)*feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices1,basis_type_u_trial,basis_index_trial,der_x_trial,der_y_trial)*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices2,basis_type_u_test,basis_index_test,der_x_test,der_y_test);
end
int_value=r;
end