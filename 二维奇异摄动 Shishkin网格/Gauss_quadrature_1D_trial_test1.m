function int_value=Gauss_quadrature_1D_trial_test1(normal,Gauss_nodes,Gauss_weights,vertices1,vertices2,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial,basis_type_test,basis_index_test,der_x_test,der_y_test,data)
%Gauss_quadrature_1D_trial_testʹ�ø�˹����������֡�
%�����ո�˹�ڵ��Ȩֵ����Ԫ�Ķ��㡢�������͡�������������Ԫ�ռ�ĵ�����������������Ͳ�������Ԫ�ռ䡣
%����1��ʾtrial���Ķ��㣬����2��ʾtest�����εĶ��㡣

% int_value=0;   %��ʼ��
Gpn=size(Gauss_nodes,2);   %Gauss�ڵ�ĸ���
int_value=0;
for k=1:Gpn     %�������и�˹�ڵ�
        %���ο�������ӳ�䵽�ο���Ԫ��
        int_value=int_value+Gauss_weights(1,k)*(normal(1)*feval(data.b1,Gauss_nodes(1,k),Gauss_nodes(2,k))+feval(data.b2,Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices1,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices2,basis_type_test,basis_index_test,der_x_test,der_y_test);
end
% int_value=r;
end