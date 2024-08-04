%-------------------------------------------------------------------------
% assemble_rhs_2D����
% f_fun:   force term
% P,T:   �������;
% Tb_test:   ����trial��test FEM�ռ���Ϣ�ľ���
% matrixsize2:  ����FEM�ռ�Ĵ�С
% basis_type_test:   test����Ԫ�ռ�����
% der_test_x,der_test_y:  test���������Ľ���
% assemble_rhs_1D����:
% b (matrix) whose entry b_{i}= \int (f_fun)((\partial_x)^der_test_x (\partial_x)^der_test_y \phi_i) 


%-------------------------------------------------------------------------


function b=assemble_rhs_2D(f_fun,P,T,Tb_test,matrixsize2,basis_type_test,der_x_test,der_y_test,order_Gauss)
b=zeros(matrixsize2,1);
number_of_elements=size(T,2);
number_of_local_basis_test=size(Tb_test,1);
for n=1:number_of_elements     %�������еĵ�Ԫ
    vertices=P(:,T((1:4),n));      %��Ԫ�Ķ���
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss);     %����Gauss����ڵ����Ȩ��
    J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);    %�ſɱȷ���任
    for beta=1:number_of_local_basis_test     
        b(Tb_test(beta,n),1)=b(Tb_test(beta,n),1)+J*Gauss_quadrature_2D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,der_x_test,der_y_test);
    end
end






