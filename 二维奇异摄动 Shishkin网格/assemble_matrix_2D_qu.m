%-------------------------------------------------------------------------
% assemble_matrix_2D����
% coe_fun:��ɢϵ��
%�������:T,P
% Tb_trial,Tb_test:����trial��test���޵�Ԫ�ռ���Ϣ�ľ���
% matrixsize1,matrixsize2:�նȾ���Ĵ�С
% basis_type_trial,basis_type_test:����Ԫ�ռ���������/����
% der_x_trial,der_x_test, der_y_trial,der_y_test: bilinea��ʽ�ĵ����Ľ���
% order_gauss:��˹����������
% assemble_matrix_2D returns: 
% A (matrix) whose entry A_{i,j}= \int (c(\partial_x)^der_x_trial (\partial_y)^der_y_trial \phi_j)((\partial_x)^der_x_test (\partial_y)^der_y_test \phi_j)

%��δ���ĺ��Ĳ�����ͨ����˹������������ÿ���������ϵĻ���ֵ�������ݼ�������װ�նȾ���

%-------------------------------------------------------------------------

function A=assemble_matrix_2D_qu(coe_fun,P,T,Tb_q_trial,Tb_test,matrixsize1,matrixsize2,basis_type_q_trial,der_x_trial,der_y_trial,basis_type_u_test,der_x_test,der_y_test,order_Gauss)
%===����ϡ����󲢶������
number_of_elements=size(T,2); %��Ԫ�ĸ���
number_of_local_basis_trial=size(Tb_q_trial,1);    %�ڵ�����Ԫ��trial�ֲ��������Ĵ���
number_of_local_basis_test=size(Tb_test,1);   %��������
A=sparse(matrixsize1,matrixsize2);

%=== ѭ������װ
for n=1:number_of_elements    %�������еĵ�Ԫ
    vertices=P(:,T((1:4),n));    %��ǰ��Ԫ�Ķ��㡣��Ԫn������T(:,n)��λ��P(:,T(:,n))��
    J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);    %�������Ա任���ſɱȾ���
   
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss);    %�ڵ�Ԫ�ϴ�����˹�ڵ㡣1 .����  2 .�������ȡ�
    for alpha=1:number_of_local_basis_trial    %ѭ��������Ԫ�ϵ�trial�����Ի��������Ԫ�ռ䡣
        for beta=1:number_of_local_basis_test      
%             int_value=0;
            %���ø�˹���ַ��� Gauss_quadrature_2D_volume_trial_test �������ڵ�ǰ�������ϻ������˻��Ļ���ֵ
            int_value=J*Gauss_quadrature_2D_volume_trial_test_qu(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_q_trial,alpha,der_x_trial,der_y_trial,basis_type_u_test,beta,der_x_test,der_y_test);
            A(Tb_test(beta,n),Tb_q_trial(alpha,n))=A(Tb_test(beta,n),Tb_q_trial(alpha,n))+int_value;

        end
    end
    
end
% A = sparse (Ivec(1:count),Jvec(1:count),Xvec(1:count),matrixsize1,matrixsize2);  %���� sparse �������� Ivec��Jvec��Xvec �е����ݹ���ϡ����� A
end