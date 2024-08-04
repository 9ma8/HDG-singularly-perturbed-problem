%-------------------------------------------------------------------------
% assemble_matrix_edges_IP����
% coe_fun:��ɢϵ��
% �������:T,P
% Tb_trial,Tb_test:����trial��test���޵�Ԫ�ռ���Ϣ�ľ���
% matrixsize1,matrixsize2:�նȾ���Ĵ�С
% basis_type_trial,basis_type_test:����Ԫ�ռ���������/����
% para: ������˹����������DG�ͷ�ϵ���Ľṹ
% assemble_matrix_edges_IP:
% A(����)��������IP��ɢ�������б���������ĸնȾ���

%-------------------------------------------------------------------------


function A=assemble_matrix_edges_I1P_uu(coe_fun,P,T,Tb_trial,Tb_test,Eb,matrixsize1,matrixsize2,basis_type_u_trial,basis_type_u_test,para,tau)

order_Gauss=para.order_Gauss;   %��˹���ֵĽ���
% alpha_coef=para.alpha_coef;    %�ͷ�������ϵ��
%===����һ��ָ����С��ϡ����� A��
%����ʼ���������洢Ԫ������������Ͳ��Կռ�ľֲ�������������
%�Լ����ڷ��ߺ��������μ��������
A=sparse(matrixsize1,matrixsize2);    
number_of_elements=size(T,2);      %��Ԫ����
number_of_local_basis_trial=size(Tb_trial,1);     %�ڵ�����Ԫ��trial�ֲ��������ĸ���
number_of_local_basis_test=size(Tb_test,1);    %�ڵ�����Ԫ��test�ֲ��������ĸ���

%===ѭ������װ�����α��ϵĹ��ס�
for n=1:number_of_elements   % Loop over the triangles.
    for j=1:4   %������ǰ�����εıߡ�����ÿ���ߣ���ȷ����Щ�����Σ�T_left, T_right����˱�����
    current_edge=T(number_of_local_basis_test+j,n);     %T��ÿ����Ԫ�ϵĵ�5��6��7�м�¼��õ�Ԫ�����ߵı������
    T_left=n;       %��Ӧ�����ǵ���n   %���ﶨ��õ�Ԫ�����б����ĵ�Ԫ�������Ǹõ�Ԫ����
    if (Eb(3,current_edge)==n)   
        T_right=Eb(4,current_edge);
    else
        T_right=Eb(3,current_edge);
    end
    
    
    %��ȡ��ǰ�ߵĶ��㲢����ߵĳ��Ⱥͷ�������
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(2,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_left=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���

    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %ʹ��generate_Gauss_1DΪ�߻��ֵ����ɸ�˹�ڵ��Ȩ��
   
%     vertices_right=P(:,T((1:4),T_right));   %���Ҳ൥Ԫ�Ķ�������
%     ���ڷǱ߽�ߣ�T_right !=-1������������������Ժ�����ռ�ľֲ����������Լ���Ծ���A�Ĺ��ס�
%     ���漰��������ϵĻ��֣�ͨ����Ҫ�ڸ�˹���������������ֵ�����ǵ��ݶ��Լ�ϵ��������ֵ
    if(T_right~=-1)    %�����Ƿ��ڱ߽���
         vertices_right=P(:,T((1:4),T_right));   %���Ҳ൥Ԫ�Ķ�������
       for alpha=1:number_of_local_basis_trial 
          for beta=1:number_of_local_basis_test 
          
                %===== Contribution of terms +\int \mu [[phi_j]][[\phi_i]]
                %��ͬһ�������ϵĹ��� \phi_j^+\phi_i^+
                int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_u_trial,alpha,0,0,basis_type_u_test,beta,0,0);
                A(Tb_test(beta,T_left),Tb_trial(alpha,T_left))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_left)) +tau*int_value; % ͬһ��������Ĺ���
%                 int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_right,vertices_left,basis_type_u_trial,alpha,0,0,basis_type_u_test,beta,0,0);
%                 A(Tb_test(beta,T_left),Tb_trial(alpha,T_right))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_right)) +tau*int_value; % ͬһ��������Ĺ���

           end
        end
        
        %�� T_right ���� -1 ʱ�����ʾ�߽��������Ҳ�û�����������εı�
        %�����Բ�ͬ�ķ�ʽ���������ıߵĹ��ף�רע��ֻ�漰��൥Ԫ����
        
    elseif(T_right==-1)
         
         for alpha=1:number_of_local_basis_trial
            for beta=1:number_of_local_basis_test
                % Contribution of terms +\int \mu [[phi_j]][[\phi_i]]
                int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_u_trial,alpha,0,0,basis_type_u_test,beta,0,0);
                A(Tb_test(beta,T_left),Tb_trial(alpha,T_left))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_left)) +tau*int_value; % Contribution from terms of the same triangle.

           end
        end
    end
%     end
    end
end
end