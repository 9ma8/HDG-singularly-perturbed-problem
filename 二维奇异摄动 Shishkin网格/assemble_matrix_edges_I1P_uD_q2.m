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


function A=assemble_matrix_edges_I1P_uD_q2(coe_fun2,P,T,Tb_q_test,Eb,matrixsize1,basis_type_q_test,para)

order_Gauss=para.order_Gauss;   %��˹���ֵĽ���
% alpha_coef=para.alpha_coef;    %�ͷ�������ϵ��
%===����һ��ָ����С��ϡ����� A��
%����ʼ���������洢Ԫ������������Ͳ��Կռ�ľֲ�������������
%�Լ����ڷ��ߺ��������μ��������
A=sparse(matrixsize1,1);    
number_of_elements=size(T,2);      %��Ԫ����
% number_of_local_basis_trial=size(Tb_trial,1);     %�ڵ�����Ԫ��trial�ֲ��������ĸ���
number_of_local_basis_test=size(Tb_q_test,1);    %�ڵ�����Ԫ��test�ֲ��������ĸ���

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
  if(T_right==-1)    %�����Ƿ��ڱ߽���
    
    %��ȡ��ǰ�ߵĶ��㲢����ߵĳ��Ⱥͷ�������
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(2,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_left=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���
    %% 
    %% 

    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %�����˸����������һ����������ÿ�����ϵĲ�û�м��㣩

    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %ʹ��generate_Gauss_1DΪ�߻��ֵ����ɸ�˹�ڵ��Ȩ��
%    vertices_ww=[P(:,vertex1),P(:,vertex2)];
    if normal(2)==0  

        for beta=1:number_of_local_basis_test 
   
                %===== Contribution of terms +\int \mu [[phi_j]][[\phi_i]]
                %��ͬһ�������ϵĹ��� \phi_j^+\phi_i^+
                int_value=normal(1)*Gauss_quadrature_2D_volume_test_uD_q1(coe_fun2,Gauss_nodes,Gauss_weights,vertices_left,basis_type_q_test,beta,0,0);
                A(Tb_q_test(beta,T_left),1)=A(Tb_q_test(beta,T_left),1) +int_value; % ͬһ��������Ĺ���
    
         end


    end
  end
    end
end
end