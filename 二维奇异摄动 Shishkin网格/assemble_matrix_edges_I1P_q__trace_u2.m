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


function A=assemble_matrix_edges_I1P_q__trace_u2(coe_fun,P,T,Tb_q_trial,Eb,matrixsize2,basis_type_q_trial,basis_type_trace_u_test,para,geo)

order_Gauss=para.order_Gauss;   %��˹���ֵĽ���
% alpha_coef=para.alpha_coef;    %�ͷ�������ϵ��
%===����һ��ָ����С��ϡ����� A��
%����ʼ���������洢Ԫ������������Ͳ��Կռ�ľֲ�������������
%�Լ����ڷ��ߺ��������μ��������

number_of_elements=size(T,2);      %��Ԫ����
number_of_local_basis_trial=size(Tb_q_trial,1);     %�ڵ�����Ԫ��trial�ֲ��������ĸ���
if basis_type_trace_u_test==104
    number_of_local_basis_test=5;
elseif basis_type_trace_u_test==103
   number_of_local_basis_test=4;
elseif basis_type_trace_u_test==102
    number_of_local_basis_test=3;
elseif basis_type_trace_u_test==101
    number_of_local_basis_test=2;
end
A=sparse(2*number_of_local_basis_test*geo.N*(geo.N-1),matrixsize2);    
% number_of_local_basis_test=size(Tb_test,1);    %�ڵ�����Ԫ��test�ֲ��������ĸ���

%===ѭ������װ�����α��ϵĹ��ס�
for n=1:number_of_elements   % Loop over the triangles.
    for j=1:4   %������ǰ�����εıߡ�����ÿ���ߣ���ȷ����Щ�����Σ�T_left, T_right����˱�����
    current_edge=T(number_of_local_basis_trial+j,n);     %T��ÿ����Ԫ�ϵĵ�5��6��7�м�¼��õ�Ԫ�����ߵı������
    T_left=n;       %��Ӧ�����ǵ���n   %���ﶨ��õ�Ԫ�����б����ĵ�Ԫ�������Ǹõ�Ԫ����
    if (Eb(3,current_edge)==n)   
        T_right=Eb(4,current_edge);
%     else
%         T_right=Eb(3,current_edge);
%     end

    
    
    %��ȡ��ǰ�ߵĶ��㲢����ߵĳ��Ⱥͷ�������
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(2,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_left=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���
    %% 
    %% 
%     edge_length=norm(P(:,vertex1)-P(:,vertex2),2);    %�����ߵĳ��ȼ���
    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %�����˸����������һ����������ÿ�����ϵĲ�û�м��㣩

if normal(1)==0
    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %ʹ��generate_Gauss_1DΪ�߻��ֵ����ɸ�˹�ڵ��Ȩ��
   
    vertices_ww=[P(:,vertex1),P(:,vertex2)];

    if(T_right~=-1)    %�����Ƿ��ڱ߽���
         vertices_right=P(:,T((1:4),T_right));   %���Ҳ൥Ԫ�Ķ�������
       for alpha=1:number_of_local_basis_trial
           for beta=1:number_of_local_basis_test 
  
                
                
                %====== Contribution of terms -\int \mu [[\phi_j]]{{\nabla phi_i}}
                %��ͬһ�������ϵĹ��� \phi_j^+\phi_i^+
                int_value=Gauss_quadrature_2D_volume_trial_test_q_trace_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_q_trial,alpha,0,0,basis_type_trace_u_test,beta,0);
                A(Eb(4+beta,current_edge),Tb_q_trial(alpha,T_left))=A(Eb(4+beta,current_edge),Tb_q_trial(alpha,T_left))+int_value; % minus because we have - integral.
 
                int_value=Gauss_quadrature_2D_volume_trial_test_q_trace_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_right,basis_type_q_trial,alpha,0,0,basis_type_trace_u_test,beta,0);
                A(Eb(4+beta,current_edge),Tb_q_trial(alpha,T_right))=A(Eb(4+beta,current_edge),Tb_q_trial(alpha,T_right))-int_value; % minus because we have - integral.
 
           end
        end
        
    end
end
    end
    end
end
end
