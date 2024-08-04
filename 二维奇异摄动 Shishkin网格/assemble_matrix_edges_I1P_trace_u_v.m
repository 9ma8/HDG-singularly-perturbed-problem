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


function A=assemble_matrix_edges_I1P_trace_u_v(coe_fun,P,T,Tb_test,Eb,matrixsize1,basis_type_trace_u_trial,basis_type_u_test,para,tau,geo)

order_Gauss=para.order_Gauss;   %��˹���ֵĽ���
% alpha_coef=para.alpha_coef;    %�ͷ�������ϵ��
%===����һ��ָ����С��ϡ����� A��
%����ʼ���������洢Ԫ������������Ͳ��Կռ�ľֲ�������������
%�Լ����ڷ��ߺ��������μ��������
   
number_of_elements=size(T,2);      %��Ԫ����

if basis_type_trace_u_trial==104
    number_of_local_basis_trial=5;
elseif basis_type_trace_u_trial==103
   number_of_local_basis_trial=4;
elseif basis_type_trace_u_trial==102
    number_of_local_basis_trial=3;
elseif basis_type_trace_u_trial==101
    number_of_local_basis_trial=2;
end
A=sparse(matrixsize1,2*number_of_local_basis_trial*geo.N*(geo.N-1)); 
% number_of_local_basis_trial=size(Tb_uu_trial,1);     %�ڵ�����Ԫ��trial�ֲ��������ĸ���
number_of_local_basis_test=size(Tb_test,1);    %�ڵ�����Ԫ��test�ֲ��������ĸ���

%===ѭ������װ�����α��ϵĹ��ס�
for n=1:number_of_elements   % Loop over the triangles.
    for j=1:4   %������ǰ�����εıߡ�����ÿ���ߣ���ȷ����Щ�����Σ�T_left, T_right����˱�����
    current_edge=T(number_of_local_basis_test+j,n);     %T��ÿ����Ԫ�ϵĵ�5��6��7�м�¼��õ�Ԫ�����ߵı������
    T_left=n;       %��Ӧ�����ǵ���n   %���ﶨ��õ�Ԫ�����б����ĵ�Ԫ�������Ǹõ�Ԫ����
    if (Eb(3,current_edge)==n)   
        T_right=Eb(4,current_edge);
%     else
%         T_right=Eb(3,current_edge);
%     end
    if(T_right~=-1)    %�����Ƿ��ڱ߽���
    
    %��ȡ��ǰ�ߵĶ��㲢����ߵĳ��Ⱥͷ�������
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(2,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_left=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���
    vertices_right=P(:,T((1:4),T_right));     % �ڱߵ���൥Ԫ�Ķ���
    %% 

    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %ʹ��generate_Gauss_1DΪ�߻��ֵ����ɸ�˹�ڵ��Ȩ��
%    if normal(1)~=0
    vertices_ww=[P(:,vertex1),P(:,vertex2)];
% vertices1=P(:,vertex1);
% vertices2=P(:,vertex2);
%      vertices_ww = sort(vertices_www, 'ascend'); 
        for alpha=1:number_of_local_basis_trial 
            for beta=1:number_of_local_basis_test 
          
                %===== Contribution of terms +\int \mu [[phi_j]][[\phi_i]]
                %��ͬһ�������ϵĹ��� \phi_j^+\phi_i^+
                %Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_left,vertices1,vertices2,basis_type_u_test,beta,0,0,basis_type_trace_u_trial,alpha,0);
                %Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u_trial,alpha,0,basis_type_u_test,beta,0,0);
                int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u_trial,alpha,0,basis_type_u_test,beta,0,0);
                A(Tb_test(beta,T_left),Eb(4+alpha,current_edge))=A(Tb_test(beta,T_left),Eb(4+alpha,current_edge)) +tau*int_value; % ͬһ��������Ĺ���
                   int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_right,basis_type_trace_u_trial,alpha,0,basis_type_u_test,beta,0,0);
                A(Tb_test(beta,T_right),Eb(4+alpha,current_edge))=A(Tb_test(beta,T_right),Eb(4+alpha,current_edge)) +tau*int_value; % ͬһ��������Ĺ���
   %Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_right,basis_type_trace_u_trial,alpha,0,basis_type_u_test,beta,0,0);
           %Gauss_quadrature_1D_trial_test_coupling(coe_fun,Gauss_nodes,Gauss_weights,vertices_right,vertices1,vertices2,basis_type_u_test,beta,0,0,basis_type_trace_u_trial,alpha,0);
            end
        end

    end
    end
    end
end
end