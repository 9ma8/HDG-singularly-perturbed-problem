%-------------------------------------------------------------------------
% compute_errors�������ڼ�������Ԫ���l2����H1���
%������:
% u:  ��ɢ����Ľ�
% uexvec:  ��ȷ������
% P��T, Tb, Eb:  ������Ϣ��FE��Ϣ��
% para:�����ṹ�Ľṹ������������������basis_type(FEM/DG�ռ�����ͣ������ԡ����ε�)����˹����������(���յ�)DG�ͷ�ϵ����
% method:ָ��������ʹ��FEM����DG���ַ��������method=='DG'����ô��Ҳ������ڱ��ϵĹ��ס�


% ע��:ֻ�б߽��ϵ�uex=0�����Ϊ0��, DG����������ȷ��
% ��Ϊ����ĿǰΪֹ������[[uex]]���ڲ�����=0�����ڱ߽��ϲ�����[[uex]]=0������uex=0


%-------------------------------------------------------------------------


function [error_L2,error_Lq,error_DG]=compute_errors_HDG(u,q1,q2,u1,uex,pex,qex,para,method,geo,tau,data,varepsilon)

error_L2=0;
error_L3=0;
error_Lq1=0;
error_Lq2=0;

[P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);

number_of_elements=size(T,2); 
% number_of_edges=size(Eb,2);   %����߽��ϵıߵ�����
number_of_local_basis_trial=size(Tb,1);    %����ֲ�������������
for n=1:number_of_elements % Loop over the triangles.
    vertices=P(:,T((1:4),n)); %��ȡ��ǰ��Ԫ�Ķ�������
    J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);% ���㵱ǰ��Ԫ���ſɱ�����ʽ
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,para.order_Gauss); %����2D��˹���ֵ��Ȩ��
    u_app_n=u(Tb((1:number_of_local_basis_trial),n));     % ��ȡ��ǰ��Ԫ�ϵ���ֵ��
 
    q1_app_n=q1(Tb((1:number_of_local_basis_trial),n));     
    
    q2_app_n=q2(Tb((1:number_of_local_basis_trial),n));     

    for k=1:size(Gauss_nodes,1) % ѭ������ÿ����˹���ֵ�
        %��ʼ���������ı���Ϊ��
        err_element=0;   
%         err_element1=0; 
        err_element_q1=0;    % construction derivative_x error inside each element
        err_element_q2=0;   % construction derivative_y error inside each element
        for alpha=1:number_of_local_basis_trial    %ѭ������ÿ���ֲ�������
            
            %�������
            err_element=err_element + u_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_u,alpha,0,0);

            err_element_q1=err_element_q1 +  q1_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_q,alpha,0,0);
            err_element_q2=err_element_q2 + q2_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_q,alpha,0,0);
           
          
        end
        %����L2������H1�������
         error_L2=error_L2+J*Gauss_weights(k)*(uex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element)^2;
        error_L3=error_L3+J*Gauss_weights(k)*(1.5*(1+Gauss_nodes(k,2)^2))*(uex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element)^2;
        error_Lq1=error_Lq1+J*Gauss_weights(k)*(pex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element_q1)^2;
        error_Lq2=error_Lq2+J*Gauss_weights(k)*(qex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element_q2)^2;
    end

end


if(strcmp(method,'DG')==1)   %����Ƿ�ʹ��DG����
    
 if para.basis_type_trace_u==104
    number_of_local_basis_testt=5;
elseif para.basis_type_trace_u==103
   number_of_local_basis_testt=4;
elseif para.basis_type_trace_u==102
    number_of_local_basis_testt=3;
elseif para.basis_type_trace_u==101
    number_of_local_basis_testt=2;
end
    
    
error_Edges2=0;   %��ʼ��DG���Ϊ��
for n=1:number_of_elements   % Loop over the triangles.
  for j=1:4   
    current_edge=T(number_of_local_basis_trial+j,n);     %T�м�¼��õ�Ԫ�����ߵı������
    T_left=n;                                                     %��Ӧ�����ǵ���n   %���ﶨ��õ�Ԫ�����б����ĵ�Ԫ�������Ǹõ�Ԫ����
     if (Eb(number_of_local_basis_testt+1,current_edge)==n)   
        T_right=Eb(number_of_local_basis_testt+2,current_edge);
    else
        T_right=Eb(number_of_local_basis_testt+1,current_edge);

     end
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(number_of_local_basis_testt,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_plus=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���
    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %�����˸����������һ����������ÿ�����ϵĲ�û�м��㣩
     [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(para.order_Gauss);
     vertices_ww1=[P(:,vertex1)';P(:,vertex2)'];
     [Gauss_nodes,Gauss_weights]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_ww1);
%     [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),para.order_Gauss);    %ʹ��generate_Gauss_1DΪ�߻��ֵ����ɸ�˹�ڵ��Ȩ��

          uapp_plus=u1(Eb((number_of_local_basis_testt+3:end),current_edge));    %��ȡ����������ϵ���ֵ��

          u_app_plus1=u(Tb((1:number_of_local_basis_trial),T_left));    %��ȡ����������ϵ���ֵ��
 
        for k=1:size(Gauss_nodes,1) % ѭ������ÿ����˹���ֵ�
            vertices1=[P(:,vertex1),P(:,vertex2)];
            %��ʼ��DG���ı���Ϊ��
            int_square_plus2=0;  % term \int_e (uapp+)^2
            int_square_minus2=0;   % term \int_e (uapp-)^2
            int_square_minus3=0;
            uapp_plus2=0;
            uapp_plus3=0;

            for alpha=1:number_of_local_basis_trial     %ѭ������ÿ���ֲ�������
                 
                uapp_plus2=uapp_plus2 +u_app_plus1(alpha)*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices_plus,para.basis_type_u,alpha,0,0);
            end
            for beta=1:number_of_local_basis_testt
                if Gauss_nodes(1,1)-Gauss_nodes(1,2)~=0

                  uapp_plus3=uapp_plus3 +uapp_plus(beta)*FE_local_basis_1D_local(Gauss_nodes(1,k),vertices1(1,:),para.basis_type_trace_u,beta,0);
%                  int_exact2= int_exact2+(u_app_plus(alpha)-u_exact_nm(alpha))*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices_plus,para.basis_type,alpha,0,0);
               elseif Gauss_nodes(1,1)-Gauss_nodes(1,2)==0
                   uapp_plus3=uapp_plus3 +uapp_plus(beta)*FE_local_basis_1D_local(Gauss_nodes(2,k),vertices1(2,:),para.basis_type_trace_u,beta,0);
               end
            end    
        end
            int_square_plus2=int_square_plus2 +Gauss_weights(k)*(tau-0.5*data.b1(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(1)-0.5*data.b2(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*(uapp_plus2).^2;
%             �����Ҳ��������ϵĻ�����
            int_square_minus2=int_square_minus2 +Gauss_weights(k)*(tau-0.5*data.b1(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(1)-0.5*data.b2(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*(uapp_plus3).^2;
           int_square_minus3=int_square_minus3 +Gauss_weights(k)*(tau-0.5*data.b1(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(1)-0.5*data.b2(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*(uapp_plus3)*(uapp_plus2);


        %����DG���
        error_Edges2=error_Edges2 + (int_square_plus2+int_square_minus2-2*int_square_minus3);  
  end  
 end




%�����ܵ�DG�������߽��ϵ�DG���
error_DG=(varepsilon^(-1))*(error_Lq1+error_Lq2)+error_L3+error_Edges2;
error_DG=sqrt(error_DG);
end
% ����L2��������ƽ�������õ����յ�L2�������ֵ
error_L2=sqrt(error_L2);
%����H1��������ƽ�������õ����յ�H1�������ֵ
error_Lq=sqrt(error_Lq1+error_Lq2);

end







