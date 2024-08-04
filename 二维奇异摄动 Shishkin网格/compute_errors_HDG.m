%-------------------------------------------------------------------------
% compute_errors函数用于计算有限元解的l2误差和H1误差
%它接收:
% u:  离散问题的解
% uexvec:  精确解向量
% P和T, Tb, Eb:  网格信息和FE信息。
% para:包含结构的结构，包含几个参数，如basis_type(FEM/DG空间的类型，如线性、二次等)、高斯正交阶数和(最终的)DG惩罚系数。
% method:指定我们是使用FEM还是DG的字符串。如果method=='DG'，那么它也会计算在边上的贡献。


% 注意:只有边界上的uex=0（真解为0）, DG范数才算正确。
% 因为，到目前为止，假设[[uex]]在内部边上=0，而在边界上不假设[[uex]]=0，除非uex=0


%-------------------------------------------------------------------------


function [error_L2,error_Lq,error_DG]=compute_errors_HDG(u,q1,q2,u1,uex,pex,qex,para,method,geo,tau,data,varepsilon)

error_L2=0;
error_L3=0;
error_Lq1=0;
error_Lq2=0;

[P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);

number_of_elements=size(T,2); 
% number_of_edges=size(Eb,2);   %计算边界上的边的数量
number_of_local_basis_trial=size(Tb,1);    %计算局部基函数的数量
for n=1:number_of_elements % Loop over the triangles.
    vertices=P(:,T((1:4),n)); %获取当前单元的顶点坐标
    J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);% 计算当前单元的雅可比行列式
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,para.order_Gauss); %生成2D高斯积分点和权重
    u_app_n=u(Tb((1:number_of_local_basis_trial),n));     % 获取当前单元上的数值解
 
    q1_app_n=q1(Tb((1:number_of_local_basis_trial),n));     
    
    q2_app_n=q2(Tb((1:number_of_local_basis_trial),n));     

    for k=1:size(Gauss_nodes,1) % 循环遍历每个高斯积分点
        %初始化构造误差的变量为零
        err_element=0;   
%         err_element1=0; 
        err_element_q1=0;    % construction derivative_x error inside each element
        err_element_q2=0;   % construction derivative_y error inside each element
        for alpha=1:number_of_local_basis_trial    %循环遍历每个局部基函数
            
            %计算误差
            err_element=err_element + u_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_u,alpha,0,0);

            err_element_q1=err_element_q1 +  q1_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_q,alpha,0,0);
            err_element_q2=err_element_q2 + q2_app_n(alpha)*FE_local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,para.basis_type_q,alpha,0,0);
           
          
        end
        %更新L2范数和H1范数误差
         error_L2=error_L2+J*Gauss_weights(k)*(uex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element)^2;
        error_L3=error_L3+J*Gauss_weights(k)*(1.5*(1+Gauss_nodes(k,2)^2))*(uex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element)^2;
        error_Lq1=error_Lq1+J*Gauss_weights(k)*(pex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element_q1)^2;
        error_Lq2=error_Lq2+J*Gauss_weights(k)*(qex(Gauss_nodes(k,1),Gauss_nodes(k,2),varepsilon)-err_element_q2)^2;
    end

end


if(strcmp(method,'DG')==1)   %检查是否使用DG方法
    
 if para.basis_type_trace_u==104
    number_of_local_basis_testt=5;
elseif para.basis_type_trace_u==103
   number_of_local_basis_testt=4;
elseif para.basis_type_trace_u==102
    number_of_local_basis_testt=3;
elseif para.basis_type_trace_u==101
    number_of_local_basis_testt=2;
end
    
    
error_Edges2=0;   %初始化DG误差为零
for n=1:number_of_elements   % Loop over the triangles.
  for j=1:4   
    current_edge=T(number_of_local_basis_trial+j,n);     %T中记录这该单元三条边的编号索引
    T_left=n;                                                     %它应该总是等于n   %这里定义该单元上所有边左侧的单元索引都是该单元本身
     if (Eb(number_of_local_basis_testt+1,current_edge)==n)   
        T_right=Eb(number_of_local_basis_testt+2,current_edge);
    else
        T_right=Eb(number_of_local_basis_testt+1,current_edge);

     end
    vertex1=Eb(1,current_edge);   %边的第一个顶点索引
    vertex2=Eb(number_of_local_basis_testt,current_edge);    %边的第二个顶点索引
    vertex_left_index=T((1:4),T_left);   %该边左侧单元的顶点索引
    vertices_plus=P(:,vertex_left_index);     % 在边的左侧单元的顶点
    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %计算了该三角形面的一个法向量（每条边上的并没有计算）
     [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(para.order_Gauss);
     vertices_ww1=[P(:,vertex1)';P(:,vertex2)'];
     [Gauss_nodes,Gauss_weights]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_ww1);
%     [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),para.order_Gauss);    %使用generate_Gauss_1D为边积分的生成高斯节点和权重

          uapp_plus=u1(Eb((number_of_local_basis_testt+3:end),current_edge));    %获取左侧三角形上的数值迹

          u_app_plus1=u(Tb((1:number_of_local_basis_trial),T_left));    %获取左侧三角形上的数值解
 
        for k=1:size(Gauss_nodes,1) % 循环遍历每个高斯积分点
            vertices1=[P(:,vertex1),P(:,vertex2)];
            %初始化DG误差的变量为零
            int_square_plus2=0;  % term \int_e (uapp+)^2
            int_square_minus2=0;   % term \int_e (uapp-)^2
            int_square_minus3=0;
            uapp_plus2=0;
            uapp_plus3=0;

            for alpha=1:number_of_local_basis_trial     %循环遍历每个局部基函数
                 
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
%             计算右侧三角形上的积分项
            int_square_minus2=int_square_minus2 +Gauss_weights(k)*(tau-0.5*data.b1(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(1)-0.5*data.b2(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*(uapp_plus3).^2;
           int_square_minus3=int_square_minus3 +Gauss_weights(k)*(tau-0.5*data.b1(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(1)-0.5*data.b2(Gauss_nodes(1,k),Gauss_nodes(2,k))*normal(2))*(uapp_plus3)*(uapp_plus2);


        %更新DG误差
        error_Edges2=error_Edges2 + (int_square_plus2+int_square_minus2-2*int_square_minus3);  
  end  
 end




%计算总的DG误差，包括边界上的DG误差
error_DG=(varepsilon^(-1))*(error_Lq1+error_Lq2)+error_L3+error_Edges2;
error_DG=sqrt(error_DG);
end
% 计算L2范数误差的平方根，得到最终的L2范数误差值
error_L2=sqrt(error_L2);
%计算H1范数误差的平方根，得到最终的H1范数误差值
error_Lq=sqrt(error_Lq1+error_Lq2);

end







