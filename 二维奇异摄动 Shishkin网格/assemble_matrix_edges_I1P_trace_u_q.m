%-------------------------------------------------------------------------
% assemble_matrix_edges_IP接收
% coe_fun:扩散系数
% 网格矩阵:T,P
% Tb_trial,Tb_test:包含trial和test有限单元空间信息的矩阵
% matrixsize1,matrixsize2:刚度矩阵的大小
% basis_type_trial,basis_type_test:有限元空间类型试验/测试
% para: 包含高斯正交阶数和DG惩罚系数的结构
% assemble_matrix_edges_IP:
% A(矩阵)，它是与IP离散化的所有边项相关联的刚度矩阵。

%-------------------------------------------------------------------------


function A=assemble_matrix_edges_I1P_trace_u_q(coe_fun,P,T,Tb_q_test,Eb,matrixsize1,basis_type_trace_u_trial,basis_type_q_test,para,geo)

order_Gauss=para.order_Gauss;   %高斯积分的阶数
% alpha_coef=para.alpha_coef;    %惩罚参数的系数
%===创建一个指定大小的稀疏矩阵 A，
%并初始化变量来存储元素数量、试验和测试空间的局部基函数数量，
%以及用于法线和其他几何计算的向量

    number_of_elements=size(T,2);      %单元个数
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
% number_of_local_basis_trial=size(Tb_trial,1);     %在单个单元上trial局部基函数的个数
number_of_local_basis_test=size(Tb_q_test,1);    %在单个单元上test局部基函数的个数

%===循环以组装三角形边上的贡献。
for n=1:number_of_elements   % Loop over the triangles.
    for j=1:4   %遍历当前三角形的边。对于每个边，它确定哪些三角形（T_left, T_right）与此边相邻
    current_edge=T(number_of_local_basis_test+j,n);     %T中每个单元上的第5，6，7行记录这该单元三条边的编号索引
    T_left=n;       %它应该总是等于n   %这里定义该单元上所有边左侧的单元索引都是该单元本身
    if (Eb(3,current_edge)==n)   
        T_right=Eb(4,current_edge);
%     else
%         T_right=Eb(3,current_edge);
%     end
%     
    
    %获取当前边的顶点并计算边的长度和法线向量
    vertex1=Eb(1,current_edge);   %边的第一个顶点索引
    vertex2=Eb(2,current_edge);    %边的第二个顶点索引
    vertex_left_index=T((1:4),T_left);   %该边左侧单元的顶点索引
    vertices_left=P(:,vertex_left_index);     % 在边的左侧单元的顶点
    %% 
    %% 

    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %计算了该三角形面的一个法向量（每条边上的并没有计算）

if normal(1)~=0
    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %使用generate_Gauss_1D为边积分的生成高斯节点和权重
   vertices_ww=[P(:,vertex1),P(:,vertex2)];
    

    if(T_right~=-1)    %检查边是否在边界上
         vertices_right=P(:,T((1:4),T_right));   %边右侧单元的顶点坐标
         for alpha=1:number_of_local_basis_trial
           for beta=1:number_of_local_basis_test 
   
                int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_q(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u_trial,alpha,0,basis_type_q_test,beta,0,0);
                A(Tb_q_test(beta,T_left),Eb(4+alpha,current_edge))=A(Tb_q_test(beta,T_left),Eb(4+alpha,current_edge))+int_value; % minus because we have - integral.
             int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_q(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_right,basis_type_trace_u_trial,alpha,0,basis_type_q_test,beta,0,0);
                A(Tb_q_test(beta,T_right),Eb(4+alpha,current_edge))=A(Tb_q_test(beta,T_right),Eb(4+alpha,current_edge))-int_value; % minus because we have - integral.
               
           end
        end
        
    end
end
    end
   end
end
end
