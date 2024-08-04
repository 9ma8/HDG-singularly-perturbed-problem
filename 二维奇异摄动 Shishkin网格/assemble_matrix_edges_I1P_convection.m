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


function A=assemble_matrix_edges_I1P_convection(P,T,Tb_trial,Tb_test,Eb,matrixsize1,matrixsize2,basis_type_trial,basis_type_test,para,data)

order_Gauss=para.order_Gauss;   %高斯积分的阶数
% alpha_coef=para.alpha_coef;    %惩罚参数的系数
%===创建一个指定大小的稀疏矩阵 A，
%并初始化变量来存储元素数量、试验和测试空间的局部基函数数量，
%以及用于法线和其他几何计算的向量
A=sparse(matrixsize1,matrixsize2);    
number_of_elements=size(T,2);      %单元个数
number_of_local_basis_trial=size(Tb_trial,1);     %在单个单元上trial局部基函数的个数
number_of_local_basis_test=size(Tb_test,1);    %在单个单元上test局部基函数的个数
% normal=zeros(2,1);
% vintern=zeros(2,1);
% tau=zeros(2,1);
%===循环以组装三角形边上的贡献。
for n=1:number_of_elements   % Loop over the triangles.
 for j=1:4   %遍历当前三角形的边。对于每个边，它确定哪些三角形（T_left, T_right）与此边相邻
    current_edge=T(number_of_local_basis_test+j,n);     %T中每个单元上的第5，6，7行记录这该单元三条边的编号索引
    T_left=n;       %它应该总是等于n   %这里定义该单元上所有边左侧的单元索引都是该单元本身
    if (Eb(3,current_edge)==n)   
        T_right=Eb(4,current_edge);
    else
        T_right=Eb(3,current_edge);
    end
    
    
    %获取当前边的顶点并计算边的长度和法线向量
    vertex1=Eb(1,current_edge);   %边的第一个顶点索引
    vertex2=Eb(2,current_edge);    %边的第二个顶点索引
    vertex_left_index=T((1:4),T_left);   %该边左侧单元的顶点索引
    vertices_left=P(:,vertex_left_index);     % 在边的左侧单元的顶点
%     edge_length=norm(P(:,vertex1)-P(:,vertex2),2);    %这条边的长度计算
% if (Eb(3,current_edge)==n)   
%     normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %计算了该三角形面的一个法向量（每条边上的并没有计算）
% else
%    normal=-find_normal(vertex1,vertex2,vertex_left_index,P);
% end
    normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %计算了该三角形面的一个法向量（每条边上的并没有计算）

    [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);    %使用generate_Gauss_1D为边积分的生成高斯节点和权重
   
   if normal(1)*feval(data.b1,(P(1,vertex1)+P(1,vertex2))/2,(P(2,vertex1)+P(2,vertex2))/2)+feval(data.b2,(P(1,vertex1)+P(1,vertex2))/2,(P(2,vertex1)+P(2,vertex2))/2)*normal(2)<0 && T_right~=-1
%     对于非边界边（T_right != -1），函数必须遍历测试和试验空间的局部基函数，以计算对矩阵A的贡献。
%     这涉及到计算边上的积分，通常需要在高斯点上求出基函数的值、它们的梯度以及系数函数的值
%     if(T_right~=-1)    %检查边是否在边界上
         vertices_right=P(:,T((1:4),T_right));   %边右侧单元的顶点坐标
        for beta=1:number_of_local_basis_test 
           for alpha=1:number_of_local_basis_trial 
                
               
                int_value=Gauss_quadrature_1D_trial_test1(normal,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0,data);
                A(Tb_test(beta,T_left),Tb_trial(alpha,T_left))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_left)) -int_value; % 同一三角形项的贡献
                %对其他三角形的贡献 \phi_j^-\phi_i^+
                %计算从相邻三角形（vertices_right）到当前三角形（vertices_left）的积分贡献。
                %这反映了边界上的函数phi_j^-和phi_i^+之间的交互
                int_value=Gauss_quadrature_1D_trial_test1(normal,Gauss_nodes,Gauss_weights,vertices_right,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0,data);
                A(Tb_test(beta,T_left),Tb_trial(alpha,T_right))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_right)) + int_value; %对邻项的贡献
                   
             
           end
        end
        
        %当 T_right 等于 -1 时，这表示边界条件或右侧没有相邻三角形的边
        %代码以不同的方式计算这样的边的贡献，专注于只涉及左侧单元的项
        
   elseif normal(1)*feval(data.b1,(P(1,vertex1)+P(1,vertex2))/2,(P(2,vertex1)+P(2,vertex2))/2)+feval(data.b2,(P(1,vertex1)+P(1,vertex2))/2,(P(2,vertex1)+P(2,vertex2))/2)*normal(2)<0 && T_right ==-1
        for beta=1:number_of_local_basis_test 
           for alpha=1:number_of_local_basis_trial             
                int_value=Gauss_quadrature_1D_trial_test1(normal,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0,data);
                A(Tb_test(beta,T_left),Tb_trial(alpha,T_left))=A(Tb_test(beta,T_left),Tb_trial(alpha,T_left)) -int_value; % Contribution from terms of the same triangle.
          
           end
        end
    end
%  else
%      disp('跳过此边');
%   end
end
end