%-------------------------------------------------------------------------
% assemble_matrix_2D接收
% coe_fun:扩散系数
%网格矩阵:T,P
% Tb_trial,Tb_test:包含trial和test有限单元空间信息的矩阵
% matrixsize1,matrixsize2:刚度矩阵的大小
% basis_type_trial,basis_type_test:有限元空间类型试验/测试
% der_x_trial,der_x_test, der_y_trial,der_y_test: bilinea形式的导数的阶数
% order_gauss:高斯正交阶数。
% assemble_matrix_2D returns: 
% A (matrix) whose entry A_{i,j}= \int (c(\partial_x)^der_x_trial (\partial_y)^der_y_trial \phi_j)((\partial_x)^der_x_test (\partial_y)^der_y_test \phi_j)

%这段代码的核心部分是通过高斯积分来计算在每个三角形上的积分值，并根据计算结果组装刚度矩阵

%-------------------------------------------------------------------------

function A=assemble_matrix_2D_qu(coe_fun,P,T,Tb_q_trial,Tb_test,matrixsize1,matrixsize2,basis_type_q_trial,der_x_trial,der_y_trial,basis_type_u_test,der_x_test,der_y_test,order_Gauss)
%===创建稀疏矩阵并定义参数
number_of_elements=size(T,2); %单元的个数
number_of_local_basis_trial=size(Tb_q_trial,1);    %在单个单元上trial局部基函数的次数
number_of_local_basis_test=size(Tb_test,1);   %返回行数
A=sparse(matrixsize1,matrixsize2);

%=== 循环到组装
for n=1:number_of_elements    %遍历所有的单元
    vertices=P(:,T((1:4),n));    %当前单元的顶点。单元n，顶点T(:,n)，位置P(:,T(:,n))。
    J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);    %计算线性变换的雅可比矩阵
   
    [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss);    %在单元上创建高斯节点。1 .重心  2 .阶数精度。
    for alpha=1:number_of_local_basis_trial    %循环遍历单元上的trial函数以获得试有限元空间。
        for beta=1:number_of_local_basis_test      
%             int_value=0;
            %调用高斯积分方法 Gauss_quadrature_2D_volume_trial_test 来计算在当前三角形上基函数乘积的积分值
            int_value=J*Gauss_quadrature_2D_volume_trial_test_qu(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_q_trial,alpha,der_x_trial,der_y_trial,basis_type_u_test,beta,der_x_test,der_y_test);
            A(Tb_test(beta,n),Tb_q_trial(alpha,n))=A(Tb_test(beta,n),Tb_q_trial(alpha,n))+int_value;

        end
    end
    
end
% A = sparse (Ivec(1:count),Jvec(1:count),Xvec(1:count),matrixsize1,matrixsize2);  %利用 sparse 函数根据 Ivec、Jvec、Xvec 中的数据构建稀疏矩阵 A
end