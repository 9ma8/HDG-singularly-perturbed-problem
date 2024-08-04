function normal=find_normal(vertex1,vertex2,vertex_left_index,P)
% 函数find_normal接收两个顶点的两个索引和一个三个顶点索引的向量。
%它找到第三个缺失的索引顶点3。
%然后计算到边的切向量，法向量计算标量积btw和法向量决定法向量的符号。
%-----------------------------------------------------------

%其目的是根据给定的顶点和三角形的顶点索引来计算三角形的法向量

normal=zeros(2,1);     %初始化一个2x1的零向量 normal，用于存储计算得到的法向量

% 根据给定的顶点索引 vertex1 和 vertex2 以及三角形的顶点索引 vertex_left_index，
% 确定三角形的第三个顶点索引 vertex3
if(vertex1==vertex_left_index(1) && vertex2==vertex_left_index(2) || vertex2==vertex_left_index(1) && vertex1==vertex_left_index(2))
         vertex3=vertex_left_index(3);
    elseif(vertex1==vertex_left_index(2) && vertex2==vertex_left_index(3) || vertex2==vertex_left_index(2) && vertex1==vertex_left_index(3)) 
         vertex3=vertex_left_index(4);
    elseif(vertex1==vertex_left_index(3) && vertex2==vertex_left_index(4) || vertex2==vertex_left_index(3) && vertex1==vertex_left_index(4)) 
         vertex3=vertex_left_index(1);
    elseif(vertex1==vertex_left_index(4) && vertex2==vertex_left_index(1) || vertex2==vertex_left_index(4) && vertex1==vertex_left_index(1)) 
         vertex3=vertex_left_index(2);
end

% 计算顶点 vertex1 到顶点 vertex2 的向量 vintern
%计算两顶点之间的切向量 tau，并将其归一化，得到单位切向量
    tau=P(:,vertex2)-P(:,vertex1);    %计算切矢量
    tau=tau/norm(tau);      %归一化切向量（单位化）
    %利用单位切向量 tau 计算法向量，其x分量为切向量y分量的负值，
    %y分量为切向量x分量的值，从而得到与切向量垂直的法向量
    
    normal(1)=-tau(2);
    normal(2)=tau(1);
    
    %计算顶点 vertex1 到顶点 vertex3 的向量 vintern，
    %并将其归一化，得到单位向量
    vintern=P(:,vertex3)-P(:,vertex1);
    vintern=vintern/norm(vintern);
    
    
    %检查向量 vintern 和法向量 normal 的点积是否为正值，
    %如果是，则将法向量取反，以确保法向量指向正确的方向
    if vintern'*normal >0
        normal=-normal;
    end
end

%这段代码的作用是计算给定顶点和三角形的顶点索引构成的三角形的法向量，
%并确保外法向量的方向正确


