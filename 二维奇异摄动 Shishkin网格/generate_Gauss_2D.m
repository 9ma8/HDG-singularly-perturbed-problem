function [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,type)
%函数'generate_Gauss_2D'在向量顶点定义的当前单元上传递高斯节点。


%这段代码的功能是计算在参考单元上的高斯积分点和权重，
%并进行必要的坐标变换，以便在实际单元上进行积分。



%=== 线性映射 
A=[(vertices(1,3)-vertices(1,1))/2,0; 0,(vertices(2,3)-vertices(2,1))/2];
%代码开始定义了一个变换矩阵 A，用于将参考单元上的高斯积分点映射到实际单元上。
%这个变换矩阵 A 是根据单元的顶点坐标计算得到的。
if type >=1  %根据积分类型 type 计算高斯积分点和权重：
    %如果 type 大于等于 1，则计算在参考单元上的高斯积分点和权重。这里使用了一维区间 ([-1, 1]) 上
    %的高斯-勒让德积分节点和权重计算函数 gauleg。
    w_2D=[];
    node_2D=[];
    % 一维区间[-1,1]上的节点
    n=type;
    [x,w]= gauleg(-1,1,n);
    for i=1:n
        for j=1:n
            node_2D=[node_2D; x(i), x(j)];
            w_2D = [w_2D,w(i).*w(j)];
        end
    end
    v=[(vertices(1,1)+vertices(1,3))/2; (vertices(2,1)+vertices(2,3))/2];
    v=repmat(v,1,size(node_2D,1));
    Gauss_nodes=node_2D;
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';  %行向量
    Gauss_weights=w_2D;
elseif type==-1  %如果 type 等于 -1，则直接使用单元的顶点作为高斯积分点，并分配相同的权重。
    Gauss_nodes=vertices';
    Gauss_weights=[1/3,1/3,1/3];
elseif type==-2    %如果 type 等于 -2，则计算特定于 P2b（二次基函数）的高斯积分点和权重，这是一个针对特定问题定制的积分方法。
    %Lumped for P2b 
    Gauss_nodes=[0,0;1,0;0,1;0.5,0;0,0.5;0.5,0.5;1/3,1/3];
    v=[vertices(1,1); vertices(2,1)];
    v=[v,v,v,v,v,v,v];
    
    %在计算高斯积分点后，代码根据变换矩阵 A 对这些积分点进行坐标变换，将它们映射到实际单元上。
    %这是因为在有限元方法中，通常需要在参考单元上进行积分，而实际单元上的积分需要通过坐标变换进行转换。
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';
    Gauss_weights=[1/20,1/20,1/20,2/15,2/15,2/15,9/20];
end
    
 
%代码定义了一个函数 gauleg，用于计算一维区间上的高斯-勒让德积分节点和权重。
function [x,w]= gauleg(a,b,n)

m=(n+1)/2;
xm=0.5*(b+a);
xl=0.5*(b-a);
xx=[];

for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
        if (abs(z-z1)<eps), break, end
    end
    xx(i)=xm-xl*z;
    xx(n+1-i)=xm+xl*z;
    ww(i)=2.0*xl/((1.0-z*z)*pp*pp);
    ww(n+1-i)=ww(i);
end

x=xx;
w=ww;
end

end
