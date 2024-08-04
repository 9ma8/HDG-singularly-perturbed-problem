%-------------------------------------------------------------------------
%generate_mesh_2D在几何函数geo.fun指定的几何体上生成具有最大网格尺寸geo.h的网格。
% generate_mesh_2D接收
% geo:包含描述几何图形和网格大小geo.h的函数句柄以及规则/非规则网格的标志的结构。
% basis_type: FEM空间的类型(线性、二次元等)。
% generate_mesh_2D返回:
% P:包含网格顶点位置的顶点矩阵。例如:P (: 1) = (X_1; Y_1)
% T:包含单元的索引顶点加上单元k所属的子域的单元矩阵。
% Eg: T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain]。如果使用DG，在T中我们也保存了三角形边的索引。
%E：边矩阵。第一行和第二行包含边的起点和终点的索引，第三行和第四行包含起点和终点参数值，
%第五行包含边段编号，第六行和第七行包含左侧和右侧子域的索引。

%Pb：自由度的顶点矩阵（如果使用P1, Pb=P）。

%Tb：在每个单元上具有独立自由度的单元矩阵（如果使用P1, Tb=T）

%Eb：边矩阵。在Eb中，第一行和第二行包含起始和结束顶点第三行和第四行根据逆时针方向在左侧和右侧包含三角形

%h：标量变量，等于边的最小长度？？？?

%-------------------------------------------------------------------------


function [P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,basis_type,basis_type_trace_u,data,varepsilon)


   N=geo.N;
   

%====DG P1 FE的数据结构
%输出：Eb第一行和第二行包含起始顶点和结束顶点的索引，
%第三行和第四行根据逆时针方向在左侧和右侧包含三角形

tau1=min(1/2,data.sigma*varepsilon*log(N)/data.beta1);
tau2=min(1/2,data.sigma*varepsilon*log(N)/data.beta2);
%输出：T有7行。3个顶点3条边
if basis_type==1
%     N1=(geo.right-geo.left)/h;
%    N2=(geo.top-geo.bottom)/h;
     h_11=2*(geo.right-tau1)/N;
     h_12=2*(tau1)/N;
     h_21=2*(geo.top-tau2)/N;
     h_22=2*(tau2)/N;
     
   tnp=(N+1)*(N+1);
   P=zeros(2,tnp);
   T=zeros(3,N*N);
   Q=zeros(N+1,N+1);

  for j=1:N+1
       for i=1:N+1
            if j<=N/2+1
                y(j)=geo.bottom+(j-1)*h_21;
            else
                y(j)=geo.top-tau2+(j-(N)/2-1)*h_22;
            end  

            if i<=(N)/2+1
               x(i)=geo.left+(i-1)*h_11;
            else
              x(i)=geo.right-tau1+(i-(N)/2-1)*h_12;
            end
       end
  end
 for j=1:tnp
         for i=1:tnp
           P(1,i)=x(1+fix((i-1)/(N+1)));
           if mod(j,N+1)==0
               P(2,j)=geo.top;
           else
            P(2,j)=y(mod(j,(N+1)));
           end
         end
 end

   for i=1:N+1
      for j=1:N+1
         Q(i,j)=(i-1)*(N+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
  for n=1:N*N
   
      if mod(n,N)==0
         row=N;
         column=n/N;
      else
         row=mod(n,N);
         column=fix(n/N)+1;
      end
   
      T(1,n)=Q(column,row);
      T(2,n)=Q(column+1,row);
      T(3,n)=Q(column+1,row+1); 
      T(4,n)=Q(column,row+1);
  

   end

%     Tb=T(1:end,:);
   Tb = zeros(size(T));
    T=[T; zeros(4,size(T,2))];  
    j=0;
    for k=1:size(T,2)    % 遍历所有的三角形单元
          for i=1:4      %将每个三角形的三个顶点的坐标添加到 Pb 中
              j=j+1;
             Pb(:, j)=P(:,T(i,k));
          end
          Tb((1:4),k)=[j-3;j-2;j-1;j];     %更新 Tb，指定每个三角形的三个顶点对应的自由度
    end
    
    
    %% ============= BUILD MATRIX E_B
    
 %Eb 是一个矩阵，用于存储三角形网格的边界信息。
% 通过遍历每个三角形，将边界信息存储在 Eb 中。
% 对于每个三角形，首先将其三个顶点的索引存储在 Eb 的前三行中，并将其三个边标记为 1。
% 然后，检查每个三角形的三条边是否已经在 Eb 中。如果不在，就将边的信息添加到 Eb 中，并更新三角形中边的索引。
    
    
max_num_edges= 2*size(P,2)- 6;     %预分配足够大的矩阵Eb的估计
num_edges=1;
Eb=-ones(4,max_num_edges);    %将所有直角三角形的初始值设定为边的“边界”。如果不是这样，我们更改值，请插入打印T_right
%根据三角形数量预分配足够大的矩阵空间，并初始化 Eb

vertex1=T(1,1);     %第一个三角形，插入所有边。
vertex2=T(2,1);
vertex3=T(3,1);
vertex4=T(4,1);
Eb(1,num_edges)=vertex1;
Eb(2,num_edges)=vertex2;
Eb(3,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex2;
Eb(2,num_edges)=vertex3;
Eb(3,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex3;
Eb(2,num_edges)=vertex4;
Eb(3,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex4;
Eb(2,num_edges)=vertex1;
Eb(3,num_edges)=1;
T(5,1)=1;
T(6,1)=2;
T(7,1)=3;   %在矩阵T中插入边的标签。
T(8,1)=4;
for k=2:size(T,2)    
    vertex1=T(1,k);     %获取当前三角形的三个顶点的索引。
    vertex2=T(2,k);
    vertex3=T(3,k);
    vertex4=T(4,k);
    m=num_edges;
    flag_1=0;
    flag_2=0;
    flag_3=0;
    flag_4=0;
for j=1:m    %检查每条边是否已经存在于 Eb 中
    if( (vertex1==Eb(2,j) && vertex2==Eb(1,j)) || (vertex1==Eb(1,j) && vertex2==Eb(2,j)) )    %检查当前边是否与 Eb 中的某条边相匹配
        Eb(4,j)=k;    %将相邻三角形设为Eb等于k
        T(5,k)=j;     % 把边放在三角形矩阵中。  %将当前边标记为相邻三角形，并在 T 中记录边的索引
        flag_1=1;
    elseif ( (vertex2==Eb(2,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(2,j)) )    %检查第二条边是否已经在里面了。
        Eb(4,j)=k;
        T(6,k)=j;
        flag_2=1;
    elseif ( (vertex3==Eb(2,j) && vertex4==Eb(1,j)) || (vertex3==Eb(1,j) && vertex4==Eb(2,j)))    %检查第三条边是否已经在里面了。
        Eb(4,j)=k;
        T(7,k)=j;
        flag_3=1;
     elseif ( (vertex4==Eb(2,j) && vertex1==Eb(1,j)) || (vertex4==Eb(1,j) && vertex1==Eb(2,j)))    %检查第三条边是否已经在里面了。
        Eb(4,j)=k;
        T(8,k)=j;
        flag_4=1;
    end
end
if(flag_1==0)     %如果当前边不存在于 Eb 中，则添加它
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex1;
    Eb(2,num_edges)=vertex2;
    Eb(3,num_edges)=k;
    T(5,k)=num_edges;
end  
if(flag_2==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex2;
    Eb(2,num_edges)=vertex3;
    Eb(3,num_edges)=k;
    T(6,k)=num_edges;
end
if(flag_3==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex3;
    Eb(2,num_edges)=vertex4;
    Eb(3,num_edges)=k;
    T(7,k)=num_edges;
end

if(flag_4==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex4;
    Eb(2,num_edges)=vertex1;
    Eb(3,num_edges)=k;
    T(8,k)=num_edges;
end
end
Eb=Eb(:,(1:num_edges));    %如果当前边不存在于 Eb 中，则添加它
% 创建两行全为 0 的矩阵
zeros_rows = zeros(2, size(Eb, 2)); % 第一个维度是添加的行数，第二个维度是 Eb 的列数

% 将全为 0 的行添加到矩阵 Eb 的末尾
Eb = [Eb; zeros_rows];
 i=1;
for j=1:size(Eb, 2)
%    if Eb(4, j)~=-1
       Eb(5,j)=i;
       Eb(6,j)=i+1;
       i=i+2;

end



elseif basis_type==2
   dN1=N*2;
     dN2=N*2;
     dh_11=2*(geo.right-tau1)/dN1;
     dh_12=2*(tau1)/dN1;
     dh_21=2*(geo.top-tau2)/dN2;
     dh_22=2*(tau2)/dN2;
   tnp=(dN1+1)*(dN2+1);
   P=zeros(2,tnp);
   T=zeros(3,N*N);
   Q=zeros(dN1+1,dN2+1);

   for j=1:dN2+1
       for i=1:dN1+1
            if j<=(dN2)/2+1
                y(j)=geo.bottom+(j-1)*dh_21;
            else
                y(j)=geo.top-tau2+(j-(dN2)/2-1)*dh_22;
            end  

            if i<=(dN1)/2+1
               x(i)=geo.left+(i-1)*dh_11;
            else
              x(i)=geo.right-tau1+(i-(dN1)/2-1)*dh_12;
            end
       end
  end
 for j=1:tnp
         for i=1:tnp
           P(1,i)=x(1+fix((i-1)/(dN2+1)));
           if mod(j,dN2+1)==0
               P(2,j)=geo.top;
           else
               P(2,j)=y(mod(j,(dN2+1)));
           end
         end
 end


   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
  for n=1:N*N
   
      if mod(n,N)==0
         row=N;
         column=n/N;
      else
         row=mod(n,N);
         column=fix(n/N)+1 ;
      end
   
      T(1,n)=Q(2*column-1,2*row-1);
      T(5,n)=Q(2*column,2*row-1); 
      T(2,n)=Q(2*column+1,2*row-1);
      T(6,n)=Q(2*column+1,2*row);
      T(3,n)=Q(2*column+1,2*row+1);
      T(7,n)=Q(2*column,2*row+1);
      T(4,n)=Q(2*column-1,2*row+1);
      T(8,n)=Q(2*column-1,2*row);
      T(9,n)=Q(2*column,2*row);

  end

    
    Tb=T(1:end,:);
    T=[T; zeros(4,size(T,2))];  %在原始网格 T 的末尾添加三行零向量，用于存储每个三角形的三条边的索引。
    j=0;
    for k=1:size(T,2)    % 遍历所有的三角形单元
          for i=1:9      %将每个三角形的三个顶点的坐标添加到 Pb 中
              j=j+1;
              Pb(:,j)=P(:,T(i,k));
          end
          Tb((1:9),k)=[j-8;j-7;j-6;j-5;j-4;j-3;j-2;j-1;j];     %更新 Tb，指定每个三角形的三个顶点对应的自由度
    end
        
    
    %% ============= BUILD MATRIX E_B
    
max_num_edges= 2*size(P,2)- 6;     %预分配足够大的矩阵Eb的估计
num_edges=1;
Eb=-ones(5,max_num_edges);    %将所有直角三角形的初始值设定为边的“边界”。如果不是这样，我们更改值，请插入打印T_right
%根据三角形数量预分配足够大的矩阵空间，并初始化 Eb
if basis_type_trace_u==104
    number_of_local_basis_trial_trace_u=5;
% number_of_local_basis_test_trace_u=5;
elseif basis_type_trace_u==103
   number_of_local_basis_trial_trace_u=4;
% number_of_local_basis_test_trace_u=4;
elseif basis_type_trace_u==102
    number_of_local_basis_trial_trace_u=3;
%  number_of_local_basis_test_trace_u=3;
elseif basis_type_trace_u==101
    number_of_local_basis_trial_trace_u=2;
%  number_of_local_basis_test_trace_u=2;
end
vertex1=T(1,1);     %第一个三角形，插入所有边。
vertex2=T(2,1);
vertex3=T(3,1);
vertex4=T(4,1);
vertex5=T(5,1);
vertex6=T(6,1);
vertex7=T(7,1);
vertex8=T(8,1);
Eb(1,num_edges)=vertex1;
Eb(2,num_edges)=vertex5;
Eb(3,num_edges)=vertex2;
Eb(4,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex2;
Eb(2,num_edges)=vertex6;
Eb(3,num_edges)=vertex3;
Eb(4,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex3;
Eb(2,num_edges)=vertex7;
Eb(3,num_edges)=vertex4;
Eb(4,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex4;
Eb(2,num_edges)=vertex8;
Eb(3,num_edges)=vertex1;
Eb(4,num_edges)=1;
T(10,1)=1;
T(11,1)=2;
T(12,1)=3;   %在矩阵T中插入边的标签。
T(13,1)=4;
for k=2:size(T,2)    
    vertex1=T(1,k);     %获取当前三角形的三个顶点的索引。
    vertex2=T(2,k);
    vertex3=T(3,k);
    vertex4=T(4,k);
    vertex5=T(5,k);
  vertex6=T(6,k);
 vertex7=T(7,k);
  vertex8=T(8,k);
    m=num_edges;
    flag_1=0;
    flag_2=0;
    flag_3=0;
    flag_4=0;
for j=1:m    %检查每条边是否已经存在于 Eb 中
    if( (vertex1==Eb(number_of_local_basis_trial_trace_u,j) && vertex2==Eb(1,j)) || (vertex1==Eb(1,j) && vertex2==Eb(number_of_local_basis_trial_trace_u,j)) )    %检查当前边是否与 Eb 中的某条边相匹配
        Eb(5,j)=k;    %将相邻三角形设为Eb等于k
        T(10,k)=j;     % 把边放在三角形矩阵中。  %将当前边标记为相邻三角形，并在 T 中记录边的索引
        flag_1=1;
    elseif ( (vertex2==Eb(number_of_local_basis_trial_trace_u,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(number_of_local_basis_trial_trace_u,j)) )    %检查第二条边是否已经在里面了。
        Eb(5,j)=k;
        T(11,k)=j;
        flag_2=1;
    elseif ( (vertex3==Eb(number_of_local_basis_trial_trace_u,j) && vertex4==Eb(1,j)) || (vertex3==Eb(1,j) && vertex4==Eb(number_of_local_basis_trial_trace_u,j)))    %检查第三条边是否已经在里面了。
        Eb(5,j)=k;
        T(12,k)=j;
        flag_3=1;
     elseif ( (vertex4==Eb(number_of_local_basis_trial_trace_u,j) && vertex1==Eb(1,j)) || (vertex4==Eb(1,j) && vertex1==Eb(number_of_local_basis_trial_trace_u,j)))    %检查第三条边是否已经在里面了。
        Eb(5,j)=k;
        T(13,k)=j;
        flag_4=1;
    end
end
if(flag_1==0)     %如果当前边不存在于 Eb 中，则添加它
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex1;
      Eb(2,num_edges)=vertex5;
    Eb(3,num_edges)=vertex2;
    Eb(4,num_edges)=k;
    T(10,k)=num_edges;
end  
if(flag_2==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex2;
    Eb(2,num_edges)=vertex6;
    Eb(3,num_edges)=vertex3;
    Eb(4,num_edges)=k;
    T(11,k)=num_edges;
end
if(flag_3==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex3;
     Eb(2,num_edges)=vertex7;
    Eb(3,num_edges)=vertex4;
    Eb(4,num_edges)=k;
    T(12,k)=num_edges;
end

if(flag_4==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex4;
     Eb(2,num_edges)=vertex8;
    Eb(3,num_edges)=vertex1;
    Eb(4,num_edges)=k;
    T(13,k)=num_edges;
end
end
Eb=Eb(:,(1:num_edges));    %如果当前边不存在于 Eb 中，则添加它
Eb=[Eb;zeros(3,size(Eb,2))];    
  
 i=1;
for j=1:size(Eb, 2)
       Eb(6,j)=i;
       Eb(7,j)=i+1;
       Eb(8,j)=i+2;
       i=i+3;

end

% for k=1:num_edges
%         T_1= Eb(3,k);
%         T_2=Eb(4,k);
%         if (T1==1)
            




end

end
