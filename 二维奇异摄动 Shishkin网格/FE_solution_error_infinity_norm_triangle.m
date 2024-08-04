function r=FE_solution_error_infinity_norm_triangle(u1,accurate_function,geo,para,data,varepsilon)

[P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);   
number_of_elements=size(T,2); 
number_of_local_basis_trial=size(Tb,1); 

if para.basis_type_trace_u==104
    number_of_local_basis_testt=5;
elseif para.basis_type_trace_u==103
   number_of_local_basis_testt=4;
elseif para.basis_type_trace_u==102
    number_of_local_basis_testt=3;
elseif para.basis_type_trace_u==101
    number_of_local_basis_testt=2;
end
 [P3,T3,Pb3,Tb3,Eb3]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);   

% [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(para.order_Gauss);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
for j=1:4   
    current_edge=T(number_of_local_basis_trial+j,n);     %T中记录这该单元三条边的编号索引
        T_left=n;                                                     %它应该总是等于n   %这里定义该单元上所有边左侧的单元索引都是该单元本身
     if (Eb3(number_of_local_basis_testt+1,current_edge)==n)   
        T_right=Eb3(number_of_local_basis_testt+2,current_edge);
    else
        T_right=Eb3(number_of_local_basis_testt+1,current_edge);

     end
     if T_right~=-1
    vertex1=Eb3(1,current_edge);   %边的第一个顶点索引
    vertex2=Eb3(number_of_local_basis_testt,current_edge);    %边的第二个顶点索
    vertices1=[P3(:,vertex1),P3(:,vertex2)];
    
     [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(para.order_Gauss);
     vertices_ww1=[P3(:,vertex1)';P3(:,vertex2)'];
     [Gauss_nodes,Gauss_weights]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_ww1);
uapp_plus3=0;
          uapp_plus=u1(Eb3((number_of_local_basis_testt+3:end),current_edge));    %获取左侧三角形上的数值迹
   for k=1:size(Gauss_nodes,1)
        for beta=1:number_of_local_basis_testt
                if Gauss_nodes(1,1)-Gauss_nodes(1,2)~=0
                 uapp_plus3=uapp_plus3 +uapp_plus(beta)*FE_local_basis_1D_local(Gauss_nodes(1,k),vertices1(1,:),para.basis_type_trace_u,beta,0);

               elseif Gauss_nodes(1,1)-Gauss_nodes(1,2)==0
                  uapp_plus3=uapp_plus3 +uapp_plus(beta)*FE_local_basis_1D_local(Gauss_nodes(2,k),vertices1(2,:),para.basis_type_trace_u,beta,0);
                end
        end 
        temp=max(abs(accurate_function(Gauss_nodes(1,k),Gauss_nodes(2,k),varepsilon)-uapp_plus));

                 if temp>r
                 r=temp;
                 end
   end
     end
end

end






