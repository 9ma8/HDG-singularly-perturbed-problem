function [maxerror_L2,maxerror_L2_trace]=compute_errors_max(u,u1,uex,uexpex,para,geo,data,varepsilon)

maxerror_L2=0;
maxerror_L2_trace=0;
[P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);
matrixsize1=size(Pb,2); 
% error_H1=0;  %初始化L2范数误差、H1范数误差和DG误差为零
number_of_elements=size(T,2); %计算三角形单元的数量

number_of_local_basis_trial=size(Tb,1);    %计算局部基函数的数量
if para.basis_type_trace_u==104
    number_of_local_basis_testt=5;
elseif para.basis_type_trace_u==103
   number_of_local_basis_testt=4;
elseif para.basis_type_trace_u==102
    number_of_local_basis_testt=3;
elseif para.basis_type_trace_u==101
    number_of_local_basis_testt=2;
end
for n=1:number_of_elements % Loop over the triangles.
   for beta=1:number_of_local_basis_trial
    temp1=uexpex(Tb(beta,n))'-u(Tb(beta,n)); 
    if abs(maxerror_L2)<abs(temp1)
         maxerror_L2=abs(temp1);
    end
    end
    for j=1:4   
    current_edge=T(number_of_local_basis_trial+j,n);     %T中记录这该单元三条边的编号索引
    if (Eb(number_of_local_basis_testt+1,current_edge)==n)   
        T_right=Eb(number_of_local_basis_testt+2,current_edge);
%     else
%         T_right=Eb(3,current_edge);
%     end
if T_right~=-1   

    for alpha=1:number_of_local_basis_testt
       x=P(1, Eb(alpha,current_edge));
       y=P(2, Eb(alpha,current_edge));
    temp2=uex(x,y,varepsilon)-u1(Eb(number_of_local_basis_testt+2+alpha,current_edge)); 
     if abs(maxerror_L2_trace)<abs(temp2)
         maxerror_L2_trace=abs(temp2);
    end
    end
end
   
    end
    end
end
end

    








