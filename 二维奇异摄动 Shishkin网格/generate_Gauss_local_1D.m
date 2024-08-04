function [Gauss_point_local_1D,Gauss_coefficient_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices)

%Generate the Gauss coefficients and Gauss points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation.
%Gauss_coefficient_local_1D,Gauss_point_local_1D:the Gauss coefficients and Gauss points on an arbitrary interval.
%Gauss_coefficient_reference,Gauss_point_reference: the Gauss coefficients and Gauss points on the reference interval [-1,1].
%将[-1, 1]上的Gauss点和权重仿射变换到任意的小区间上有如下变换
% edge_length=norm(vertices1-vertices2,2);
 edge_length = sqrt((vertices(2,1) - vertices(1,1))^2 + (vertices(2,2) - vertices(1,2))^2);
Gauss_coefficient_local_1D = edge_length * Gauss_coefficient_reference_1D / 2;
Gauss_point_local_1D = zeros(2,length(Gauss_point_reference_1D));

    
 for i = 1:length(Gauss_point_reference_1D)   
if vertices(2,2)-vertices(1,2)==0
Gauss_point_local_1D(1,i)=(vertices(2,1)-vertices(1,1))*Gauss_point_reference_1D(i)/2+(vertices(2,1)+vertices(1,1))/2;
Gauss_point_local_1D(2,i)=vertices(2,2);
else
Gauss_point_local_1D(1,i)=vertices(1,1);
Gauss_point_local_1D(2,i)=(vertices(2,2)-vertices(1,2))*Gauss_point_reference_1D(i)/2+(vertices(2,2)+vertices(1,2))/2;
end
% Gauss_point_local_1D=Gauss_point_local_1D';
 end
 
 
 
 
 
 
 
 
 
 
 
 