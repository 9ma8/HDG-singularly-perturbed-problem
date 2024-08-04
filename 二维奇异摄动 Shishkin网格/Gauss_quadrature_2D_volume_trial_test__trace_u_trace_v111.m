function int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_trace_v111(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trace_u_trial,basis_index_trial,der_trial,basis_type_trace_u_test,basis_index_test,der_test)
%Gauss_quadrature_2D_volume_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space both for trial
%and test FE spaces.

 Gpn=size(Gauss_nodes,2);    %Gauss��ĸ���


 r=0;    %��ʼ����һ��ֵ

% c=zeros(Gpn,1);
for k=1:Gpn
%     c(k)=feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k));
% end
if Gauss_nodes(1,1)-Gauss_nodes(1,2)~=0
r=r+Gauss_weights(1,k)*(feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_111D(Gauss_nodes(1,k),vertices(1,:),basis_type_trace_u_trial,basis_index_trial,der_trial).*FE_local_basis_111D(Gauss_nodes(1,k),vertices(1,:),basis_type_trace_u_test,basis_index_test,der_test));
elseif Gauss_nodes(1,1)-Gauss_nodes(1,2)==0
r=r+Gauss_weights(1,k)*(feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_111D(Gauss_nodes(2,k),vertices(2,:),basis_type_trace_u_trial,basis_index_trial,der_trial).*FE_local_basis_111D(Gauss_nodes(2,k),vertices(2,:),basis_type_trace_u_test,basis_index_test,der_test));
end

end
int_value=r;
end
% int_value=int_value(1);








