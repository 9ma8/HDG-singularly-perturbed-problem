function int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices1,vertices,basis_type_trace_u_trial,basis_index_trial,der_trial,basis_type_u_test,basis_index_test,der_x_test,der_y_test)
%Gauss_quadrature_2D_volume_trial_test computes the integral using Gaussian
%quadrature. It receives the gauss nodes and weigths, the vertices of the element, the basis_type, 
%the basis index and the derivative order for the FEM space both for trial
%and test FE spaces.

 Gpn=size(Gauss_nodes,2);    %Gauss��ĸ���


%  int_value=0;    %��ʼ����һ��ֵ
r=0;
% c=zeros(Gpn,1);
for k=1:Gpn
%     c(k)=feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k));

if Gauss_nodes(1,1)-Gauss_nodes(1,2)~=0
r=r+Gauss_weights(1,k)*(feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_1D_local(Gauss_nodes(1,k),vertices1(1,:),basis_type_trace_u_trial,basis_index_trial,der_trial).*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_u_test,basis_index_test,der_x_test,der_y_test));
elseif Gauss_nodes(1,1)-Gauss_nodes(1,2)==0
%     Gauss_weights(1,k);
r=r+Gauss_weights(1,k)*(feval(coe_fun,Gauss_nodes(1,k),Gauss_nodes(2,k))*FE_local_basis_1D_local(Gauss_nodes(2,k),vertices1(2,:),basis_type_trace_u_trial,basis_index_trial,der_trial).*FE_local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_u_test,basis_index_test,der_x_test,der_y_test));
end
end
int_value=r;
end







