%------------------------------------------------------------------------

% example_DG��һ���ڱ߽��ϴ���Dirichlet�߽������ļ���ͼ�����ü��٤�ɽ𷽷���Ⲵ������Ľű�ʾ����
%�û�����ͨ��ѡ����ʵġ�Poisson_solver����ʹ���ڲ��ͷ���(IP)���ӽ��ڲ��ͷ���(IPH)��

% ��һ��������������Ľ�����

%ע��:Ŀǰ������generate_boundarynodesֻ�ڼ�������������(0,1)^2ʱ������



%-------------------------------------------------------------------------
clear all;

geo.right=1;
geo.left=0;
geo.top=1;
geo.bottom=0;

% �ṹ���ݰ����й��������ݵ���Ϣ������Դ���ɢϵ�����߽�����
data.c=@(x,y,El) 1;     %��ɢ��
data.c1=@(x,y,El) 1;     %��Ӧ��
data.Dirichlet_fun=@(x,y) 0;       
data.sigma=3;
data.beta1=1;
data.beta2=2;

data.a=@(x, y) 2+3*y^2;
data.b1=@(x, y) 2-x;
data.b2=@(x,y) 3-y^3;

para.order_Gauss=3;     %�����Ľ�
% varepsilon=1e-4;
% para.alpha_coef=varepsilon;      %DG�ͷ�����
para.basis_type_u=2;    
para.basis_type_q=2;
para.basis_type_trace_u=102;
% tau=1;

% �������

uex=@(x,y,varepsilon) (y.^3).*sin(x).*(1-exp(-(1-x)/varepsilon)).*(1-exp(-2*(1-y)/varepsilon));     %��ȷ��
ux=@(x,y,varepsilon) (cos(x).*(1-exp(-(1-x)/varepsilon))-sin(x)*exp(-(1-x)/varepsilon)/varepsilon).*(y^3).*(1-exp(-2*(1-y)/varepsilon));
uy=@(x,y,varepsilon) (3*(y^2).*(1-exp(-2*(1-y)/varepsilon))-2*(y^3)*exp(-2*(1-y)/varepsilon)/varepsilon).*sin(x).*(1-exp(-(1-x)/varepsilon));
uxx=@(x,y,varepsilon) (-sin(x).*(1-exp(-(1-x)/varepsilon))-2*cos(x)*exp(-(1-x)/varepsilon)/varepsilon-sin(x)*exp(-(1-x)/varepsilon)/(varepsilon^2)).*(y^3).*(1-exp(-2*(1-y)/varepsilon));
uyy=@(x,y,varepsilon) (6*y.*(1-exp(-2*(1-y)/varepsilon))-12*(y^2)*exp(-2*(1-y)/varepsilon)/varepsilon-4*(y^3)*exp(-2*(1-y)/varepsilon)/(varepsilon^2)).*sin(x).*(1-exp(-(1-x)/varepsilon));
data.f=@(x,y,varepsilon) -varepsilon*(uxx(x,y,varepsilon)+uyy(x,y,varepsilon))+(2-x).*ux(x,y,varepsilon)+(3-y^3).*uy(x,y,varepsilon)+uex(x,y,varepsilon);     %�Ҷ���
% hh=1/2;
pex=@(x,y,varepsilon) -varepsilon*ux(x,y,varepsilon);
qex=@(x,y,varepsilon) -varepsilon*uy(x,y,varepsilon);

epsilon=[1e-5;1e-6;1e-7;1e-8];
N=[4;8;16;32;64;128;256];

outputFile = fopen('superclosenes2.txt', 'w');

%=================================
%   ���������
%=================================

% ѭ������ÿ����ֵ
for j = 1:length(epsilon)
    varepsilon = epsilon(j);
   fprintf(outputFile, 'Results for varepsilon= %e\n', varepsilon);
    fprintf(outputFile, '%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n', 'N', 'supercloseness_L2 Error', 'supercloseness_Lq Error', 'supercloseness_DG Error', 'L2 Conv Order', 'Lq Conv Order','DG Conv Order');
    
    error_L2 = zeros(length(N), 1);
    error_Lq = zeros(length(N), 1);
    error_DG = zeros(length(N), 1);
    convergence_orders_L2 = zeros(length(N)-1, 1);
    convergence_orders_Lq = zeros(length(N)-1, 1);
    convergence_orders_DG = zeros(length(N)-1, 1);

for i=1:length(N)
geo.N=N(i);
[P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,para.basis_type_u,para.basis_type_trace_u,data,varepsilon);
tau=3;

[uhat, local] = hdg_Poisson_GlobalSystem(P,T,Tb,Eb,data.c,data.f,para.basis_type_trace_u,para.basis_type_u,tau,para.order_Gauss,geo,data,varepsilon);
[u,q1,q2] = hdg_Poisson_LocalProblem(Pb,T,Tb,Eb,uhat,local,para.basis_type_trace_u);

% hvec(i)=hmax;    %�ߵ���󳤶�
for k=1:size(Pb,2)     %�����ɶȼ��㾫ȷ��
	x=Pb(1,k); y=Pb(2,k);
	uexvec(k)=uex(x,y,varepsilon);
    pexvec(k)=pex(x,y,varepsilon);
    qexvec(k)=qex(x,y,varepsilon);
end

%  [error_L2(i),error_Lq(i),error_DG(i)]=compute_errors_HDG(u,q1,q2,uhat,uex,pex,qex,para,'DG',geo,tau,data,varepsilon)   %����L2,H1��DG���
 [error_supercloseness_L2(i),error_supercloseness_Lq(i),error_supercloseness_DG(i)]=compute_supercloseness_errors_HDG(u,q1,q2,uhat,uexvec,pexvec,qexvec,uex,para,'DG',geo,tau,data,varepsilon)   %����L2,H1��DG���ƽ����
% [maxerror_L2(i),maxerror_L2_trace(i)]=compute_errors_max(u,uhat,uex,uexvec,para,geo,data,varepsilon)
% infinity_norm_trace_u(i)=FE_solution_error_infinity_norm_triangle(uhat,uex,geo,para,data,varepsilon)
% end



if i>1
%     convergence_orders_L2(i-1) = log(error_L2(i-1)/error_L2(i))/log(2)
%     convergence_orders_Lq(i-1) = log(error_Lq(i-1)/error_Lq(i))/log(2*log(geo.N)/log(2*geo.N))
%     convergence_orders_DG(i-1) = log(error_DG(i-1)/error_DG(i))/log(2*log(geo.N)/log(2*geo.N))
    convergence_orders_supercloseness_L2(i-1) = log(error_supercloseness_L2(i-1)/error_supercloseness_L2(i))/log(2)
    convergence_orders_supercloseness_Lq(i-1) = log(error_supercloseness_Lq(i-1)/error_supercloseness_Lq(i))/log(2*log(geo.N)/log(2*geo.N))
    convergence_orders_supercloseness_DG(i-1) = log(error_supercloseness_DG(i-1)/error_supercloseness_DG(i))/log(2*log(geo.N)/log(2*geo.N))
%     convergence_orders_maxerror_L2(i-1) = log(maxerror_L2(i-1)/maxerror_L2(i))/log(2*log(geo.N)/log(2*geo.N))
%     convergence_orders_maxerror_L2_trace(i-1) = log(maxerror_L2_trace(i-1)/maxerror_L2_trace(i))/log(2*log(geo.N)/log(2*geo.N))

end

if i == 1
      fprintf(outputFile, '%d\t%.5e\t%.5e\t%.5e\t%s\t%s\n', N(i), error_supercloseness_L2(i), error_supercloseness_Lq(i), error_supercloseness_DG(i), 'N/A', 'N/A');
 else
      fprintf(outputFile, '%d\t\t%.5e\t\t%.5e\t\t%.5e\t\t%.5f\t\t%.5f\t\t%.5f\n', N(i), error_supercloseness_L2(i), error_supercloseness_Lq(i), error_supercloseness_DG(i), convergence_orders_supercloseness_L2(i-1), convergence_orders_supercloseness_Lq(i-1), convergence_orders_supercloseness_DG(i-1));
end
end

  fprintf(outputFile, '\n'); % ��ӿ����Էָ���ͬ�ŵĽ��
end
% �ر��ļ�
fclose(outputFile);





