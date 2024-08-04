function [Zql1,Zql2,Zul,zqf1,zqf2,zuf,Ke,fe] = hdg_Poisson_ElementalMatrices(P,T,Tb,Eb,coe_fun,f_fun,basis_type_test,basis_type_trial,basis_type_trace_u,tau,vertices,order_Gauss,i,data,varepsilon)

%��ȡԪ�ص���������������tau������
%����Ԫ�صĲ������Ժ����������й�
nOfFaces = 4;

number_of_local_basis_trial=size(Tb,1);  
number_of_local_basis_test=size(Tb,1);  

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
%����Ԫ�صĽڵ�������λ�����ɶȺ��������ɶȡ�
nOfElementNodes = size(Tb,1);
ndofU = nOfElementNodes;
ndofQ = ndofU;

%����߽��ϵ������ɶ�����
%����ÿ���棬������Ķ���ʽ�����ۼ�ÿ����Ľڵ���
ndofUHat =4*number_of_local_basis_trial_trace_u;


% Initialisation
fq = zeros(ndofQ,1);
fu = zeros(ndofU,1);
fl = zeros(ndofUHat,1);
Auu = zeros(ndofU,ndofU);
Auu1= zeros(ndofU,ndofU);
Auu2= zeros(ndofU,ndofU);
Auu3= zeros(ndofU,ndofU);
Auq1 = zeros(ndofU,ndofQ);
Auq2 = zeros(ndofU,ndofQ);
Aul = zeros(ndofU,ndofUHat);
Aul1= zeros(ndofU,ndofUHat);
Aul2= zeros(ndofU,ndofUHat);
Aqq1 = zeros(ndofQ,ndofQ);
Aqq2 = zeros(ndofQ,ndofQ);
Aql1 = zeros(ndofQ,ndofUHat);
Aql2 = zeros(ndofQ,ndofUHat);
All = zeros(ndofUHat,ndofUHat);

All1= zeros(ndofUHat,ndofUHat);
All2= zeros(ndofUHat,ndofUHat);
Q=zeros(ndofU,ndofU);

% Element computation -----------------------------------------------------
  J=abs((vertices(1,3)-vertices(1,1))*(vertices(2,3)-vertices(2,1))/4);    %�������Ա任���ſɱȾ���
   
 [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,order_Gauss);    %�ڵ�Ԫ�ϴ�����˹�ڵ㡣1 .����  2 .�������ȡ�

    for alpha=1:number_of_local_basis_trial    %ѭ��������Ԫ�ϵ�trial�����Ի��������Ԫ�ռ䡣
        for beta=1:number_of_local_basis_test   %ע�⣬���������Ӧ��ΪV_{m}�л��������������������� 
%             int_value=0;
            %���ø�˹���ַ��� Gauss_quadrature_2D_volume_trial_test �������ڵ�ǰ�������ϻ������˻��Ļ���ֵ
            int_value1=J*Gauss_quadrature_2D_volume_trial_test_uq(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,1,0,basis_type_test,beta,0,0);
            Auq1(beta,alpha)=Auq1(beta,alpha)+int_value1;
            int_value2=J*Gauss_quadrature_2D_volume_trial_test_uq(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,1,basis_type_test,beta,0,0);
            Auq2(beta,alpha)=Auq2(beta,alpha)+int_value2;
            int_value3=J*Gauss_quadrature_2D_volume_q_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
            Aqq1(beta,alpha)=Aqq1(beta,alpha)-(varepsilon^(-1))*int_value3;
           int_value4=J*Gauss_quadrature_2D_volume_q_trial_test(coe_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
            Aqq2(beta,alpha)=Aqq2(beta,alpha)-(varepsilon^(-1))*int_value4;
            int_value5=J*Gauss_quadrature_2D_volume_trial_test(data.a,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
            Auu1(beta,alpha)=Auu1(beta,alpha)+int_value5;
              int_value6=J*Gauss_quadrature_2D_volume_trial_test(data.b1,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,1,0);
            Auu2(beta,alpha)=Auu2(beta,alpha)-int_value6;
            int_value7=J*Gauss_quadrature_2D_volume_trial_test(data.b2,Gauss_nodes,Gauss_weights,vertices,basis_type_trial,alpha,0,0,basis_type_test,beta,0,1);
            Auu3(beta,alpha)=Auu3(beta,alpha)-int_value7;
        end
    end
for beta=1:number_of_local_basis_test   
 fu(beta,1)=fu(beta,1)+J*Gauss_quadrature_2D_volume_test(f_fun,Gauss_nodes,Gauss_weights,vertices,basis_type_test,beta,0,0,varepsilon);
end


% Face computation --------------------------------------------------------
indexFlip = zeros(1,ndofUHat);
indexFaceIni = 1;
% nOfNodesPreviousFaces = 0;
for iFace = 1:nOfFaces
    current_edge=T(number_of_local_basis_trial+iFace, i);     %T��ÿ����Ԫ�ϵĵ�5��6��7�м�¼��õ�Ԫ�����ߵı������
    T_left=i;       %��Ӧ�����ǵ���n   %���ﶨ��õ�Ԫ�����б����ĵ�Ԫ�������Ǹõ�Ԫ����
       indexFaceEnd = indexFaceIni + number_of_local_basis_trial_trace_u  - 1;
        indexFaceV = indexFaceIni:indexFaceEnd;  
 if (Eb(number_of_local_basis_trial_trace_u+1,current_edge)==i)   
        T_right=Eb(number_of_local_basis_trial_trace_u+2,current_edge);
       indexFlip(indexFaceV)=indexFaceV;
 else
        T_right=Eb(number_of_local_basis_trial_trace_u+1,current_edge);     
       indexFaceV1=indexFaceV;
        temp = indexFaceV1(1);    % ����һ��Ԫ�ش洢����ʱ������
        indexFaceV1(1) = indexFaceV1(end); % ���ڶ���Ԫ�ظ�ֵ����һ��Ԫ��
        indexFaceV1(end) = temp;    % ����ʱ�����е�ֵ��ֵ���ڶ���Ԫ��
         indexFlip(indexFaceV)=indexFaceV1;
%       indexFlip(indexFaceV) = flip(indexFaceV);
 end
    %��ȡ��ǰ�ߵĶ��㲢����ߵĳ��Ⱥͷ�������
    vertex1=Eb(1,current_edge);   %�ߵĵ�һ����������
    vertex2=Eb(number_of_local_basis_trial_trace_u,current_edge);    %�ߵĵڶ�����������
    vertex_left_index=T((1:4),T_left);   %�ñ���൥Ԫ�Ķ�������
    vertices_left=P(:,vertex_left_index);     % �ڱߵ���൥Ԫ�Ķ���
     normal=find_normal(vertex1,vertex2,vertex_left_index,P);   %�����˸����������һ����������ÿ�����ϵĲ�û�м��㣩
     [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(order_Gauss);
     vertices_ww1=[P(:,vertex1)';P(:,vertex2)'];
     [Gauss_nodes,Gauss_weights]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_ww1);
%      [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(P(:,vertex1),P(:,vertex2),order_Gauss);   
    vertices_ww=[P(:,vertex1),P(:,vertex2)];
 if(T_right~=-1)    %�����Ƿ��ڱ߽���

       for alpha=1:number_of_local_basis_trial 
          for beta=1:number_of_local_basis_test 
  
                int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
%                 int_value=int_value+Gauss_quadrature_1D_trial_test_u(data.a,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
                Auu(beta,alpha)=Auu(beta,alpha) +tau*int_value; % ͬһ��������Ĺ���

           end
       end
       
       for alpha=1:number_of_local_basis_trial_trace_u 
            for beta=1:number_of_local_basis_test 
                 int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u,alpha,0,basis_type_test,beta,0,0);
                Aul(beta,indexFaceV(alpha))=Aul(beta,indexFaceV(alpha)) +tau*int_value; % ͬһ��������Ĺ���
             end
       end
        for alpha=1:number_of_local_basis_trial_trace_u 
            for beta=1:number_of_local_basis_test 
                if normal(1)~=0
                 int_value=normal(1)*Gauss_quadrature_2D_volume_trial_test__trace_u_v(data.b1,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u,alpha,0,basis_type_test,beta,0,0);
                Aul1(beta,indexFaceV(alpha))=Aul1(beta,indexFaceV(alpha)) -int_value; % ͬһ��������Ĺ���
                else
                  int_value=normal(2)*Gauss_quadrature_2D_volume_trial_test__trace_u_v(data.b2,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u,alpha,0,basis_type_test,beta,0,0);
                Aul2(beta,indexFaceV(alpha))=Aul2(beta,indexFaceV(alpha)) -int_value; % ͬһ��������Ĺ���   
                end
             end
       end
       
       
       for alpha=1:number_of_local_basis_trial_trace_u
           for beta=1:number_of_local_basis_test 
      if normal(1)~=0
                int_value=normal(1)*Gauss_quadrature_2D_volume_trial_test__trace_u_q(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u,alpha,0,basis_type_test,beta,0,0);
               Aql1(beta,indexFaceV(alpha))=Aql1(beta,indexFaceV(alpha))+int_value; % minus because we have - integral.
      else
               int_value=normal(2)*Gauss_quadrature_2D_volume_trial_test__trace_u_q(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,vertices_left,basis_type_trace_u,alpha,0,basis_type_test,beta,0,0);
               Aql2(beta,indexFaceV(alpha))=Aql2(beta,indexFaceV(alpha))+int_value; % minus because we have - integral.
      end

           end

       end
  for alpha=1:number_of_local_basis_trial_trace_u
                for delta=1:number_of_local_basis_trial_trace_u
               int_value=Gauss_quadrature_2D_volume_trial_test__trace_u_trace_v(coe_fun,Gauss_nodes,Gauss_weights,vertices_ww,basis_type_trace_u,alpha,0,basis_type_trace_u,delta,0);
                All(indexFaceV(delta),indexFaceV(alpha))=All(indexFaceV(delta),indexFaceV(alpha)) -tau*int_value; % ͬһ��������Ĺ���
                end
  end
  for alpha=1:number_of_local_basis_trial_trace_u
                for delta=1:number_of_local_basis_trial_trace_u
                  if normal(1)~=0   
               int_value=normal(1)*Gauss_quadrature_2D_volume_trial_test__trace_u_trace_v(data.b1,Gauss_nodes,Gauss_weights,vertices_ww,basis_type_trace_u,alpha,0,basis_type_trace_u,delta,0);
                All1(indexFaceV(delta),indexFaceV(alpha))=All1(indexFaceV(delta),indexFaceV(alpha)) +int_value; % ͬһ��������Ĺ���
                  else
                int_value=normal(2)*Gauss_quadrature_2D_volume_trial_test__trace_u_trace_v(data.b2,Gauss_nodes,Gauss_weights,vertices_ww,basis_type_trace_u,alpha,0,basis_type_trace_u,delta,0);
                All2(indexFaceV(delta),indexFaceV(alpha))=All2(indexFaceV(delta),indexFaceV(alpha)) +int_value; % ͬһ��������Ĺ���
                  end
                end
  end
         elseif(T_right==-1)
         
         for alpha=1:number_of_local_basis_trial
            for beta=1:number_of_local_basis_test
                % Contribution of terms +\int \mu [[phi_j]][[\phi_i]]
                int_value=Gauss_quadrature_1D_trial_test_u(coe_fun,Gauss_nodes,Gauss_weights,vertices_left,vertices_left,basis_type_trial,alpha,0,0,basis_type_test,beta,0,0);
                Auu(beta,alpha)=Auu(beta,alpha) +tau*int_value; % Contribution from terms of the same triangle.

           end
         end
 end  
    % Global indexing
    indexFaceIni = indexFaceEnd + 1;
end

% Elemental mapping -------------------------------------------------------
% B=[Auu Auq1 Auq2; Auq1' Aqq1 Q; Auq2' Q Aqq2];
% b=[Aul fu; Aql1 fq;Aql2 fq];
Z = [Auu+Auu1+Auu2+Auu3 Auq1 Auq2; Auq1' Aqq1 Q; Auq2' Q Aqq2]\[Aul+Aul1+Aul2 fu; Aql1 fq; Aql2 fq];
vU = 1:ndofU;
vQ = ndofU+1:ndofU+ndofQ;
vP=ndofU+ndofQ+1: ndofU+2*ndofQ;
Zul = Z(vU,1:ndofUHat);
zuf = Z(vU,ndofUHat+1);
Zql1 = Z(vQ,1:ndofUHat);
Zql2 = Z(vP,1:ndofUHat);
zqf1 = Z(vQ,ndofUHat+1);
zqf2 = Z(vP,ndofUHat+1);

Alq1 = Aql1';
Alq2 = Aql2';
Alu = Aul';
% Elemental matrices
Ke = Alq1*Zql1+Alq2*Zql2 + Alu*Zul + All+All1+All2;
fe = fl - (Alq1*zqf1+Alq2*zqf2 + Alu*zuf);






