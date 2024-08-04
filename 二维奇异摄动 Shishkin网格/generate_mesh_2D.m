%-------------------------------------------------------------------------
%generate_mesh_2D�ڼ��κ���geo.funָ���ļ����������ɾ����������ߴ�geo.h������
% generate_mesh_2D����
% geo:������������ͼ�κ������Сgeo.h�ĺ�������Լ�����/�ǹ�������ı�־�Ľṹ��
% basis_type: FEM�ռ������(���ԡ�����Ԫ��)��
% generate_mesh_2D����:
% P:�������񶥵�λ�õĶ����������:P (: 1) = (X_1; Y_1)
% T:������Ԫ������������ϵ�Ԫk����������ĵ�Ԫ����
% Eg: T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain]�����ʹ��DG����T������Ҳ�����������αߵ�������
%E���߾��󡣵�һ�к͵ڶ��а����ߵ������յ�������������к͵����а��������յ����ֵ��
%�����а����߶α�ţ������к͵����а��������Ҳ������������

%Pb�����ɶȵĶ���������ʹ��P1, Pb=P����

%Tb����ÿ����Ԫ�Ͼ��ж������ɶȵĵ�Ԫ�������ʹ��P1, Tb=T��

%Eb���߾�����Eb�У���һ�к͵ڶ��а�����ʼ�ͽ�����������к͵����и�����ʱ�뷽���������Ҳ����������

%h���������������ڱߵ���С���ȣ�����?

%-------------------------------------------------------------------------


function [P,T,Pb,Tb,Eb]=generate_mesh_2D(geo,basis_type,basis_type_trace_u,data,varepsilon)


   N=geo.N;
   

%====DG P1 FE�����ݽṹ
%�����Eb��һ�к͵ڶ��а�����ʼ����ͽ��������������
%�����к͵����и�����ʱ�뷽���������Ҳ����������

tau1=min(1/2,data.sigma*varepsilon*log(N)/data.beta1);
tau2=min(1/2,data.sigma*varepsilon*log(N)/data.beta2);
%�����T��7�С�3������3����
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
    for k=1:size(T,2)    % �������е������ε�Ԫ
          for i=1:4      %��ÿ�������ε����������������ӵ� Pb ��
              j=j+1;
             Pb(:, j)=P(:,T(i,k));
          end
          Tb((1:4),k)=[j-3;j-2;j-1;j];     %���� Tb��ָ��ÿ�������ε����������Ӧ�����ɶ�
    end
    
    
    %% ============= BUILD MATRIX E_B
    
 %Eb ��һ���������ڴ洢����������ı߽���Ϣ��
% ͨ������ÿ�������Σ����߽���Ϣ�洢�� Eb �С�
% ����ÿ�������Σ����Ƚ�����������������洢�� Eb ��ǰ�����У������������߱��Ϊ 1��
% Ȼ�󣬼��ÿ�������ε��������Ƿ��Ѿ��� Eb �С�������ڣ��ͽ��ߵ���Ϣ��ӵ� Eb �У��������������бߵ�������
    
    
max_num_edges= 2*size(P,2)- 6;     %Ԥ�����㹻��ľ���Eb�Ĺ���
num_edges=1;
Eb=-ones(4,max_num_edges);    %������ֱ�������εĳ�ʼֵ�趨Ϊ�ߵġ��߽硱������������������Ǹ���ֵ��������ӡT_right
%��������������Ԥ�����㹻��ľ���ռ䣬����ʼ�� Eb

vertex1=T(1,1);     %��һ�������Σ��������бߡ�
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
T(7,1)=3;   %�ھ���T�в���ߵı�ǩ��
T(8,1)=4;
for k=2:size(T,2)    
    vertex1=T(1,k);     %��ȡ��ǰ�����ε����������������
    vertex2=T(2,k);
    vertex3=T(3,k);
    vertex4=T(4,k);
    m=num_edges;
    flag_1=0;
    flag_2=0;
    flag_3=0;
    flag_4=0;
for j=1:m    %���ÿ�����Ƿ��Ѿ������� Eb ��
    if( (vertex1==Eb(2,j) && vertex2==Eb(1,j)) || (vertex1==Eb(1,j) && vertex2==Eb(2,j)) )    %��鵱ǰ���Ƿ��� Eb �е�ĳ������ƥ��
        Eb(4,j)=k;    %��������������ΪEb����k
        T(5,k)=j;     % �ѱ߷��������ξ����С�  %����ǰ�߱��Ϊ���������Σ����� T �м�¼�ߵ�����
        flag_1=1;
    elseif ( (vertex2==Eb(2,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(2,j)) )    %���ڶ������Ƿ��Ѿ��������ˡ�
        Eb(4,j)=k;
        T(6,k)=j;
        flag_2=1;
    elseif ( (vertex3==Eb(2,j) && vertex4==Eb(1,j)) || (vertex3==Eb(1,j) && vertex4==Eb(2,j)))    %�����������Ƿ��Ѿ��������ˡ�
        Eb(4,j)=k;
        T(7,k)=j;
        flag_3=1;
     elseif ( (vertex4==Eb(2,j) && vertex1==Eb(1,j)) || (vertex4==Eb(1,j) && vertex1==Eb(2,j)))    %�����������Ƿ��Ѿ��������ˡ�
        Eb(4,j)=k;
        T(8,k)=j;
        flag_4=1;
    end
end
if(flag_1==0)     %�����ǰ�߲������� Eb �У��������
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
Eb=Eb(:,(1:num_edges));    %�����ǰ�߲������� Eb �У��������
% ��������ȫΪ 0 �ľ���
zeros_rows = zeros(2, size(Eb, 2)); % ��һ��ά������ӵ��������ڶ���ά���� Eb ������

% ��ȫΪ 0 ������ӵ����� Eb ��ĩβ
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
    T=[T; zeros(4,size(T,2))];  %��ԭʼ���� T ��ĩβ������������������ڴ洢ÿ�������ε������ߵ�������
    j=0;
    for k=1:size(T,2)    % �������е������ε�Ԫ
          for i=1:9      %��ÿ�������ε����������������ӵ� Pb ��
              j=j+1;
              Pb(:,j)=P(:,T(i,k));
          end
          Tb((1:9),k)=[j-8;j-7;j-6;j-5;j-4;j-3;j-2;j-1;j];     %���� Tb��ָ��ÿ�������ε����������Ӧ�����ɶ�
    end
        
    
    %% ============= BUILD MATRIX E_B
    
max_num_edges= 2*size(P,2)- 6;     %Ԥ�����㹻��ľ���Eb�Ĺ���
num_edges=1;
Eb=-ones(5,max_num_edges);    %������ֱ�������εĳ�ʼֵ�趨Ϊ�ߵġ��߽硱������������������Ǹ���ֵ��������ӡT_right
%��������������Ԥ�����㹻��ľ���ռ䣬����ʼ�� Eb
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
vertex1=T(1,1);     %��һ�������Σ��������бߡ�
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
T(12,1)=3;   %�ھ���T�в���ߵı�ǩ��
T(13,1)=4;
for k=2:size(T,2)    
    vertex1=T(1,k);     %��ȡ��ǰ�����ε����������������
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
for j=1:m    %���ÿ�����Ƿ��Ѿ������� Eb ��
    if( (vertex1==Eb(number_of_local_basis_trial_trace_u,j) && vertex2==Eb(1,j)) || (vertex1==Eb(1,j) && vertex2==Eb(number_of_local_basis_trial_trace_u,j)) )    %��鵱ǰ���Ƿ��� Eb �е�ĳ������ƥ��
        Eb(5,j)=k;    %��������������ΪEb����k
        T(10,k)=j;     % �ѱ߷��������ξ����С�  %����ǰ�߱��Ϊ���������Σ����� T �м�¼�ߵ�����
        flag_1=1;
    elseif ( (vertex2==Eb(number_of_local_basis_trial_trace_u,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(number_of_local_basis_trial_trace_u,j)) )    %���ڶ������Ƿ��Ѿ��������ˡ�
        Eb(5,j)=k;
        T(11,k)=j;
        flag_2=1;
    elseif ( (vertex3==Eb(number_of_local_basis_trial_trace_u,j) && vertex4==Eb(1,j)) || (vertex3==Eb(1,j) && vertex4==Eb(number_of_local_basis_trial_trace_u,j)))    %�����������Ƿ��Ѿ��������ˡ�
        Eb(5,j)=k;
        T(12,k)=j;
        flag_3=1;
     elseif ( (vertex4==Eb(number_of_local_basis_trial_trace_u,j) && vertex1==Eb(1,j)) || (vertex4==Eb(1,j) && vertex1==Eb(number_of_local_basis_trial_trace_u,j)))    %�����������Ƿ��Ѿ��������ˡ�
        Eb(5,j)=k;
        T(13,k)=j;
        flag_4=1;
    end
end
if(flag_1==0)     %�����ǰ�߲������� Eb �У��������
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
Eb=Eb(:,(1:num_edges));    %�����ǰ�߲������� Eb �У��������
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
