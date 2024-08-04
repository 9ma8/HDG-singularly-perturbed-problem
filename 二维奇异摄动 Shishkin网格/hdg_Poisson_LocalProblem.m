function [u,q1,q2] = hdg_Poisson_LocalProblem(Pb,T,Tb,Eb,uHat,local,basis_type_trace_u)

% Count the number of unknowns
nOfDOFU=size(Pb,2);
nOfElements=size(Tb,2);
number_of_local_basis_trial=size(Tb,1); 
% Initialisation
u = zeros(nOfDOFU,1);
q1 = zeros(nOfDOFU,1);   %��ʼ�������� u ��ͨ������ q
q2 = zeros(nOfDOFU,1);
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
indexUIni = 1; 
indexQ1Ini = 1;    %��ʼ�������� u ��ͨ������ q ������
indexQ2Ini = 1; 
for iElem = 1:nOfElements    % �ٴ�ѭ���������е�Ԫ
%     pElem = mesh.pElem(iElem);    %��ȡ��ǰ��Ԫ�Ĳο���Ԫ����
    
    %���ݵ�Ԫ������Ϣ��HDG���������Ϣ����ȫ��������
    indexGlobalFace1 = Eb(number_of_local_basis_trial_trace_u+3:end,T(number_of_local_basis_trial+1:number_of_local_basis_trial+4,iElem));
    indexGlobalFace =  indexGlobalFace1(:)';
    %���㵱ǰ��Ԫ�Ľ��ͨ�������ɶ�����
    nOfElementNodes = number_of_local_basis_trial;
    nDOFsElemU = nOfElementNodes;
    nDOFsElemQ1 = nDOFsElemU;
    nDOFsElemQ2 = nDOFsElemU;
    
    %���㵱ǰ��Ԫ���ͨ�������Ľ�������
    indexUEnd = indexUIni + nDOFsElemU - 1;
    indexQ1End = indexQ1Ini + nDOFsElemQ1 - 1;
    indexQ2End = indexQ2Ini + nDOFsElemQ2-1;
    
    %���㵱ǰ��Ԫ�Ľ�������ͨ��������
    %���� Zul �� Zql �Ǿֲ��նȾ���uHat �Ǳ߽�������zuf �� zqf �Ǿֲ��غ�����
    u(indexUIni:indexUEnd) = local(iElem).Zul*uHat(indexGlobalFace) + local(iElem).zuf;
    q1(indexQ1Ini:indexQ1End) = local(iElem).Zql1*uHat(indexGlobalFace) + local(iElem).zqf1;
    q2(indexQ2Ini:indexQ2End) = local(iElem).Zql2*uHat(indexGlobalFace) + local(iElem).zqf2;
    
    %���½��ͨ����������ʼ����
    indexUIni = indexUEnd + 1;
    indexQ1Ini = indexQ1End + 1;
    indexQ2Ini = indexQ2End + 1;
end