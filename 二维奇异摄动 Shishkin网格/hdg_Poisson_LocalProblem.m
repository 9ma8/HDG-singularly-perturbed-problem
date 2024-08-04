function [u,q1,q2] = hdg_Poisson_LocalProblem(Pb,T,Tb,Eb,uHat,local,basis_type_trace_u)

% Count the number of unknowns
nOfDOFU=size(Pb,2);
nOfElements=size(Tb,2);
number_of_local_basis_trial=size(Tb,1); 
% Initialisation
u = zeros(nOfDOFU,1);
q1 = zeros(nOfDOFU,1);   %初始化解向量 u 和通量向量 q
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
indexQ1Ini = 1;    %初始化解向量 u 和通量向量 q 的索引
indexQ2Ini = 1; 
for iElem = 1:nOfElements    % 再次循环遍历所有单元
%     pElem = mesh.pElem(iElem);    %获取当前单元的参考单元索引
    
    %根据单元的面信息和HDG方法相关信息计算全局面索引
    indexGlobalFace1 = Eb(number_of_local_basis_trial_trace_u+3:end,T(number_of_local_basis_trial+1:number_of_local_basis_trial+4,iElem));
    indexGlobalFace =  indexGlobalFace1(:)';
    %计算当前单元的解和通量的自由度数量
    nOfElementNodes = number_of_local_basis_trial;
    nDOFsElemU = nOfElementNodes;
    nDOFsElemQ1 = nDOFsElemU;
    nDOFsElemQ2 = nDOFsElemU;
    
    %计算当前单元解和通量向量的结束索引
    indexUEnd = indexUIni + nDOFsElemU - 1;
    indexQ1End = indexQ1Ini + nDOFsElemQ1 - 1;
    indexQ2End = indexQ2Ini + nDOFsElemQ2-1;
    
    %计算当前单元的解向量和通量向量，
    %其中 Zul 和 Zql 是局部刚度矩阵，uHat 是边界条件，zuf 和 zqf 是局部载荷向量
    u(indexUIni:indexUEnd) = local(iElem).Zul*uHat(indexGlobalFace) + local(iElem).zuf;
    q1(indexQ1Ini:indexQ1End) = local(iElem).Zql1*uHat(indexGlobalFace) + local(iElem).zqf1;
    q2(indexQ2Ini:indexQ2End) = local(iElem).Zql2*uHat(indexGlobalFace) + local(iElem).zqf2;
    
    %更新解和通量向量的起始索引
    indexUIni = indexUEnd + 1;
    indexQ1Ini = indexQ1End + 1;
    indexQ2Ini = indexQ2End + 1;
end