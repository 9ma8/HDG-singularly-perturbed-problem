function [uhat, local] = hdg_Poisson_GlobalSystem(P,T,Tb,Eb,coe_fun,f_fun,basis_type_trace_u,basis_type_u,tau,order_Gauss,geo,data,varepsilon)



number_of_elements=size(T,2);      %单元个数
number_of_local_basis_trialll=size(Tb,1);     %在单个单元上trial局部基函数的个数
if basis_type_trace_u==104
    number_of_local_basis_trial=5;
% number_of_local_basis_test=5;
elseif basis_type_trace_u==103
   number_of_local_basis_trial=4;
% number_of_local_basis_test=4;
elseif basis_type_trace_u==102
    number_of_local_basis_trial=3;
%  number_of_local_basis_test=3;
elseif basis_type_trace_u==101
    
    number_of_local_basis_trial=2;
%  number_of_local_basis_test=2;
else
    % 这里可以处理所有其他情况或输出一条错误信息
    fprintf('Warning: basis_type_trace_u has an unsupported value (%d).\n', basis_type_trace_u);

end
% Initialisation
mat.i = zeros(1,2*number_of_local_basis_trial*geo.N*(geo.N+1));
mat.j = zeros(1,2*number_of_local_basis_trial*geo.N*(geo.N+1));
mat.Kij = zeros(1,2*number_of_local_basis_trial*geo.N*(geo.N+1));
f = zeros(2*number_of_local_basis_trial*geo.N*(geo.N+1),1);

local(number_of_elements).Zql1 = [];
local(number_of_elements).Zql2 = [];
local(number_of_elements).Zul = [];
local(number_of_elements).zqf1= [];
local(number_of_elements).zqf2 = [];
local(number_of_elements).zuf = [];

indexIni = 0;
for iElem = 1:number_of_elements
%     for j=1:4   %遍历当前三角形的边。对于每个边，它确定哪些三角形（T_left, T_right）与此边相邻
   vertices=P(:,T((1:4),iElem)); 
    [Zqle1,Zqle2,Zule,zqfe1,zqfe2,zufe,Ke,fe] = hdg_Poisson_ElementalMatrices(P,T,Tb,Eb,coe_fun,f_fun,basis_type_u,basis_type_u,basis_type_trace_u,tau,vertices,order_Gauss,iElem,data,varepsilon);
   
    
    % Store the local solver
%     tic;
    local(iElem).Zql1 = Zqle1;
    local(iElem).Zql2 = Zqle2;
    local(iElem).Zul = Zule;
    local(iElem).zqf1 = zqfe1;
     local(iElem).zqf2 = zqfe2;
    local(iElem).zuf = zufe;
    
%     elapsed_time = toc;

% fprintf('Execution time: %.4f seconds\n', elapsed_time);
    
    % Assembly
%     [indexGlobalFace, nOfFaceDOF] = hdgElemToFaceIndex(mesh.indexTf,refFace,hdg.faceInfo(iElem),refElem(pElem).nOfFaces);
nOfFaceDOF=4*number_of_local_basis_trial;   
index=Eb(number_of_local_basis_trial+3:end,T(number_of_local_basis_trialll+1:end,iElem));
indexGlobalFace=reshape(index,1,[]);
nOfFaceDOF2 = nOfFaceDOF^2;
    currentIndex = indexIni + (1:nOfFaceDOF);
    currentIndex2 = indexIni + (1:nOfFaceDOF2);
    for i=1:nOfFaceDOF
        mat.i(currentIndex) = indexGlobalFace(i);
        mat.j(currentIndex) = indexGlobalFace;
        currentIndex = currentIndex + nOfFaceDOF;
    end
    mat.Kij(currentIndex2) = reshape(Ke',1,nOfFaceDOF2);
    
    % RHS
    f(indexGlobalFace) = f(indexGlobalFace) + fe;
    
    % Update indices
    indexIni = indexIni + nOfFaceDOF2;
%     end
end

K = sparse(mat.i,mat.j,mat.Kij);   %使用行索引mat.i、列索引mat.j和值mat.Kij构建全局稀疏刚度矩阵K


% u=K\f;
% Solution
%这部分代码初始化未知数向量uhat，然后求解线性系统以找到未知数uhat。
%这里只求解那些需要解决的自由度，即hdg.vDOFtoSolv
selected_columns = find(Eb(number_of_local_basis_trial+2,:) ~= -1);
selected_values = Eb(number_of_local_basis_trial+3:end, selected_columns);
vDOFtoSolve = selected_values(:)';

% tic;
uhat = zeros(2*number_of_local_basis_trial*geo.N*(geo.N+1),1);
uhat(vDOFtoSolve) = K(vDOFtoSolve,vDOFtoSolve)\f(vDOFtoSolve);
% elapsed_time1 = toc;
% fprintf('Execution time: %.4f seconds\n', elapsed_time1);










