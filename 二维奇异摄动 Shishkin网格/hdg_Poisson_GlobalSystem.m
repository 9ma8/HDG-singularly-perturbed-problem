function [uhat, local] = hdg_Poisson_GlobalSystem(P,T,Tb,Eb,coe_fun,f_fun,basis_type_trace_u,basis_type_u,tau,order_Gauss,geo,data,varepsilon)



number_of_elements=size(T,2);      %��Ԫ����
number_of_local_basis_trialll=size(Tb,1);     %�ڵ�����Ԫ��trial�ֲ��������ĸ���
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
    % ������Դ�������������������һ��������Ϣ
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
%     for j=1:4   %������ǰ�����εıߡ�����ÿ���ߣ���ȷ����Щ�����Σ�T_left, T_right����˱�����
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

K = sparse(mat.i,mat.j,mat.Kij);   %ʹ��������mat.i��������mat.j��ֵmat.Kij����ȫ��ϡ��նȾ���K


% u=K\f;
% Solution
%�ⲿ�ִ����ʼ��δ֪������uhat��Ȼ���������ϵͳ���ҵ�δ֪��uhat��
%����ֻ�����Щ��Ҫ��������ɶȣ���hdg.vDOFtoSolv
selected_columns = find(Eb(number_of_local_basis_trial+2,:) ~= -1);
selected_values = Eb(number_of_local_basis_trial+3:end, selected_columns);
vDOFtoSolve = selected_values(:)';

% tic;
uhat = zeros(2*number_of_local_basis_trial*geo.N*(geo.N+1),1);
uhat(vDOFtoSolve) = K(vDOFtoSolve,vDOFtoSolve)\f(vDOFtoSolve);
% elapsed_time1 = toc;
% fprintf('Execution time: %.4f seconds\n', elapsed_time1);










