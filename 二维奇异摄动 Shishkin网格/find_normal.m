function normal=find_normal(vertex1,vertex2,vertex_left_index,P)
% ����find_normal�����������������������һ����������������������
%���ҵ�������ȱʧ����������3��
%Ȼ����㵽�ߵ������������������������btw�ͷ����������������ķ��š�
%-----------------------------------------------------------

%��Ŀ���Ǹ��ݸ����Ķ���������εĶ������������������εķ�����

normal=zeros(2,1);     %��ʼ��һ��2x1�������� normal�����ڴ洢����õ��ķ�����

% ���ݸ����Ķ������� vertex1 �� vertex2 �Լ������εĶ������� vertex_left_index��
% ȷ�������εĵ������������� vertex3
if(vertex1==vertex_left_index(1) && vertex2==vertex_left_index(2) || vertex2==vertex_left_index(1) && vertex1==vertex_left_index(2))
         vertex3=vertex_left_index(3);
    elseif(vertex1==vertex_left_index(2) && vertex2==vertex_left_index(3) || vertex2==vertex_left_index(2) && vertex1==vertex_left_index(3)) 
         vertex3=vertex_left_index(4);
    elseif(vertex1==vertex_left_index(3) && vertex2==vertex_left_index(4) || vertex2==vertex_left_index(3) && vertex1==vertex_left_index(4)) 
         vertex3=vertex_left_index(1);
    elseif(vertex1==vertex_left_index(4) && vertex2==vertex_left_index(1) || vertex2==vertex_left_index(4) && vertex1==vertex_left_index(1)) 
         vertex3=vertex_left_index(2);
end

% ���㶥�� vertex1 ������ vertex2 ������ vintern
%����������֮��������� tau���������һ�����õ���λ������
    tau=P(:,vertex2)-P(:,vertex1);    %������ʸ��
    tau=tau/norm(tau);      %��һ������������λ����
    %���õ�λ������ tau ���㷨��������x����Ϊ������y�����ĸ�ֵ��
    %y����Ϊ������x������ֵ���Ӷ��õ�����������ֱ�ķ�����
    
    normal(1)=-tau(2);
    normal(2)=tau(1);
    
    %���㶥�� vertex1 ������ vertex3 ������ vintern��
    %�������һ�����õ���λ����
    vintern=P(:,vertex3)-P(:,vertex1);
    vintern=vintern/norm(vintern);
    
    
    %������� vintern �ͷ����� normal �ĵ���Ƿ�Ϊ��ֵ��
    %����ǣ��򽫷�����ȡ������ȷ��������ָ����ȷ�ķ���
    if vintern'*normal >0
        normal=-normal;
    end
end

%��δ���������Ǽ����������������εĶ����������ɵ������εķ�������
%��ȷ���ⷨ�����ķ�����ȷ


