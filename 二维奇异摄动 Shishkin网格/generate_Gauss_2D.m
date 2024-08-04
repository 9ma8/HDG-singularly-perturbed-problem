function [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,type)
%����'generate_Gauss_2D'���������㶨��ĵ�ǰ��Ԫ�ϴ��ݸ�˹�ڵ㡣


%��δ���Ĺ����Ǽ����ڲο���Ԫ�ϵĸ�˹���ֵ��Ȩ�أ�
%�����б�Ҫ������任���Ա���ʵ�ʵ�Ԫ�Ͻ��л��֡�



%=== ����ӳ�� 
A=[(vertices(1,3)-vertices(1,1))/2,0; 0,(vertices(2,3)-vertices(2,1))/2];
%���뿪ʼ������һ���任���� A�����ڽ��ο���Ԫ�ϵĸ�˹���ֵ�ӳ�䵽ʵ�ʵ�Ԫ�ϡ�
%����任���� A �Ǹ��ݵ�Ԫ�Ķ����������õ��ġ�
if type >=1  %���ݻ������� type �����˹���ֵ��Ȩ�أ�
    %��� type ���ڵ��� 1��������ڲο���Ԫ�ϵĸ�˹���ֵ��Ȩ�ء�����ʹ����һά���� ([-1, 1]) ��
    %�ĸ�˹-���õ»��ֽڵ��Ȩ�ؼ��㺯�� gauleg��
    w_2D=[];
    node_2D=[];
    % һά����[-1,1]�ϵĽڵ�
    n=type;
    [x,w]= gauleg(-1,1,n);
    for i=1:n
        for j=1:n
            node_2D=[node_2D; x(i), x(j)];
            w_2D = [w_2D,w(i).*w(j)];
        end
    end
    v=[(vertices(1,1)+vertices(1,3))/2; (vertices(2,1)+vertices(2,3))/2];
    v=repmat(v,1,size(node_2D,1));
    Gauss_nodes=node_2D;
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';  %������
    Gauss_weights=w_2D;
elseif type==-1  %��� type ���� -1����ֱ��ʹ�õ�Ԫ�Ķ�����Ϊ��˹���ֵ㣬��������ͬ��Ȩ�ء�
    Gauss_nodes=vertices';
    Gauss_weights=[1/3,1/3,1/3];
elseif type==-2    %��� type ���� -2��������ض��� P2b�����λ��������ĸ�˹���ֵ��Ȩ�أ�����һ������ض����ⶨ�ƵĻ��ַ�����
    %Lumped for P2b 
    Gauss_nodes=[0,0;1,0;0,1;0.5,0;0,0.5;0.5,0.5;1/3,1/3];
    v=[vertices(1,1); vertices(2,1)];
    v=[v,v,v,v,v,v,v];
    
    %�ڼ����˹���ֵ�󣬴�����ݱ任���� A ����Щ���ֵ��������任��������ӳ�䵽ʵ�ʵ�Ԫ�ϡ�
    %������Ϊ������Ԫ�����У�ͨ����Ҫ�ڲο���Ԫ�Ͻ��л��֣���ʵ�ʵ�Ԫ�ϵĻ�����Ҫͨ������任����ת����
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';
    Gauss_weights=[1/20,1/20,1/20,2/15,2/15,2/15,9/20];
end
    
 
%���붨����һ������ gauleg�����ڼ���һά�����ϵĸ�˹-���õ»��ֽڵ��Ȩ�ء�
function [x,w]= gauleg(a,b,n)

m=(n+1)/2;
xm=0.5*(b+a);
xl=0.5*(b-a);
xx=[];

for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
        if (abs(z-z1)<eps), break, end
    end
    xx(i)=xm-xl*z;
    xx(n+1-i)=xm+xl*z;
    ww(i)=2.0*xl/((1.0-z*z)*pp*pp);
    ww(n+1-i)=ww(i);
end

x=xx;
w=ww;
end

end
