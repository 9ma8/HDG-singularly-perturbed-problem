                                                      
% % %%*-----------------------------------------------------------------------*
% % %%*
% % %%*  2020年5月3日 
% % %%*- ----------------------------------------------------------------------*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
N=[4 8 16 32 64 128 256]


%%Shishkin mesh k=1, CD=10^{-8} 
% error=[0.161E+0 0.787E-1 0.326E-1 0.120E-1 0.413E-2 0.135E-2 0.427E-3 0.136E-3]
error=[8.442e+0   1.549e+0   2.770e-1   5.115e-2   1.129e-2  3.639e-3   1.455e-3]
% error=[0.702E-1 0.180E-1 0.457E-2 0.115E-2 0.288E-3 0.721E-4 0.180E-4  0.448E-5]
h1=loglog(N,abs(error),'b:s');
hold on

%%Shishkin mesh k=2, CD=10^{-8} 
% error=[0.556E-1 0.193E-1 0.520E-2 0.117E-2 0.237E-3 0.445E-4 0.797E-5 0.663E-4]
% error=[0.216E-1 0.329E-2 0.516E-3 0.840E-4 0.141E-4 0.242E-5 0.421E-6 0.781E-7]
error=[4.389e-1  4.086e-2  5.857e-03   1.553e-03  4.534e-04  1.260e-4    3.282e-5]
h2=loglog(N,abs(error),'b--o');
hold on 

% %%Shishkin mesh k=3, CD=10^{-8} 
% % error=[0.164E-1 0.394E-2 0.683E-3 0.940E-4 0.112E-4 0.124E-5 0.186E-5 0.123E-4]\
% % error=[0.382E-2 0.246E-3 0.157E-4 0.100E-5 0.637E-7 0.169E-7 0.310E-7 0.619E-7]
% error=[0.382E-2 0.246E-3 0.157E-4 0.100E-5 0.637E-7 0.169E-7]
% h3=loglog(N,abs(error),'b--s');
% hold on 
% 

% %%Shishkin mesh k=4, CD=10^{-8}
% % error=[0.476E-2 0.811E-3 0.916E-4 0.775E-5 0.514E-6 0.120E-6 0.383E-5 0.183E-4]
% error=[0.476E-2 0.811E-3 0.916E-4 0.775E-5]
% h4=loglog(N,abs(error),'b:o');
% hold on
%     
% 
% 
% %%Shishkin mesh k=5, CD=10^{-8}
% % error=[0.136E-2 0.165E-3 0.123E-4 0.727E-6 0.211E-6 0.440E-6 0.195E-5 0.222E-4]
% error=[0.136E-2 0.165E-3 0.123E-4 0.727E-6]
% h5=loglog(N,abs(error),'b:*');
% hold on
    
%%Bakhvalov-type mesh k=1, CD=10^{-13} 
% error=[0.161E+0 0.786E-1 0.325E-1 0.120E-1 0.414E-2 0.135E-2 0.430E-3 0.152E-3]
% h6=loglog(N,abs(error),'b--*');
% hold on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%**   
%%**   
%%** 描绘参考曲线y=N^{-2}
%%** 
%%** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:7
%  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
%    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
cf_plot(i)=  N(i)^(-1.5)* (log( N(i) )^(2));
% cf_plot(i)= 10e-4* N(i)^( -1)+N(i)^(-2)
end 
r1=loglog(N,cf_plot,'r-s');
hold on 



for i=1:1:7
%  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
%    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
cf_plot(i)=  N(i)^(-2.5)* (log( N(i) )^(3));
% cf_plot(i)= 10e-4* N(i)^( -2)+N(i)^(-3)
end 
r2=loglog(N,cf_plot,'r-o');
hold on 

% for i=1:1:6
% %  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
% %    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
% cf_plot(i)=  N(i)^(-3.5)
% % cf_plot(i)= 10e-4* N(i)^( -2)+N(i)^(-3)
% end  
% r3=loglog(N,cf_plot,'r--s');
% hold on 

% for i=1:1:4
% %  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
% %    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
% cf_plot(i)=  N(i)^( -5)* (log( N(i) )^(5.5))
% % cf_plot(i)= 10e-4* N(i)^( -2)+N(i)^(-3)
% end  
% r4=loglog(N,cf_plot,'r--o');
% hold on 
% 
% for i=1:1:4
% %  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
% %    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
% cf_plot(i)=  N(i)^( -6)* (log( N(i) )^(6.5))
% % cf_plot(i)= 10e-4* N(i)^( -2)+N(i)^(-3)
% end  
% r5=loglog(N,cf_plot,'r-*');
% hold on 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N1=[8 16 32 64]
% %%Bakhvalov-Shishkin mesh k=4, CD=10^{-8}
% error1=[0.7286E-3 0.5709E-4 0.3979E-5 0.2758E-6]
% h5=loglog(N1,abs(error1),'b:d');
% hold on
%     
% %%Bakhvalov-type mesh k=4, CD=10^{-8} 
% error1=[0.1727E-2 0.8239E-4 0.4656E-5 0.3034E-6]
% h6=loglog(N1,abs(error1),'b--d');
% hold on 
% 
% for i=1:1:4
% %  cf_plot(i)= N(i)^(-1.5)* (log( N(i) )^(3/2));
% %    cf_plot(i)= N(i)^( -(p+1) )* (log( N(i) )^(p+1));
% cf_plot1(i)=  N(i)^( -4)
% % cf_plot(i)= 10e-4* N(i)^( -1)+N(i)^(-2)
% end 
% r3=loglog(N1,cf_plot1,'r-d');
% hold on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%**   
%%**   
%%** 添加各条曲线的legends 
%%** 
%%** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend([r1,r2,h1,h2],'\fontsize{14} N^{-1.5} (log(N(i))^2)','\fontsize{14} N^{-2.5} (log(N(i))^3)','\fontsize{14} k=1','\fontsize{14} k=2');
xlabel('N');
ylabel('Errors $|||\Pi_{1}q-q_{h}, \Pi_{2}u-u_{h}, Pu-\hat{u}_{h})|||$','Interpreter','Latex');
title('$|||(\Pi_{1}q-q_{h}, \Pi_{2}u-u_{h}, Pu-\hat{u}_{h})|||$ in the case of $\varepsilon = 10^{-6}$ on Shishkin mesh','Interpreter','Latex');
% title('$|||(\bm{q}-\bm{q}_{h}, u-u_{h}, u-\hat{u}_{h})|||$', 'Interpreter', 'Latex');
% title(['in the case of $\varepsilon = 10^{-6}$ on Shishkin mesh'], 'Interpreter', 'Latex');

 
clear all

% % %%*-----------------------------------------------------------------------*
% % %%*
% % %%*  2020年5月3日 
% % %%*-----------------------------------------------------------------------*
% clear all
% figure
% N=[8 16 32 64 128]
% 
% % %%Bakhvalov-type mesh k=4, CD=10^{-4} 
% error=[0.463E-3 0.1829E-4 0.9326E-6 0.5492E-7 0.339E-8]
% h1=loglog(N,abs(error),'b--o');
% hold on 
% 
% % %%Bakhvalov-type mesh k=4, CD=10^{-5} 
% error=[0.1594E-3 0.6187E-5 0.3033E-6 0.1750E-7 0.1072E-8]
% h2=loglog(N,abs(error),'b--+');
% hold on 
% 
% % %%Bakhvalov-type mesh k=4, CD=10^{-6} 
% error=[0.5366E-4 0.2046E-5 0.9787E-7 0.5571E-8 0.3396E-9]
% h3=loglog(N,abs(error),'b--*');
% hold on 
% 
% % %%Bakhvalov-type mesh k=4, CD=10^{-7} 
% error=[0.1774E-4 0.6679E-6 0.3141E-7 0.1770E-8 0.1075E-9]
% h4=loglog(N,abs(error),'b--s');
% hold on 
% 
% % %%Bakhvalov-type mesh k=4, CD=10^{-8} 
% error=[0.5799E-5 0.2164E-6 0.1005E-7 0.5619E-9 0.3406E-10]
% h5=loglog(N,abs(error),'b--d');
% hold on 
%  
% 
% for i=1:1:5
% cf_plot(i)= N(i)^( -4)
% % cf_plot(i)= 10e-4* N(i)^( -1)+N(i)^(-2)
% end 
% r1=loglog(N,cf_plot,'r-');
% hold on  
% 
% hh=legend([r1,h1,h2,h3,h4,h5],'$N^{-4 }$','$\varepsilon=10^{-4}$','$\varepsilon=10^{-5}$','$\varepsilon=10^{-6}$','$\varepsilon=10^{-7}$','$\varepsilon=10^{-8}$');
% set(hh,'Interpreter','Latex')
% xlabel('N');
% ylabel('Errors $\Vert u-u^N \Vert_{\varepsilon}$','Interpreter','Latex');
% clear all


%%%%%%%%%下述为备用数据
% % %%Bakhvalov-type mesh k=2, CD=10^{-4} 
% error=[0.302E-2 0.580E-3 0.132E-3 0.321E-4 0.799E-5 0.200E-5]
% h1=loglog(N,abs(error),'b-o');
% hold on 
% 
% % %%Bakhvalov-type mesh k=2, CD=10^{-5} 
% error=[0.982E-3 0.186E-3 0.419E-4 0.102E-4 0.253E-5 0.631E-6]
% h2=loglog(N,abs(error),'b-+');
% hold on 
% 
% % %%Bakhvalov-type mesh k=2, CD=10^{-6} 
% error=[0.315E-3 0.592E-4 0.133E-4 0.322E-5 0.799E-6 0.200E-6]
% h3=loglog(N,abs(error),'b-*');
% hold on 
% 
% % %%Bakhvalov-type mesh k=2, CD=10^{-7} 
% error=[0.101E-3 0.188E-4 0.421E-5 0.102E-5 0.253E-6 0.631E-7]
% h4=loglog(N,abs(error),'b-s');
% hold on 
% 
% % %%Bakhvalov-type mesh k=2, CD=10^{-8} 
% error=[0.322E-4 0.597E-5 0.133E-5 0.322E-6 0.799E-7 0.200E-7]
% h5=loglog(N,abs(error),'b-d');
% hold on 
%  
% 
% for i=1:1:6
% cf_plot(i)=  N(i)^( -2)
% % cf_plot(i)= 10e-4* N(i)^( -1)+N(i)^(-2)
% end 
% r1=loglog(N,cf_plot,'r-');
% hold on  
% 
% 
% % %%Bakhvalov-type mesh k=3, CD=10^{-4} 
% error=[0.114E-2 0.998E-4 0.108E-4 0.130E-5 0.161E-6 0.201E-7]
% h6=loglog(N,abs(error),'b--o');
% hold on 
% 
% % %%Bakhvalov-type mesh k=3, CD=10^{-5} 
% error=[0.379E-3 0.326E-4 0.347E-5 0.413E-6 0.510E-7 0.636E-8]
% h7=loglog(N,abs(error),'b-+');
% hold on 
% 
% % %%Bakhvalov-type mesh k=3, CD=10^{-6} 
% error=[0.124E-3 0.105E-4 0.110E-5 0.131E-6 0.161E-7 0.201E-8]
% h8=loglog(N,abs(error),'b-*');
% hold on 
% 
% % %%Bakhvalov-type mesh k=3, CD=10^{-7} 
% error=[0.420E-4 0.337E-5 0.351E-6 0.415E-7 0.511E-8 0.637E-9]
% h9=loglog(N,abs(error),'b-s');
% hold on 
% 
% % %%Bakhvalov-type mesh k=3, CD=10^{-8} 
% error=[0.130E-4 0.108E-5 0.112E-6 0.131E-7 0.162E-8 0.201E-9]
% h10=loglog(N,abs(error),'b-d');
% hold on 
%  
% 
% for i=1:1:6
% cf_plot(i)=  N(i)^( -3)
% % cf_plot(i)= 10e-4* N(i)^( -1)+N(i)^(-2)
% end 
% r1=loglog(N,cf_plot,'r-');
% hold on  