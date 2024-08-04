function [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number)

%Generate the Gauss coefficients and Gauss points on the reference interval [-1,1].
%Gauss_point_number:the number of Gauss points in the formula. The Gauss formula depends on it.
%Gauss_coefficient_reference_1D,Gauss_point_reference_1D: the Gauss coefficients and Gauss points on the reference interval.

if Gauss_point_number==4
    Gauss_coefficient_reference_1D=[0.3478548451,0.3478548451,0.6521451549,0.6521451549];
    Gauss_point_reference_1D=[0.8611363116,-0.8611363116,0.3399810436,-0.3399810436];
elseif Gauss_point_number==3
    Gauss_coefficient_reference_1D=[0.5555555556,0.8888888889,0.5555555556];
    Gauss_point_reference_1D=[0.7745966692,0,-0.7745966692];
elseif Gauss_point_number==5
    Gauss_coefficient_reference_1D=[0.5688888889,0.4786286705,0.4786286705,0.2369268851,0.2369268851];
    Gauss_point_reference_1D=[0,0.5384693101,-0.5384693101,0.9061798459,-0.9061798459];
elseif Gauss_point_number==6
    Gauss_coefficient_reference_1D=[0.4679139346,0.4679139346,0.3607615730,0.3607615730,0.1713244924,0.1713244924];
    Gauss_point_reference_1D=[0.2386191861,-0.2386191861,0.6612093865,-0.6612093865,0.9324695142,-0.9324695142];
elseif Gauss_point_number==8
    Gauss_coefficient_reference_1D=[0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834];
    Gauss_point_reference_1D=[0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425];
elseif Gauss_point_number==2
    Gauss_coefficient_reference_1D=[1,1];
    Gauss_point_reference_1D=[-1/sqrt(3),1/sqrt(3)];
end