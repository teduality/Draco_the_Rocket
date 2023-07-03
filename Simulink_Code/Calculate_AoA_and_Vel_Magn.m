function [Vel_Magn, Inertial_Vel, AoA] = Calculate_AoA_and_Vel_Magn(Cart_Vel, DCM, Angles)

Vel_Magn = sqrt(Cart_Vel(1)^2+Cart_Vel(2)^2+Cart_Vel(3)^2);

Inertial_Vel = Cart_Vel*DCM';

AoA = [Angles(1) - unwrap(2*atan2(Inertial_Vel(1),Inertial_Vel(3)))/2; 
        Angles(2) - unwrap(2*atan2(Inertial_Vel(2),Inertial_Vel(3)))/2];
end