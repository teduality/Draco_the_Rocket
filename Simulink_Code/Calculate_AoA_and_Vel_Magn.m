function [Vel_Magn_Rel_Wind, AoA] = Calculate_AoA_and_Vel_Magn(Wind, Cart_Vel, DCM, Angles)

Inertial_Vel = (Cart_Vel'*DCM')';

Inertial_Vel_Rel_Wind = [Inertial_Vel(1)+Wind(1);Inertial_Vel(2)+Wind(2);Inertial_Vel(3)+Wind(3)];

Vel_Magn_Rel_Wind = sqrt(Inertial_Vel_Rel_Wind(1)^2+Inertial_Vel_Rel_Wind(2)^2+Inertial_Vel_Rel_Wind(3)^2);

AoA = [Angles(1) - unwrap(2*atan2(Inertial_Vel_Rel_Wind(1),Inertial_Vel_Rel_Wind(3)))/2;
    Angles(2) - unwrap(2*atan2(Inertial_Vel_Rel_Wind(2),Inertial_Vel_Rel_Wind(3)))/2];
end