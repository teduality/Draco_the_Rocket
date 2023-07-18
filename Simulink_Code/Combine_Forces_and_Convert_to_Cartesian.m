function [Cart_Accel, Inertial_Accel] = Combine_Forces_and_Convert_to_Cartesian(XYZ_Inertial_Force, Impulse, Mass, DCM)

XYZ_Inertial_Force(3) = Impulse-XYZ_Inertial_Force(3);

XYZ_Inertial_Accel = (XYZ_Inertial_Force)/(Mass/1000);

Cart_Accel = (XYZ_Inertial_Accel'*DCM)';

Cart_Accel(3) = Cart_Accel(3) - 9.81;

Inertial_Accel = (Cart_Accel'*DCM')';
end