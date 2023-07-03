function [Angles_real, change] = If_on_rail_dont_change(XYZ_Coords, Angles)
Dist_magn = sqrt(XYZ_Coords(1)^2+XYZ_Coords(2)^2+XYZ_Coords(3)^2);
if Dist_magn > 10
    Angles_real = Angles;
    change=1;
else
    Angles_real = [deg2rad(5); deg2rad(0); deg2rad(0)];
    change=0;
end
end