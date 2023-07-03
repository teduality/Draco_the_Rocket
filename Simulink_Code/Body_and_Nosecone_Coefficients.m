function [NC_Base_A, C_N_B, C_N_NC, C_A_BN, C_P_NC, C_P_B] = Body_and_Nosecone_Coefficients(Re, L_N, L_B, D_B, D_D, L_C, AoA, Nu, K)
%{
L_N = Length of Nosecone
L_B = Length of Bodytube
L_C = Length of Conical Change
D_B = Diameter of Bodytube
D_D = Diameter of Endtip
Re = Reynolds Number
AoA = angle of attack
%} 

Re_C = 5.*10.^5; %critical reynolds number

if Re < Re_C
    C_f = 1.328./(sqrt(Re));
else
    B = Re_C.*((0.074./(Re.^(1./5)))-(1.328./(sqrt(Re))));
    C_f = (0.074./(Re.^(1./5)))-(B./Re);
end

L_TR = L_N + L_B;
C_D_fb = (1 + (60./((L_TR./D_B).^3)) + 0.0025.*(L_B./D_B))...
        .*(2.7.*(L_N./D_B) + 4.*(L_B./D_B) + 2.*(1-(D_D./D_B)).*(L_C./D_B)).*C_f;
%C_D full body
C_D_b = 0.029.*((D_D./D_B).^3)./sqrt(C_D_fb); %C_D base


delta = 0.9; %estimate from graphs later

C_D_AoA = 2..*delta.*(AoA.^2)...
           + ((3.6.*Nu.*(1.36.*L_TR-0.55.*L_N))./(pi.*D_B)).*(AoA.^3); %Additional drag at angle of attack

C_D_B = (C_D_fb + C_D_b + C_D_AoA);

%Coefficients Normal

B_Plan_A = D_B.*L_B; %body planform area m.^2

NC_Plan_A = (2./3).*D_B.*L_N; %Nosecone planform area m.^2

NC_Base_A = (pi./4).*D_B.^2; %Area base of nosecone, reference geometry

C_N_B = (K.*abs(AoA).*AoA.*(B_Plan_A./NC_Base_A)); %Coefficient for Body

C_N_NC = (K.*abs(AoA).*AoA.*(NC_Plan_A./NC_Base_A)); %Coefficient for Nosecone

%Coefficients Axial

C_N = C_N_NC + C_N_B; %total coefficient normal of nosecone and body

C_A_BN = ((C_D_B.*cos(AoA)-(1./2).*C_N.*sin(2.*AoA))./(1-sin(AoA).^2));

%Centers of Pressure

C_P_NC = 0.466.*L_N; %Nosecone (L>6R)

C_P_B = 0.5.*L_B+L_N; %Bodytube
end