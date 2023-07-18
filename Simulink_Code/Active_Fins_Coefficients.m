function [C_N_flaps, C_N_roll, C_P_flaps, C_A_flaps] = Active_Fins_Coefficients(d_bt, AoA, AoA_flaps, Re)
n = 1; %number of flaps
dn = d_bt; %diameter at base of nosecone
df = d_bt; %diameter of body tube at fins
ls = 0.05; %fin span
lr = .04; %fin root chord
lt = .04; %fin tip-chord
lm = lt + (lr-lt)./2; %fin mid chord
Tf = 0.003175; %fin thickness m
lTS = 2.*ls + df; %fin total span
K_fb = 1 + ((df./2)./(ls+(df./2)));

Re_C = 5.*10.^5; %critical reynolds number

%Normal Coefficients

Dimension_Calc = K_fb.*(4.*n.*(ls./dn).^2)./(1+sqrt(1+((2.*lm)./(lr+lt)).^2));

C_N_flap_1 = [((AoA(1)+AoA_flaps(1)).*Dimension_Calc);...
              0]; %coefficient of stability

C_N_flap_2 = [0;...
              ((AoA(2)+AoA_flaps(2)).*Dimension_Calc)]; %coefficient of stability

C_N_flap_3 = [((AoA(1)+AoA_flaps(3)).*Dimension_Calc);...
              0]; %coefficient of stability

C_N_flap_4 = [0;...
              ((AoA(2)+AoA_flaps(4)).*Dimension_Calc)]; %coefficient of stability

C_N_flaps =  C_N_flap_1 + C_N_flap_2 + C_N_flap_3 + C_N_flap_4;

C_N_roll = C_N_flap_1(1) + C_N_flap_1(2) + C_N_flap_2(1) + C_N_flap_2(2) - C_N_flap_3(1) - C_N_flap_3(2) - C_N_flap_4(1) - C_N_flap_4(2); 

%Center of Pressure

Xf = 1.2; %tip of nosecone to fin leading edge

C_P_flaps = Xf + (lm.*(lr+2.*lt))./(3.*(lr+lt)) + (1./6).*(lr+lt-(lr.*lt)./(lr+lt)); %center of pressure

%Drag Coefficients

A_fe = (1./2).*(lr+lt).*ls; %m.^2fin exposed area

A_fp = A_fe + (1./2).*df.*lr; %m.^2 fin planform area

if Re < 1
    C_f = 0;
elseif Re < Re_C
    C_f = 1.328./(sqrt(Re)); 
else
    B = Re_C.*((0.074./(Re.^(1./5)))-(1.328./(sqrt(Re))));
    C_f = (0.074./(Re.^(1./5)))-(B./Re);
end

C_D_0 = 2.*C_f.*(1+2.*(Tf./lm)).*(4.*n.*A_fp); %coefficient of fin drag zero angle of attack

C_D_i = 2.*C_f.*(1+2.*Tf./lm).*(4.*n.*(A_fp-A_fe))./(pi.*df.^2); %interference drag

Rs = lTS./df;

k_fb = 0.8065.*Rs.^2 + 1.1553.*Rs; %fin-body interference coefficient

k_bf = 0.1935.*Rs.^2 + 0.8174.*Rs + 1; %body-fin inteference coefficient

C_D_AoA = AoA.^2.*(1.2.*(A_fp.*4)./(pi.*df.^2)+3.12.*(k_fb+k_bf-1).*((A_fe.*4)./(pi.*df.^2))); %Coefficient of alpha drag

C_D_flaps = C_D_0 + C_D_i + C_D_AoA; %Total drag coefficient

%Axial Coefficient

C_A_flaps = (C_D_flaps.*cos(AoA)-(1./2).*C_N_flaps.*sin(2.*AoA))./(1-sin(AoA).^2);
end