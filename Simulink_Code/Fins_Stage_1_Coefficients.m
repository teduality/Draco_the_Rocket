function [C_N_f1, C_P_f1, C_A_f1] = Fins_Stage_1_Coefficients(d_bt, AoA, Re)
n = 2; %number of fins in a direction
dn = d_bt; %diameter at base of nosecone m
df = d_bt; %diameter of body tube at fins m
ls = 0.05; %fin span m
lr = .1; %fin root chord m
lt = .06; %fin tip-chord m
lm = lt + (lr-lt)./2; %fin mid chord m
Tf = 0.003175; %fin thickness m
lTS = 2.*ls + df; %fin total span
K_fb = 1 + ((df./2)./(ls+(df./2)));

Re_C = 5.*10.^5; %critical reynolds number

%Normal Coefficients

C_N_f1 = ((AoA).*K_fb.*(2.*n.*(ls./dn).^2)./(1+sqrt(1+((2.*lm)./(lr+lt)).^2))); %coefficient of stability

%Center of Pressure

%Xf = .9; %tip of nosecone to fin leading edge

%C_P_fin = Xf + (lm.*(lr+2.*lt))./(3.*(lr+lt)) + (1./6).*(lr+lt-(lr.*lt)./(lr+lt)); %center of pressure

C_P_f1 = 1.451;

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

C_D_f1 = (C_D_0 + C_D_i + C_D_AoA); %Total drag coefficient

%Axial Coefficient
C_A_f1 = ((C_D_f1.*cos(AoA)-(1./2).*C_N_f1.*sin(2.*AoA))./(1-sin(AoA).^2));
end