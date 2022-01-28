function [equilibrium_ceiling] = capture_ceiling_CFA(phi_H2O_KOH,TinK)
%capture_ceiling_CFA New fit function with KOH and H2O for coal fly ash, T
%in K. This one is a bit more complex than the kaolin one due to higher impurities. 

TinK = min(max(TinK,773),2073);     % making sure it's not outside the T-range

logphi = log(phi_H2O_KOH);

%% First Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

firststep_p1 = 0.000000054615827538531813749026524442906;
firststep_p2 =  -0.00019740621473815381323874129648743;
firststep_p3 = 0.26427760353353318967251084359305;
firststep_p4 = -142.97149894904464417777489870787;

logphi_firststep = (firststep_p1*TinK^3 + firststep_p2*TinK^2 + firststep_p3*TinK + firststep_p4);

%% Second Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

secondstep_p1 = 0.000000028270839537495622985630299318222;
secondstep_p2 =  -0.00012101917980949013383101114937901;
secondstep_p3 = 0.19281543612387025410370711142605;
secondstep_p4 = -120.61398859884849343870882876217;

logphi_secondstep = (secondstep_p1*TinK^3 + secondstep_p2*TinK^2 + secondstep_p3*TinK + secondstep_p4);

%% Third Step value
%  f(x) = p1*T^5 + p2*T^4 + p3*T^3 + p4*T^2 + p5*T + p6

thirdstep_p1 = 0.0000000000015464229238253555274888500476537;
thirdstep_p2 = -0.0000000083837940390443407254879870931126;
thirdstep_p3 = 0.000017979237250463776771154320033297;
thirdstep_p4 = -0.01908416062303403340494334372579;
thirdstep_p5 = 10.070820086542463300816052651498;
thirdstep_p6 = -2149.0799174794105965702328830957;

logphi_thirdstep = thirdstep_p1*TinK^5 + thirdstep_p2*TinK^4 + thirdstep_p3*TinK^3 + thirdstep_p4*TinK^2 + thirdstep_p5*TinK + thirdstep_p6;

%% Fourth Step value
% linear, because third order would have a negative p1

logphi_fourthstep = 0.03994 * TinK - 58.83;

%% Fifth Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

fifthstep_p1 = 0.000000026647702266242045367618834664052;
fifthstep_p2 =  -0.00010682806921626427233262129901092;
fifthstep_p3 = 0.16080818678252062348121853574412;
fifthstep_p4 = -99.030101554522659057511191349477;

logphi_fifthstep = (fifthstep_p1*TinK^3 + fifthstep_p2*TinK^2 + fifthstep_p3*TinK + fifthstep_p4);

%% Sixth Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

sixthstep_p1 = 0.000000011936618450002664765068061125675;
sixthstep_p2 =  -0.000055956337213561115387881750304899;
sixthstep_p3 = 0.10684392433048719994737041361077;
sixthstep_p4 =  -80.490341091955343699737568385899;

logphi_sixthstep = (sixthstep_p1*TinK^3 + sixthstep_p2*TinK^2 + sixthstep_p3*TinK + sixthstep_p4);

%% Seventh Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

seventhstep_p1 = 0.000000022176414423781928264672768844577;
seventhstep_p2 =  -0.00010623122335248450398993008692017;
seventhstep_p3 = 0.18174445719047493952302829711698;
seventhstep_p4 =  -112.79980134763265198216686258093;

logphi_seventhstep = (seventhstep_p1*TinK^3 + seventhstep_p2*TinK^2 + seventhstep_p3*TinK + seventhstep_p4);

%% Eighth Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

eighthstep_p1 = 0.0000000210516635028017922041590611823;
eighthstep_p2 =  -0.000095567769844479941057473204590877;
eighthstep_p3 = 0.1578561815233889276122170031158;
eighthstep_p4 =  -96.438471870390543472240096889436;

logphi_eighthstep = (eighthstep_p1*TinK^3 + eighthstep_p2*TinK^2 + eighthstep_p3*TinK + eighthstep_p4);

%% Slag Fit
% f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p21*x^2*y + p12*x*y^2 + p03*y^3
                    
slagfit_cfa_p00 =       1.377;
slagfit_cfa_p10 =  -0.0009797;
slagfit_cfa_p01 =      0.3022;
slagfit_cfa_p20 =   1.757e-07;
slagfit_cfa_p11 =  -0.0002233;
slagfit_cfa_p02 =     0.01755;
slagfit_cfa_p21 =   4.471e-08;
slagfit_cfa_p12 =  -6.955e-06;
slagfit_cfa_p03 =   0.0003096;

slagfit_cfa = heaviside(TinK-1370)*max(slagfit_cfa_p00 + slagfit_cfa_p10*TinK + slagfit_cfa_p01*logphi + slagfit_cfa_p20*TinK^2 + ...
        	  slagfit_cfa_p11*TinK*logphi + slagfit_cfa_p02*logphi^2 + slagfit_cfa_p21*TinK^2*logphi + ...
              slagfit_cfa_p12*TinK*logphi^2 + slagfit_cfa_p03*logphi^3, 0);

%% Hardcoded ceilings (from stoichiometry or graph)

ck_afterstep1 = 0.0508; % from graph
ck_afterstep2 = 0.063; % from graph
ck_afterstep3 = 0.068; % from graph
ck_afterstep4 = min(max(-3.5e-5*TinK+0.10356,ck_afterstep3),0.08); % from graph
ck_afterstep5 = min(max(1.71e-5*TinK+0.08019,0.09),0.1); % from graph
ck_afterstep6 = 0.1155; % from graph
ck_afterstep7 = 0.1295; % from graph
ck_afterstep8 = 0.194; % from graph

%% Getting the fit value
% heaviside step functions to switch between the plateaus

equilibrium_ceiling =   min(max(heaviside(logphi - min([logphi_firststep,logphi_secondstep])) * ck_afterstep1 + ...
                        heaviside(logphi - min([logphi_secondstep,logphi_thirdstep])) * (ck_afterstep2 - ck_afterstep1) + ...
                        heaviside(logphi - min([logphi_thirdstep,logphi_fourthstep,logphi_fifthstep])) * (ck_afterstep3 - ck_afterstep2) + ...
                        heaviside(logphi - min([logphi_fourthstep,logphi_fifthstep,logphi_sixthstep])) * (ck_afterstep4 - ck_afterstep3) + ...
                        heaviside(logphi - min([logphi_fifthstep,logphi_sixthstep,logphi_seventhstep])) * (ck_afterstep5 - ck_afterstep4) + ...
                        heaviside(logphi - min(logphi_sixthstep,logphi_seventhstep)) * (ck_afterstep6 - ck_afterstep5) + ...
                        heaviside(logphi - min(logphi_seventhstep,logphi_eighthstep)) * (ck_afterstep7 - ck_afterstep6) + ...
                        heaviside(logphi - logphi_eighthstep) * (ck_afterstep8 - ck_afterstep7), slagfit_cfa),  ...
                        ck_afterstep8);

end

%% Getting the fit data 
% resultslagslope = [];
% for T = 8:1:size(T_vec,2)   % set only to data with the slag curve
%     resultslagslope = vertcat(resultslagslope,[repmat(T_vec(T),[size(withSlagA.Cor.KOH_H2O(1:end,T),1),1]),withSlagA.Cor.KOH_H2O(1:end,T),withSlagA.result.mol.cK_cor_withSlagA(1:end,T)]);
% end
% POSTPROCESS THE DATA HERE TO ONLY INCLUDE THE RELEVANT CURVE
% resultslagslopeT = resultslagslope(:,1);
% resultslagslopeCk = resultslagslope(:,3);
% resultslagslopePhi = log(resultslagslope(:,2));