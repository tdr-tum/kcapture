function [equilibrium_ceiling] = capture_ceiling_kaolin(phi_H2O_KOH,TinK)
%capture_ceiling_kaolin, new fit function with KOH and H2O for Kaolin, T in K. 

TinK = min(max(TinK,773),2073);     % making sure it's not outside the T-range

logphi = log(phi_H2O_KOH);

%% Model Fit Parameters 
% Slag model: a*ln(phi)^2+b*ln(phi)+c*T+d

slagfit_a = 0.00078713829256592433847306011784895; 
slagfit_b = 0.028872204811572651367290731627691;
slagfit_c = -0.00018936911996406542409501772006308;
slagfit_d = 0.53845652199360249934301236862666;

% slagfit_phizero = exp(-slagfit_b/(2*slagfit_a));    % bottom of the function where d/dphi = 0.

ck_pred_slag = slagfit_a*logphi^2+slagfit_b*logphi+slagfit_c*TinK+slagfit_d;

%% First Step value
%  f(x) = p1*T^3 + p2*T^2 + p3*T + p4

firststep_p1 = 0.000000022635833943081541601832350641549;
firststep_p2 =  -0.00010181388485737569453339762537425;
firststep_p3 = 0.17043235418623378651936661754007;
firststep_p4 = -112.71757703646878212566662114114;

logphi_firststep = (firststep_p1*TinK^3 + firststep_p2*TinK^2 + firststep_p3*TinK + firststep_p4);

%% Second Step value
%  f(x) = p1*T^5 + p2*T^4 + p3*T^3 + p4*T^2 + p5*T + p6

secondstep_p1 = 0.00000000000030585550010907762617367149231012;
secondstep_p2 = -0.0000000016151975123924336386191202938768;
secondstep_p3 = 0.0000034127593207868893626416497560916;
secondstep_p4 =  -0.0036331316581537506302246320899485;
secondstep_p5 = 1.9917273982835677070823976464453;
secondstep_p6 = -482.53159511635249145911075174809;

logphi_secondstep = (secondstep_p1*TinK^5 + secondstep_p2*TinK^4 + secondstep_p3*TinK^3 + secondstep_p4*TinK^2 + secondstep_p5*TinK + secondstep_p6);

%% Third Step value
%  f(x) = p1*T^7 + p2*T^6 + p3*T^5 + p4*T^4 + p5*T^3 + p6*T^2 + p7*T + p8

thirdstep_p1 = 0.00000000000000000050773185184892504036144270137178;
thirdstep_p2 = -0.0000000000000045113225487202448014831506741288;
thirdstep_p3 = 0.000000000016943843554207104265615099427831;
thirdstep_p4 = -0.000000034872534425354132457511752455834;
thirdstep_p5 = 0.000042504539982575694184296216660002;
thirdstep_p6 = -0.030733929355301690955482030176427;
thirdstep_p7 = 12.266825769019444436480625881813;
thirdstep_p8 = -2123.5869726042783440789207816124;

logphi_thirdstep = (thirdstep_p1*TinK^7 + thirdstep_p2*TinK^6 + thirdstep_p3*TinK^5 + thirdstep_p4*TinK^4 + thirdstep_p5*TinK^3 + thirdstep_p6*TinK^2 + thirdstep_p7*TinK + thirdstep_p8);

%% Fourth Step value
%  f(x) = p1*T^5 + p2*T^4 + p3*T^3 + p4*T^2 + p5*T + p6

fourthstep_p1 = 0.000000000000022365587521137420158148209567174;
fourthstep_p2 = -0.0000000001511136868657877722147262319064;
fourthstep_p3 = 0.00000041521074761431613105839254218576;
fourthstep_p4 =  -0.00059339443624422945171292820276676;
fourthstep_p5 = 0.46371885336917761932795656321105;
fourthstep_p6 = -170.72352391532336923773982562125;

logphi_fourthstep = (fourthstep_p1*TinK^5 + fourthstep_p2*TinK^4 + fourthstep_p3*TinK^3 + fourthstep_p4*TinK^2 + fourthstep_p5*TinK + fourthstep_p6);

%% Hardcoded ceilings (from stoichiometry or graph)

ck_afterstep1 = 0.0545; % from graph
ck_afterstep2 = 0.0903; % from graph
ck_afterstep3 = 0.141; % from graph
ck_afterstep4 = 0.2626; % from graph

%% Getting the fit value
% heaviside step functions to switch between the plateaus

equilibrium_ceiling =   min(max(max(min( heaviside(logphi - logphi_firststep) * (ck_afterstep1),  (ck_pred_slag)),0) + ...
                        heaviside(logphi - min(logphi_secondstep,logphi_thirdstep) ) * (ck_afterstep2 - ck_afterstep1) + ...
                        heaviside(logphi - min(logphi_thirdstep,logphi_fourthstep)) * (ck_afterstep3 - ck_afterstep2) + ...
                        heaviside(logphi - logphi_fourthstep) * (ck_afterstep4 - ck_afterstep3), ...
                        heaviside(phi_H2O_KOH - 1e-5) * heaviside(TinK - 1372) * ck_pred_slag), ck_afterstep4);

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