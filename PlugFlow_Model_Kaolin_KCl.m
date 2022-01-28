function [CK,comp1_mat,comp2_mat,comp3_mat,comp4_mat,comp5_mat]...
    = PlugFlow_Model_Kaolin_KCl(steps,KSpecies_in,EtOH_in,H2O_in,O2_in,N2_in,Temp,...
    Kaolin_in,residenceTime,A_input,E_a_input,n_input,m_input,species)

if floor(steps) == steps && steps > 0 %CHECKING IF STEPS IS A POSITIVE INT
    args{1} = steps;
else
    error('steps is not a positive integer')
    return
end

MWA_Temps = load('MWA_Temps.mat').MWA_Temps;

%% define the reactor first / 1D = ixo surface of the reactor
reactorLength = 2; % in m

%% define input values

kin.A = A_input; kin.Ea = E_a_input;
kin.n = n_input; kin.m = m_input;
% A_start, E_A_start, n_start, m_start

%% Adding the KOH adsorption
% latest values
kinKOH.A = 11917.0341; kinKOH.Ea =55896.15;
kinKOH.n = -2.62866; kinKOH.m = 1.3749;

%%
% A_start = pre-exponential factor
% E_a_start = activation energy
% n_start = Temperature exponent
% m_start = reaction order

%Arrhenius for surface area degredation
Sp.A = 389;
Sp.E_a = 94540.1;
R = 8.314;
surface_area_init = 12700;  % m^2/kg

%% Some useful species quantities
%Molar Masses in g/mol
mw.K = 39.1;
mw.H2O = 18; mw.KOH = 56.1;
mw.KCL = 74.6; mw.HCL = 36.5;
mw.K2SO4 =174.3; mw.SO3 = 80.1;
%% Cantera Calculation
%% Get phase from Chemkin file
solu = Solution('Mortensen-plus-Ethanol.xml','gas');

%% Parameters
%Defining the temporal change according to the position in the reactor
dt = residenceTime/steps; %in s
totalMassInput = KSpecies_in + EtOH_in + H2O_in + O2_in + N2_in; %total Mass as Input but neglected Kaolin

%% Matrix containing the compositions
% data like time position etc.

comp1_mat = array2table([[0:steps]',...
    linspace(0,residenceTime,steps+1)',...
    linspace(0,reactorLength,steps+1)',...
    arrayfun(@(x) TemperatureInsideReactor(MWA_Temps, Temp, x),linspace(0, 2, steps+1))',...
    Temp*ones(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    [surface_area_init;zeros(steps,1)],...
    [surface_area_init;zeros(steps,1)],...
    [totalMassInput;zeros(steps,1)],...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),zeros(steps+1,1),zeros(steps+1,1),...
    zeros(steps+1,1),zeros(steps+1,1),zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),...
    zeros(steps+1,1),zeros(steps+1,1),zeros(steps+1,1),...
    zeros(steps+1,1),zeros(steps+1,1)],...
    'VariableNames',{'steps','resTime','xPosition','Temp','aimedTemp','reactorTime','reactorTemp',...
    'Sp_FE','Sp_BE','GasMass','dMass','m_Kgas','m_Kadd',...
    'molfracKSpeciesStepEnd','molfracKOH','molfracKCL','molfracHCL','molfracKO2','molfracK','molfracO2',...
    'corK_O2','CKmax','rr','dCK','dm_K','CK','CK_fromKCL','CK_fromKOH','shareOfKCl'});
% totalMass, GasMass, and MassCaptured take into account the Gasphase
% Material (no Kaolin) after each Adsorption step

for i=1:steps
    % Next specific Surface from Forward Euler
    comp1_mat.Sp_FE(i+1) = comp1_mat.Sp_FE(i)*(1 - Sp.A*exp(-Sp.E_a/...
        (R*comp1_mat.Temp(i)))*dt); %Area Degradation

    % Next specific Surface from  Backward Euler
    comp1_mat.Sp_BE(i+1) = comp1_mat.Sp_BE(i) / (1 + Sp.A*exp(-Sp.E_a/...
        (R*comp1_mat.Temp(i+1)))*dt); %Area Degradation

end
%% Mass after reactor step
comp2_mat = array2table(zeros(steps+1,nSpecies(solu)),...
    'VariableNames',speciesNames(solu));
comp2_mat.KCL(1) = KSpecies_in;
comp2_mat.O2(1) = O2_in;
comp2_mat.N2(1) = N2_in;
comp2_mat.H2O(1) = H2O_in;
comp2_mat.C2H5OH(1) = EtOH_in;

%% Mass fractions after reactor step
comp3_mat = array2table(zeros(steps+1,nSpecies(solu)),...
    'VariableNames',speciesNames(solu));
comp3_mat.KCL(1) = KSpecies_in/totalMassInput;
comp3_mat.O2(1) = O2_in/totalMassInput;
comp3_mat.N2(1) = N2_in/totalMassInput;
comp3_mat.H2O(1) = H2O_in/totalMassInput;
comp3_mat.C2H5OH(1) = EtOH_in/totalMassInput;

%% Mass after adsorption step
comp4_mat = array2table(zeros(steps+1,nSpecies(solu)),...
    'VariableNames',speciesNames(solu));
comp4_mat.KCL(1) = KSpecies_in;
comp4_mat.O2(1) = O2_in;
comp4_mat.N2(1) = N2_in;
comp4_mat.H2O(1) = H2O_in;
comp4_mat.C2H5OH(1) = EtOH_in;

%% Mass fractions after adsorption step
comp5_mat = array2table(zeros(steps+1,nSpecies(solu)),...
    'VariableNames',speciesNames(solu));
comp5_mat.KCL(1) = KSpecies_in/totalMassInput;
comp5_mat.O2(1) = O2_in/totalMassInput;
comp5_mat.N2(1) = N2_in/totalMassInput;
comp5_mat.H2O(1) = H2O_in/totalMassInput;
comp5_mat.C2H5OH(1) = EtOH_in/totalMassInput;

%% Set chem-Solution at the inlet
% This writes the new massfractions into a format importable for a
% solution. It automaticlally adapts to the number of species

Y = reshape([char(speciesNames(solu)'),repmat(':',nSpecies(solu),1),num2str(comp5_mat{1,:}'),repmat(' ',nSpecies(solu),1)]',1,[]);
set(solu, 'Y', Y);
set(solu, 'T', comp1_mat.Temp(1));
set(solu, 'P', 101325);
comp1_mat.molfracKSpeciesStepEnd(1) = moleFraction(solu, 'KOH') +moleFraction(solu, 'K') +moleFraction(solu, 'KO2') +moleFraction(solu, 'KCL');
comp1_mat.molfracO2(1) = moleFraction(solu, 'O2') ;
comp1_mat.molfracK(1) = moleFraction(solu, 'K');
comp1_mat.molfracKCL(1) = moleFraction(solu, 'KCL');

comp1_mat.m_Kgas(1) =  massFraction(solu, 'KCL')*comp1_mat.GasMass(1)*(mw.K/mw.KOH);


%% Looping through the whole reactor with Cantera
for i = 1:steps

    % activate for kinetics mode
    % time step in reactor
    reac = IdealGasConstPressureReactor(solu);
    network = ReactorNet({reac});
    advance(network, dt);

    comp1_mat.reactorTime(i+1)= time(network);
    comp1_mat.reactorTemp(i+1)= temperature(solu);

    % save reactor results
    comp3_mat{i+1,:}= massFractions(reac);

    %         % activate for equilibrium mode 
    %         equSolu= equilibrate(solu, 'TP');
    %         comp3_mat{i+1,:}= massFractions(equSolu)';
    %         comp1_mat.reactorTemp(i+1)= temperature(solu);

    comp2_mat{i+1,:}=comp3_mat{i+1,:}*comp1_mat.GasMass(i);

    % compute capture ceiling
    comp1_mat.corK_O2(i+1) = power(moleFraction(solu,'O2')*pressure(solu),0.25) *moleFraction(solu,'K')*pressure(solu);
    comp1_mat.corKOH_H2O(i+1) = moleFraction(solu,'KOH')*pressure(solu) *power(moleFraction(solu,'H2O')*pressure(solu),-0.5);

    comp1_mat.CKmax(i+1) = capture_ceiling_kaolin(comp1_mat.corKOH_H2O(i+1),comp1_mat.Temp(i)); % in [kg/kg]
    
    %% Integrating the capture equation numerically

    %%%% Euler Semi-Backwards v2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This should be the solution in case you include the function of
    % p(t+1) = f(CK(t+1)), while leaving CK_max, m_gas and M_gas
    % explicit to avoid expensive re-calculation with cantera/the fit. This
    % saves you from taking extremely small time steps

    meanmolarweight = density(solu)/molarDensity(solu);

    % First KOH
    lambda1_KOH = dt*kinKOH.A*exp(-kinKOH.Ea/(R*comp1_mat.Temp(i+1)))*power(comp1_mat.Temp(i+1),kinKOH.n)*comp1_mat.Sp_BE(i+1)*power(Kaolin_in*pressure(solu)*meanmolarweight/(comp1_mat.GasMass(i)*mw.KOH),kinKOH.m);
    lambda2_KOH = comp4_mat.KOH(i)/Kaolin_in+comp1_mat.CK(i);

    % solving with fzero for speed
    if comp1_mat.CKmax(i+1) > 0 % check whether max capture is zero.
        eulerResult_KOH = fzero(@(testereuler) comp1_mat.CK(i) + lambda1_KOH * power(max(0,(lambda2_KOH - testereuler)),kinKOH.m)*(1-testereuler/comp1_mat.CKmax(i+1))-testereuler, comp1_mat.CK(i));
    else
        eulerResult_KOH = 0;
    end

    %         % iteratively solving with a fixed point method
    %         eulerResult_KOH = comp1_mat.CK(i); % starting value
    %         error = 1;  % starting value
    %         URF_KOH = 0.4;  % under-relaxation-factor, 1 = fast, unstable, 0 = slow, stable
    %         while error > 1e-9      % performing the iterative loop
    %             eulerOld = eulerResult_KOH; % value before iteration
    %             if comp1_mat.CKmax(i+1) == 0 % check whether max capture is zero.
    %                 break
    %             end
    %             eulerResult_KOH = (1-URF_KOH) * eulerResult_KOH + URF_KOH * (comp1_mat.CK(i) + lambda1_KOH * power(max(0,(lambda2_KOH - eulerResult_KOH)),kinKOH.m)*(1-eulerResult_KOH/comp1_mat.CKmax(i+1))); % iteration
    %             if eulerResult_KOH == 0
    %                 break
    %             end
    %             error = abs((eulerResult_KOH-eulerOld)/eulerResult_KOH);     % arbitrary error definition
    %         end

    % Next KCl
    lambda1_KCL = dt*kin.A*exp(-kin.Ea/(R*comp1_mat.Temp(i+1)))*power(comp1_mat.Temp(i+1),kin.n)*comp1_mat.Sp_BE(i+1)*power(Kaolin_in*pressure(solu)*meanmolarweight/(comp1_mat.GasMass(i)*mw.KCL),kin.m);
    lambda2_KCL = comp4_mat.KCL(i)/Kaolin_in+comp1_mat.CK(i);

    % solving with fzero for speed
    if comp1_mat.CKmax(i+1) > 0 % check whether max capture is zero.
        eulerResult_KCL = fzero(@(testereuler) comp1_mat.CK(i) + lambda1_KCL * power(max(0,(lambda2_KCL - testereuler)),kin.m)*(1-testereuler/comp1_mat.CKmax(i+1))-testereuler, comp1_mat.CK(i));
    else
        eulerResult_KCL = 0;
    end

    %         % iteratively solving with a fixed point method
    %         eulerResult_KCL = comp1_mat.CK(i); % starting value
    %         error = 1;  % starting value
    %         URF_KCL = 0.4;  % under-relaxation-factor, 1 = fast, unstable, 0 = slow, stable
    %         while error > 1e-9      % performing the iterative loop
    %             eulerOld = eulerResult_KCL; % value before iteration
    %             if comp1_mat.CKmax(i+1) == 0 % check whether max capture is zero.
    %                 break
    %             end
    %             eulerResult_KCL = (1-URF_KCL) * eulerResult_KCL + URF_KCL * (comp1_mat.CK(i) + lambda1_KCL * power(max(0,(lambda2_KCL - eulerResult_KCL)),kin.m)*(1-eulerResult_KCL/comp1_mat.CKmax(i+1))); % iteration
    %             if eulerResult_KCL == 0
    %                 break
    %             end
    %             error = abs((eulerResult_KCL-eulerOld)/eulerResult_KCL);     % arbitrary error definition
    %         end

    deltaCk_KOH = max(0,eulerResult_KOH - comp1_mat.CK(i));         % preventing release of K
    deltaCk_KCL = max(0,eulerResult_KCL - comp1_mat.CK(i));         % preventing release of K
    deltaCk_tot = deltaCk_KCL+deltaCk_KOH;

    comp1_mat.CK(i+1) = comp1_mat.CK(i) + deltaCk_tot;
    % Calculating the change in capture value within the step
    comp1_mat.dCK(i+1) = comp1_mat.CK(i+1)-comp1_mat.CK(i); % in kg_K/kg_Add

    comp1_mat.rr(i+1) = comp1_mat.dCK(i+1) / dt;
    % Calculating the mass of K consumed by the additive within the step

    comp1_mat.dm_K(i+1) = comp1_mat.dCK(i+1)*Kaolin_in; %mass of K consumed in kg_K
    comp1_mat.m_Kgas(i+1) = comp1_mat.m_Kgas(i)- comp1_mat.dm_K(i+1);
    comp1_mat.m_Kadd(i+1) = comp1_mat.m_Kadd(i)+comp1_mat.dm_K(i+1);
   
    %% Adding changes to the gas phase due to capture
    % Kaolin + 2KCl + H2O --> K-Kaolin + 2HCl
    comp4_mat{i+1,:} = comp2_mat{i+1,:};
    comp4_mat.KCL(i+1) = comp4_mat.KCL(i+1) - (eulerResult_KCL-comp1_mat.CK(i)) * Kaolin_in * (mw.KCL/mw.K);
    comp4_mat.HCL(i+1) = comp4_mat.HCL(i+1) + (eulerResult_KCL-comp1_mat.CK(i)) * Kaolin_in * (mw.HCL/mw.K);
    comp4_mat.H2O(i+1) = comp4_mat.H2O(i+1) - (eulerResult_KCL-comp1_mat.CK(i)) * Kaolin_in * (mw.H2O/mw.K) * 0.5;

    % Kaolin + 2KOH --> K-Kaolin + H20
    comp4_mat.KOH(i+1) = comp4_mat.KOH(i+1) - (eulerResult_KOH-comp1_mat.CK(i)) * Kaolin_in * (mw.KOH/mw.K);
    comp4_mat.H2O(i+1) = comp4_mat.H2O(i+1) + (eulerResult_KOH-comp1_mat.CK(i)) * Kaolin_in * (mw.H2O/mw.K) * 0.5;

    comp1_mat.GasMass(i+1) = sum(comp4_mat{i+1,:});
    comp1_mat.dMass(i+1) = comp1_mat.GasMass(i+1)-comp1_mat.GasMass(i);
    comp5_mat{i+1,:}= comp4_mat{i+1,:}/comp1_mat.GasMass(i+1);

    %% Setting some monitors 
    comp1_mat.molfracO2(i+1) = moleFraction(solu, 'O2') ;
    comp1_mat.molfracK(i+1) = moleFraction(solu, 'K');
    comp1_mat.molfracKO2(i+1) = moleFraction(solu, 'KO2');
    comp1_mat.molfracKOH(i+1) = moleFraction(solu, 'KOH');
    comp1_mat.molfracKCL(i+1) = moleFraction(solu, 'KCL');
    comp1_mat.molfracHCL(i+1) = moleFraction(solu, 'HCL');

    %% Set Chem-Solution for next step
    Y = reshape([char(speciesNames(solu)'),repmat(':',nSpecies(solu),1),num2str(comp5_mat{i+1,:}'),repmat(' ',nSpecies(solu),1)]',1,[]);
    set(solu, 'Y', Y);
    set(solu, 'T', comp1_mat.Temp(i+1));
    set(solu, 'P', 101325);

    % Molefrac of K-Species at the end of step
    comp1_mat.molfracKSpeciesStepEnd(i+1) = moleFraction(solu, 'KOH') +moleFraction(solu, 'K') +moleFraction(solu, 'KO2')+moleFraction(solu, 'KCL');

end

% Return Capture at Outlet
CK = comp1_mat.CK(end);
end



