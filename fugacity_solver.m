function [fugacity_coef_vdw, fugacity_vdw, fugacity_coef_rk, fugacity_rk, fugacity_coef_preos, fugacity_preos] = calc_fugacity(Tc,Pc,T,P, omega)
    R = 8.134; % J/(mol*K)
    a_vdw = (27/64)*((R*Tc)^2)/(Pc);
    b_vdw = (R*Tc)/(8*(Pc));
    syms v
    assume(v,'real')
    vdw_eqn = (R*T)/(v-b_vdw) - (a_vdw/v^2);
    v_vdw = vpasolve(vdw_eqn == P, [-Inf Inf]);
    fugacity_coef_vdw = exp(-log((v_vdw-b_vdw)*(P)/(R*T)) + b_vdw/(v_vdw-b_vdw) - 2*a_vdw/(R*T*v_vdw));
    fugacity_vdw = fugacity_coef_vdw * P;
    
    a_rk = (0.42748*(R^2)*(Tc^(2.5)))/Pc;
    b_rk = 0.08664*R*Tc/Pc;
    syms v_rk
    assume(v_rk,'real')
    rk_eqn = (R*T)/(v_rk-b_rk) - (a_rk/((sqrt(T)*v_rk*(v_rk+b_rk))));
    v_rk = vpasolve(rk_eqn == P, [-Inf Inf]);
    z_rk = (P*v_rk)/(R*T);
    fugacity_coef_rk = exp(z_rk - 1 - log((v_rk-b_rk)*P/(R*T)) - ((a_rk/(b_rk*R*(T^1.5)))*log(1+(b_rk/v_rk))));
    fugacity_rk = fugacity_coef_rk*P;

    a_preos = 0.45624*(R^2)*(Tc^2)/Pc;
    b_preos = 0.07780*R*Tc/Pc;
    kappa = 0.37464 + 1.54226*omega - 0.26992*(omega^2);
    alpha = (1 + kappa*(1- sqrt(T/Tc)))^2;
    syms v_preos
    assume(v_preos,'real')
    preos_eqn = (R*T)/(v_preos-b_preos) - ((a_preos*alpha)/((v_preos*(v_preos+b_preos))+b_preos*(v_preos-b_preos)));
    v_preos = vpasolve(preos_eqn == P, [-Inf Inf]);
    z_preos = (P*v_preos)/(R*T);
    fugacity_coef_preos = exp(z_preos - 1 - log((v_preos-b_preos)*P/(R*T)) - (a_preos*alpha)/((2*sqrt(2)*b_preos*R*T))*log(((v_preos+(1+sqrt(2))*b_preos))/((v_preos+(1-sqrt(2))*b_preos))));
    fugacity_preos = fugacity_coef_preos*P;
end
Tc_prompt = "Input your compounds critical temperature in units of Kelvin: ";
input_Tc = input(Tc_prompt);
Pc_prompt = "Input your compounds critical pressure in units of bar: ";
input_Pc_bar = input(Pc_prompt);
P_prompt = "Input your system pressure: ";
input_P = input(P_prompt);
P_unit_prompt = "Input the units of your system pressure ('bar','atm', 'torr', 'Pa', 'kPa', or 'MPa'): ";
input_P_unit = input(P_unit_prompt, 's');
T_prompt = "Input your system temperature: ";
input_T = input(T_prompt);
T_unit_prompt = "Input the units of your system temperature ('C', 'K', 'T', or 'R'): ";
input_T_unit = input(T_unit_prompt, 's');
omega_prompt = "Input your compounds accentricity factor (omega): ";
input_omega = input(omega_prompt);
if strcmpi(input_T_unit, 'C')
    input_T = input_T + 273.15;  % Celsius to Kelvin
elseif strcmpi(input_T_unit, 'T')
    input_T = ((input_T - 32) * (5/9)) + 273.15; % Fahrenheit to Kelvin
elseif strcmpi(input_T_unit, 'R')
    input_T = input_T * (5/9); % Rankine to Kelvin
elseif strcmpi(input_T_unit, 'K')
    % No conversion needed for Kelvin
else
    error(['Invalid temperature unit. Please input either ''C'', ''K'', or ''R''.']);
end

% Convert pressure to bar if necessary
if strcmpi(input_P_unit, 'atm')
    input_P = input_P * 101325;  % atm to Pa
elseif strcmpi(input_P_unit, 'torr')
    input_P = input_P * 133.322;
elseif strcmpi(input_P_unit, 'bar')
    input_P = input_P*10^5;
elseif strcmpi(input_P_unit, 'kPa')
    input_P = input_P * 10^3;
elseif strcmpi(input_P_unit, 'MPa')
    input_P = input_P * 10^6;
elseif strcmpi(input_P_unit, 'Pa')
    input_P = input_P*1.0;
else
    error(['Invalid pressure unit. Please input either ''bar'', ''atm'', ''torr''' ...
        '''Pa'', ''kPa'', or ''MPa'' '])
end
input_Pc = input_Pc_bar * 10^5;
[vdw_coef, vdw_fugacity, rk_coef, rk_fugacity, preos_coef, preos_fugacity] = calc_fugacity(input_Tc,input_Pc,input_T,input_P,input_omega);

fprintf('Van der Waals:\n\tFugacity Coefficient: %.6f\n\tFugacity: %.6f Pa\nRedlich-Kwong:\n\tFugacity Coefficient: %.6f\n\tFugacity: %.6f Pa\nPeng-Robinson\n\tFugacity Coefficient: %.6f\n\tFugacity: %.6f Pa\n', vdw_coef, vdw_fugacity, rk_coef, rk_fugacity, preos_coef, preos_fugacity)