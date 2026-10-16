function [theta_l, theta_u, delta_st_l, delta_st_u] = extension_thwaites_turbulent_v2(nu, s_l, s_u, c_l, c_u, theta_l_0, theta_u_0, delta_st_l_0, delta_st_u_0) 
     
    % Lower surface resolution
    [theta_l, delta_st_l] = solve_turbulent_side(nu, s_l, c_l, theta_l_0, delta_st_l_0); 
     
    % Upper surface resolution
    [theta_u, delta_st_u] = solve_turbulent_side(nu, s_u, c_u, theta_u_0, delta_st_u_0); 
     
    function [theta_sol, delta_st_sol] = solve_turbulent_side(nu_val, s, Ue, th_0, dst_0) 
 
        % Regression model constants 
        Cc = 1.45;         
        CRe_inf = 0.0024;  
        Cm = 7.23;         
 
        s = s(:); Ue = Ue(:); 
 
        % 1. External velocity gradient
        dUe_ds = gradient(Ue, s); 
 
        % 2. Interpolation functions for ode45
        Ue_fun = @(sq) interp1(s, Ue, sq, 'makima'); 
        dUe_ds_fun = @(sq) interp1(s, dUe_ds, sq, 'makima'); 
 
        % 3. Differential Equation for theta
        ode_theta = @(sq, th) (CRe_inf / 2) ... 
                            - (Cm / 2) * (th / (Ue_fun(sq) + eps)) * dUe_ds_fun(sq) ... 
                            + (Cc / 2) * (nu_val / (Ue_fun(sq) * th + eps)); 
 
        % 4. ODE Resolution
        [~, theta_sol] = ode45(ode_theta, s, th_0); 
        theta_sol = theta_sol(:);  
 
        Ue_safe = max(Ue, 1e-3); 
        dUe_ds_lim = max(min(dUe_ds, 50), -50);  
 
        % Ideal turbulent shape factor H_eq
        Re_theta = (Ue_safe .* theta_sol) ./ nu_val; 
        Gamma = - (theta_sol ./ Ue_safe) .* dUe_ds_lim .* (Re_theta.^0.25); 
        H_eq = 1.4 + 0.9 * max(Gamma, 0); 
        H_eq = min(H_eq, 2.8);  
 
        delta_st_sol = zeros(size(s)); 
        delta_st_sol(1) = dst_0; 
 
        % Zero Pressure Gradient (ZPG) minimum growth rate
        % Divided by 8 to convert delta to delta* (0.05 offset avoids s=0 singularity)
        Re_s = (Ue_safe .* (s - s(1) + 0.05)) ./ nu_val;  
        ddeltastar_zpg_ds = (0.296 / 8) .* (Re_s .^ -0.2);  
 
        for i = 2:length(s) 
            ds = s(i) - s(i-1); 
 
            % Ideal steady-state target
            dst_target = theta_sol(i) * H_eq(i); 
 
            % Minimum physical growth guaranteed by ZPG term
            min_growth = ddeltastar_zpg_ds(i) * ds; 
 
            % ZPG Transition
            if dst_target <= delta_st_sol(i-1) 
                delta_st_sol(i) = delta_st_sol(i-1) + min_growth; 
            else 
                relax = 0.15; % 15% target, 85% previous inertia 
                delta_st_sol(i) = relax * dst_target + (1 - relax) * (delta_st_sol(i-1) + min_growth); 
            end 
        end 
    end 
 
end