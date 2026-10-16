function [theta_l, theta_u, delta_st_l, delta_st_u] = extension_thwaites_turbulent_v1(nu, s_l, s_u, c_l, c_u, theta_l_0, theta_u_0, delta_st_l_0, delta_st_u_0)
    % =====================================================================
    % MODELLO THWAITES TURBOLENTO (Agrawal et al., 2024)
    % =====================================================================
    
    % Risoluzione Ventre (Lower)
    [theta_l, delta_st_l] = solve_turbulent_side(nu, s_l, c_l, theta_l_0, delta_st_l_0);
    
    % Risoluzione Dorso (Upper)
    [theta_u, delta_st_u] = solve_turbulent_side(nu, s_u, c_u, theta_u_0, delta_st_u_0);

    function [theta_sol, delta_st_sol] = solve_turbulent_side(nu_val, s, Ue, th_0, dst_0)
        % =====================================================================
        % MODELLO DI STRATO LIMITE TURBOLENTO - AGRAWAL ET AL. (2024)
        % =====================================================================

        % Coefficienti fittati dagli autori (Sezione 6.1, Paper Principale)
        Cc = 1.45;        
        CRe_inf = 0.0024; 
        Cm = 7.23;        

        % Assicuriamo che gli array siano vettori colonna
        s = s(:); 
        Ue = Ue(:);

        % 1. Calcolo del gradiente di velocità esterna
        dUe_ds = gradient(Ue, s);

        % Funzioni interpolanti per i sub-step numerici del solutore ode45
        Ue_fun = @(sq) interp1(s, Ue, sq, 'makima');
        dUe_ds_fun = @(sq) interp1(s, dUe_ds, sq, 'makima');

        % =====================================================================
        % FASE 1: Calcolo dello spessore di quantità di moto (theta)
        % Basato sull'Eq. 4.5 (e 4.3) del Paper Principale
        % =====================================================================

        % Definizione dell'Equazione Differenziale Ordinaria per theta
        ode_theta = @(sq, th) (CRe_inf / 2) ...
                            - (Cm / 2) * (th / (Ue_fun(sq) + eps)) * dUe_ds_fun(sq) ...
                            + (Cc / 2) * (nu_val / (Ue_fun(sq) * th + eps));

        % Risoluzione dell'ODE (ode45 restituisce la soluzione nei punti esatti 's')
        [~, theta_sol] = ode45(ode_theta, s, th_0);
        theta_sol = theta_sol(:); % Sicurezza colonna

        % =====================================================================
        % FASE 2: Calcolo dello spessore di spostamento (delta_star)
        % Basato sulle Eq. 4 e Eq. 6 del Materiale Supplementare
        % =====================================================================

        % Parametro di Holstein-Bohlen (m) per strato limite sottile
        m_sol = -(theta_sol.^2 ./ nu_val) .* dUe_ds;

        % Derivata dello spessore di quantità di moto al quadrato rispetto a s
        theta2_sol = theta_sol.^2;
        dtheta2_ds = gradient(theta2_sol, s);

        % Calcolo del termine di correzione ZPG (Zero Pressure Gradient)
        % Usiamo la legge classica turbolenta (aggiungiamo un piccolo offset 
        % a 's' per evitare singolarità fisiche all'inizio del dominio)
        Re_s = (Ue .* (s - s(1) + 1e-4)) ./ nu_val; 
        ddelta_zpg_ds = 0.296 .* (Re_s .^ -0.2);

        % Eq. 6 (Supplemento): Derivata dello spessore fisico dello strato limite (delta)
        ddelta_ds = (5/200) + (8/200) * m_sol - (Ue ./ (200 * nu_val)) .* dtheta2_ds + ddelta_zpg_ds;

        % Calcolo dell'integrale contenuto nell'Eq. 4
        % Integrale da s0 a s di: (d(delta)/dr * Ue)
        integrand = ddelta_ds .* Ue;
        integral_term = cumtrapz(s, integrand);

        % Eq. 4 (Supplemento): Derivata dello spessore di spostamento (delta*)
        ddeltastar_ds = ddelta_ds - (1 ./ (Ue.^2 + eps)) .* dUe_ds .* integral_term;

        % Integrazione finale per ottenere delta* lungo il profilo
        delta_st_sol = dst_0 + cumtrapz(s, ddeltastar_ds);

    end

    % function [theta_sol, delta_st_sol] = solve_turbulent_side(nu_val, s, Ue, th_0, dst_0)
    % 
    %     % Costanti del modello di regressione spostate QUI DENTRO
    %     Cc = 1.45;        
    %     CRe_inf = 0.0024; 
    %     Cm = 7.23;        
    % 
    %     % Forziamo array a colonna per sicurezza
    %     s = s(:); Ue = Ue(:);
    % 
    %     % 1. Calcolo del gradiente di velocità esterna
    %     dUe_ds = gradient(Ue, s);
    % 
    %     % 2. Funzioni interpolanti per i sub-step numerici di ode45.
    %     Ue_fun = @(sq) interp1(s, Ue, sq, 'makima');
    %     dUe_ds_fun = @(sq) interp1(s, dUe_ds, sq, 'makima');
    % 
    %     % 3. Definizione dell'Equazione Differenziale per theta.
    %     ode_theta = @(sq, th) (CRe_inf / 2) ...
    %                         - (Cm / 2) * (th / (Ue_fun(sq) + eps)) * dUe_ds_fun(sq) ...
    %                         + (Cc / 2) * (nu_val / (Ue_fun(sq) * th + eps));
    % 
    %     % 4. Risoluzione ODE
    %     [~, theta_sol] = ode45(ode_theta, s, th_0);
    %     theta_sol = theta_sol(:); % Sicurezza colonna
    % 
    %     Ue_safe = max(Ue, 1e-3);
    %     dUe_ds_lim = max(min(dUe_ds, 50), -50); 
    % 
    %     % Calcolo del Fattore di Forma Turbolento ideale (H_eq di Buri)
    %     Re_theta = (Ue_safe .* theta_sol) ./ nu_val;
    %     Gamma = - (theta_sol ./ Ue_safe) .* dUe_ds_lim .* (Re_theta.^0.25);
    %     H_eq = 1.4 + 0.9 * max(Gamma, 0);
    %     H_eq = min(H_eq, 2.8); 
    % 
    %     delta_st_sol = zeros(size(s));
    %     delta_st_sol(1) = dst_0;
    % 
    %     % ZPG
    %     % Calcoliamo il rateo di crescita minimo. Diviso per 8 per passare da delta a delta*
    %     % (Offset di 0.05 per evitare singolarità a s=0)
    %     Re_s = (Ue_safe .* (s - s(1) + 0.05)) ./ nu_val; 
    %     ddeltastar_zpg_ds = (0.296 / 8) .* (Re_s .^ -0.2); 
    % 
    %     for i = 2:length(s)
    %         ds = s(i) - s(i-1);
    % 
    %         % Il valore ideale a regime che la fisica "cruda" vorrebbe raggiungere
    %         dst_target = theta_sol(i) * H_eq(i);
    % 
    %         % La crescita minima fisica garantita dal termine ZPG
    %         crescita_minima = ddeltastar_zpg_ds(i) * ds;
    % 
    %         % --- COSTRUZIONE DELLA PASSERELLA ANTI-SALTO ---
    %         if dst_target <= delta_st_sol(i-1)
    %             % FASE 1 (La Buca): H è crollato, dst_target vorrebbe saltare in basso.
    %             % Invece di cadere, avanziamo dritti usando il dolce rateo ZPG!
    %             delta_st_sol(i) = delta_st_sol(i-1) + crescita_minima;
    %         else
    %             % FASE 2 (La Risalita): L'ala ha ingrossato il target oltre la passerella.
    %             % Ci agganciamo alla curva ideale usando un filtro passa-basso 
    %             % per mantenere la curva "burrosa" senza creare angoli (kink).
    %             relax = 0.15; % 15% target, 85% inerzia precedente
    %             delta_st_sol(i) = relax * dst_target + (1 - relax) * (delta_st_sol(i-1) + crescita_minima);
    %         end
    %     end
    % end

end