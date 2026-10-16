clear all
% 0. Caricamento parametri costanti 
caso_studio = 'AIRFOIL';    % Scegli: 'AIRFOIL' oppure 'BLASIUS'

Re = 3 * 10^6;
rho = 1.225;                % [kg / m^2]
visc_d = 1.78 * 10^-5;      % [Pa * s]
visc_c = visc_d / rho;      % [m^2 / s]
corda = 1.0;                  % [m]
U_inf =  Re * visc_d / (rho * corda);   % [m/s]
alpha = deg2rad(0);
N = 160;
file_coordinate = 'n64a410_coordinate.dat';
load('Tabella_Thwaites.mat');

switch caso_studio
    case 'AIRFOIL'
        max_iter = 40;        
        tol_delta = 1e-4;     
        relax = 0.2; 
        file_corrente = file_coordinate; 
        iter = 0;
        err = 1;
        profilo = nan((2*N)+1, 2);
        delta_st_full_old = zeros(2*N, 1);
        
        while err > tol_delta && iter < max_iter
            iter = iter + 1;

            % SOLUTORE POTENZIALE 
            [Qt, Cp, mids, n_vec, profilo_old, t_vec] = panels_method_v3(N, alpha, U_inf, file_corrente, profilo);

            if iter == 1
                c_iniziale = max(profilo_old(:,1));
                coordinate = profilo_old .* (corda / c_iniziale);
                mesh_coords = mids .* (corda / c_iniziale);
                
                figure(14);
                plot(coordinate(:,1), coordinate(:,2));
                axis equal;
                hold on;
            end
            
            c_max = max(coordinate(:,1)); 
            x_vento = coordinate(:,1) * cos(alpha) + coordinate(:,2) * sin(alpha);
            [~, idx_naso_geom] = min(x_vento);
            N_tot = length(coordinate(:,1));
            ampiezza_finestra = floor(N_tot * 0.2);
            indici_ricerca = mod((idx_naso_geom - ampiezza_finestra : idx_naso_geom + ampiezza_finestra) - 1, N_tot) + 1;
            [~, idx_locale] = max(Cp(indici_ricerca));
            idx_ris = indici_ricerca(idx_locale);
            
            x1  = flip(mids(1:idx_ris, 1));
            y1  = flip(mids(1:idx_ris, 2));
            c_l = flip(abs(Qt(1:idx_ris))); 
            
            x2  = mids(idx_ris+1:end, 1);
            y2  = mids(idx_ris+1:end, 2);
            c_u = abs(Qt(idx_ris+1:end));
            
            s1 = [0; cumsum(sqrt(diff(x1).^2 + diff(y1).^2))];      
            s2 = [0; cumsum(sqrt(diff(x2).^2 + diff(y2).^2))];
            
            % =========================================================================
            % 1. SOLUTORE LAMINARE (Thwaites Classico)
            % =========================================================================
            theta_sq_l = (((0.45 * visc_c) ./ (c_l .^ 6 + eps)) .* cumtrapz(s1, c_l .^ 5));
            theta_sq_u = (((0.45 * visc_c) ./ (c_u .^ 6 + eps)) .* cumtrapz(s2, c_u .^ 5));
            
            theta_l = sqrt(theta_sq_l);
            theta_u = sqrt(theta_sq_u);
            
            lambda_l = (theta_sq_l ./ visc_c) .* gradient(c_l, s1);
            lambda_u = (theta_sq_u ./ visc_c) .* gradient(c_u, s2);
            
            H_lambda_l = interp1(lambda_table, H_table, lambda_l, 'linear', 'extrap');
            H_lambda_u = interp1(lambda_table, H_table, lambda_u, 'linear', 'extrap');

            l_lambda_l = interp1(lambda_table, l_table, lambda_l, 'linear', 'extrap');
            l_lambda_u = interp1(lambda_table, l_table, lambda_u, 'linear', 'extrap');
            
            delta_st_l = H_lambda_l .* theta_l;
            delta_st_u = H_lambda_u .* theta_u;

            idx_trans_l = find(lambda_l <= -0.09, 1);
            idx_trans_u = find(lambda_u <= -0.09, 1);
            
            % =========================================================================
            % 3. SOLUTORE TURBOLENTO
            % =========================================================================
            theta_l_0 = theta_l(idx_trans_l);
            delta_st_l_0 = delta_st_l(idx_trans_l);
            
            theta_u_0 = theta_u(idx_trans_u);
            delta_st_u_0 = delta_st_u(idx_trans_u);
            
            [theta_turb_l, theta_turb_u, delta_st_turb_l, delta_st_turb_u] = ...
                extension_thwaites_turbulent_v2(visc_c, ...
                s1(idx_trans_l:end), s2(idx_trans_u:end), ...
                c_l(idx_trans_l:end), c_u(idx_trans_u:end), ...
                theta_l_0, theta_u_0, delta_st_l_0, delta_st_u_0);
                
            theta_l(idx_trans_l:end) = theta_turb_l;
            delta_st_l(idx_trans_l:end) = delta_st_turb_l;
            
            theta_u(idx_trans_u:end) = theta_turb_u;
            delta_st_u(idx_trans_u:end) = delta_st_turb_u;
            % =========================================================================
            % 3.5 CONTROLLO DISTACCO TURBOLENTO 
            % =========================================================================
            % Calcolo parametro di Alber (Lambda) su tutto il profilo
            % Lambda = - (theta / Ue) * (dUe/ds)
            Lambda_l = - (theta_l ./ (c_l + eps)) .* gradient(c_l, s1);
            Lambda_u = - (theta_u ./ (c_u + eps)) .* gradient(c_u, s2);
            
            % Soglia critica di distacco dal paper (0.004)
            soglia_distacco = 0.004; 
            
            % --- Ventre ---
            idx_sep_turb_l = find(Lambda_l(idx_trans_l:end) >= soglia_distacco, 1);
            if ~isempty(idx_sep_turb_l)
                idx_sep_turb_l = idx_sep_turb_l + idx_trans_l - 1; % Indice assoluto
                
                % spessore 'congelato'
                delta_st_l(idx_sep_turb_l:end) = delta_st_l(idx_sep_turb_l); 
            end
            
            % --- Dorso ---
            idx_sep_turb_u = find(Lambda_u(idx_trans_u:end) >= soglia_distacco, 1);
            if ~isempty(idx_sep_turb_u)
                idx_sep_turb_u = idx_sep_turb_u + idx_trans_u - 1;
                
                % spessore 'congelato'
                delta_st_u(idx_sep_turb_u:end) = delta_st_u(idx_sep_turb_u);
            end

            delta_st_l_f = delta_st_l; delta_st_u_f = delta_st_u;
            theta_l_f = theta_l;       theta_u_f = theta_u;
            % =========================================================================
            % 4. CONGELAMENTO GEOMETRIA IN CODA
            % =========================================================================
           
            x_target_geom = 0.95 * c_max;
            
            % Ventre congelato al 95% della corda 
            [~, idx_95_l] = min(abs(x1 - x_target_geom));
            delta_st_l(idx_95_l:end) = delta_st_l(idx_95_l);
            
            % Dorso congelato
            [~, idx_95_u] = min(abs(x2 - x_target_geom));
            delta_st_u(idx_95_u:end) = delta_st_u(idx_95_u);

            delta_st_full = zeros(2*N, 1);
            delta_st_full(1:idx_ris) = flip(delta_st_l); 
            delta_st_full(idx_ris+1:end) = delta_st_u;
            
            % Valvole di sicurezza e Smoothing
            max_delta = 0.03 * corda;                                   % Hp. strato limite <= al 3 % della corda
            delta_st_full(~isfinite(delta_st_full)) = max_delta;        % gestione NaN e inf
            delta_st_full = min(max(delta_st_full, 0), max_delta);      % clamping
            delta_st_full = smoothdata(delta_st_full, 'gaussian', 9);   % gestione gradini

            % --- SOTTO-RILASSAMENTO E CONVERGENZA SULL'ARRAY FISSO ---
            delta_st_full = relax * delta_st_full + (1 - relax) * delta_st_full_old;
            
            err = max(abs(delta_st_full - delta_st_full_old));
            fprintf('Iterazione %d | Max variazione Spessore: %.6e m\n', iter, err);
            
            if err <= tol_delta
                fprintf('>>> Convergenza Viscosa Raggiunta in %d iterazioni! <<<\n', iter);
                break;
            end
            
            delta_st_full_old = delta_st_full;
           
            % AGGIORNAMENTO GEOMETRIA
            delta_corners = zeros((2*N)+1,1);
            n_corners = zeros((2*N)+1,2);
            
            for i = 2 : (2*N) 
                delta_corners(i) =  (delta_st_full(i) + delta_st_full(i-1)) / 2;
                n_sum = n_vec(i,:) + n_vec(i-1,:);
                n_corners(i,:) = n_sum / norm(n_sum);
            end
            
            % Inizializzazione estremi corners
            delta_corners(1) = delta_st_full(1);
            delta_corners(end) = delta_st_full(end);
            n_corners(1,:) = n_vec(1,:);
            n_corners(end,:) = n_vec(end,:);
            
            delta_corners = smoothdata(delta_corners, 'movmean', 5);
            profilo = coordinate + (delta_corners .* (-n_corners));

            TE_ventre = profilo(1, :);
            TE_dorso = profilo(end, :);
            
            % Punto medio d'intersezione 
            TE_medio = (TE_ventre + TE_dorso) / 2;
            
            % Vettori di spostamento necessari per chiudere i due lembi
            shift_ventre = TE_medio - TE_ventre;
            shift_dorso  = TE_medio - TE_dorso;
            
            inizio_pinch = 0.90 * c_max; % Inizio pinch dal 90% della corda
            
            for i = 1 : size(profilo, 1)
                x_curr = coordinate(i, 1);
                
                if x_curr > inizio_pinch
                    % Progresso lineare da 0 (al 90%) a 1 (alla coda)
                    f_lin = (x_curr - inizio_pinch) / (c_max - inizio_pinch);
                    f_lin = min(1, max(0, f_lin));
                    
                    % chiusura 'morbida' con coseno
                    f = 0.5 * (1 - cos(pi * f_lin)); 
                    
                    if i <= N+1 
                        % Lembo del Ventre: verso TE_medio
                        profilo(i, :) = profilo(i, :) + f * shift_ventre;
                    else        
                        % Lembo del Dorso: verso TE_medio
                        profilo(i, :) = profilo(i, :) + f * shift_dorso;
                    end
                end
            end

            figure(14);
            plot(profilo(:,1), profilo(:,2));
            axis equal;
        end
     
    case 'BLASIUS'
        n_punti = 200; 
        x_flat = (linspace(0, 1, n_punti)'.^2) * corda; 
        y_flat = zeros(n_punti, 1);
        u_flat = ones(n_punti, 1) * U_inf;
        x1 = x_flat; y1 = y_flat; c_l = u_flat;
        x2 = x_flat; y2 = y_flat; c_u = u_flat;
        s1 = [0; cumsum(sqrt(diff(x1).^2 + diff(y1).^2))];      
        s2 = [0; cumsum(sqrt(diff(x2).^2 + diff(y2).^2))];
        theta_sq_l = (((0.45 * visc_c) ./ (c_l .^ 6 + eps)) .* cumtrapz(s1, c_l .^ 5));
        theta_sq_u = (((0.45 * visc_c) ./ (c_u .^ 6 + eps)) .* cumtrapz(s2, c_u .^ 5));
        theta_l = sqrt(theta_sq_l);
        theta_u = sqrt(theta_sq_u);
        lambda_l = (theta_sq_l ./ visc_c) .* gradient(c_l, s1);
        lambda_u = (theta_sq_u ./ visc_c) .* gradient(c_u, s2);
        
        l_lambda_l = interp1(lambda_table, l_table, lambda_l, 'linear', 'extrap');
        H_lambda_l = interp1(lambda_table, H_table, lambda_l, 'linear', 'extrap');
        
        l_lambda_u = interp1(lambda_table, l_table, lambda_u, 'linear', 'extrap');
        H_lambda_u = interp1(lambda_table, H_table, lambda_u, 'linear', 'extrap');
        
        tau_w_l = visc_d .* (c_l ./ theta_l) .* l_lambda_l;
        tau_w_u = visc_d .* (c_u ./ theta_u) .* l_lambda_u;
        
        delta_st_l = H_lambda_l .* theta_l;
        delta_st_u = H_lambda_u .* theta_u;
        
        idx_sep_l = find(lambda_l <= -0.09, 1);
        idx_sep_u = find(lambda_u <= -0.09, 1);
        
end

%% Calcolo Forze 
if strcmp(caso_studio, 'AIRFOIL')
    delta_st_l_f(1:4) = delta_st_l_f(5);
    delta_st_u_f(1:4) = delta_st_u_f(5);

   % 1. PORTANZA DI PRESSIONE
    C_l_press = 0;
    c = max(profilo(:,1));
    dir_L = [-sin(alpha), cos(alpha)];
    
    dx_pannelli = diff(profilo(:, 1));
    dz_pannelli = diff(profilo(:, 2));
    ds = sqrt(dx_pannelli.^2 + dz_pannelli.^2);
    
    for j = 1 : 2*N
        dF = -Cp(j) * (ds(j) / c) * -n_vec(j, :);
        
        C_l_press = C_l_press + dot(dF, dir_L);
    end
    
    % 2. SQUIRE-YOUNG (Calcolato al 95% della corda)
    x_target = 0.95 * c;
    
    [~, idx_95_u] = min(abs(x2 - x_target)); 
    theta_sqy_u = theta_u_f(idx_95_u);
    Ue_sqy_u = abs(c_u(idx_95_u));
    H_sqy_u = delta_st_u_f(idx_95_u) / (theta_sqy_u + eps);
    
    [~, idx_95_l] = min(abs(x1 - x_target)); 
    theta_sqy_l = theta_l_f(idx_95_l);
    Ue_sqy_l = abs(c_l(idx_95_l));
    H_sqy_l = delta_st_l_f(idx_95_l) / (theta_sqy_l + eps);
    
    theta_TE_tot = theta_sqy_u + theta_sqy_l;
    Ue_TE_medio = (Ue_sqy_u + Ue_sqy_l) / 2;
    H_TE_medio = (H_sqy_u + H_sqy_l) / 2;
    
    Cd_sqy = 2 * (theta_TE_tot / c) * ((Ue_TE_medio / U_inf) ^ ((H_TE_medio+5)/2));
    
    fprintf('\nCoefficiente di portanza (Cl Integrato): %.4f\n', C_l_press);
    fprintf('Coefficiente di resistenza aerodinamica (Squire-Young): %.5f\n', Cd_sqy); 

    %% Plot
    % --- 8.1 Pulizia dei dati post-distacco ---
    color_st_u = delta_st_u_f;
    color_st_l = delta_st_l_f;
    
    % 1. Cancellazione mesh se si stacca in fase LAMINARE 
    if exist('idx_sep_u', 'var') && ~isempty(idx_sep_u) && idx_sep_u < idx_trans_u
        color_st_u(idx_sep_u:end) = NaN;
    end
    if exist('idx_sep_l', 'var') && ~isempty(idx_sep_l) && idx_sep_l < idx_trans_l
        color_st_l(idx_sep_l:end) = NaN;
    end
    
    % 2. Cancellazione mesh se si stacca in fase TURBOLENTA 
    if exist('idx_sep_turb_u', 'var') && ~isempty(idx_sep_turb_u)
        color_st_u(idx_sep_turb_u:end) = NaN;
    end
    if exist('idx_sep_turb_l', 'var') && ~isempty(idx_sep_turb_l)
        color_st_l(idx_sep_turb_l:end) = NaN;
    end
    
    % --- 8.2 PROFILO BASE E NORMALI UNIFICATE ---
    x_tot = mesh_coords(:, 1);
    y_tot = mesh_coords(:, 2);
    delta_tot = [flipud(delta_st_l_f); delta_st_u_f(1:end)];
    color_tot = [flipud(color_st_l); color_st_u(1:end)];
    
    % Estrazione normali da n_vec 
    nx_tot = -n_vec(:, 1);
    ny_tot = -n_vec(:, 2);
    
    if color_tot(1) <= 1e-6, color_tot(1) = color_tot(2); delta_tot(1) = delta_tot(2); end
    if color_tot(end) <= 1e-6, color_tot(end) = color_tot(end-1); delta_tot(end) = delta_tot(end-1); end
    
    c_max = max(x_tot);
    x_target = 0.90 * c_max;
    [~, idx_naso] = min(x_tot); 
    
    % Congelamento angolo normali Ventre (da TE Lower al 90%)
    [~, idx_90_l] = min(abs(x_tot(1:idx_naso) - x_target));
    for i = 1 : idx_90_l
        nx_tot(i) = nx_tot(idx_90_l);
        ny_tot(i) = ny_tot(idx_90_l);
    end
    
    % Congelamento angolo normali Dorso (dal 90% a TE Upper)
    [~, idx_90_u_rel] = min(abs(x_tot(idx_naso:end) - x_target));
    idx_90_u = idx_naso + idx_90_u_rel - 1;
    for i = idx_90_u : length(x_tot)
        nx_tot(i) = nx_tot(idx_90_u);
        ny_tot(i) = ny_tot(idx_90_u);
    end
    
    % Rinomormalizzazione vettori 
    norm_n = sqrt(nx_tot.^2 + ny_tot.^2);
    nx_tot = nx_tot ./ norm_n;
    ny_tot = ny_tot ./ norm_n;
    
    % --- 8.3 MESH DELLO STRATO LIMITE UNIFICATA ---
    z_mesh = linspace(0, 5, 50);
    scala_spessore_plot = 50; % Fattore di ingrandimento visivo
    offset_visivo = 0.02;    % Spostamento guscio dall'ala per evitare compenetrazioni
     
    x_mesh_unificata = x_tot + (nx_tot .* (delta_tot * scala_spessore_plot + offset_visivo));
    y_mesh_unificata = y_tot + (ny_tot .* (delta_tot * scala_spessore_plot + offset_visivo));
    
    % Matrici 3D airfoil
    X_base = repmat(x_tot, 1, 50); 
    Y_base = repmat(y_tot, 1, 50); 
    Z_base = repmat(linspace(-0.5, 5.5, 50), length(x_tot), 1);
    C_grigio = repmat(reshape([0.4 0.4 0.4], 1, 1, 3), size(X_base, 1), size(X_base, 2));
    
    % Matrici 3D strato limite
    X_mesh = repmat(x_mesh_unificata, 1, 50); 
    Y_mesh = repmat(y_mesh_unificata, 1, 50);
    Z_mesh = repmat(z_mesh, length(x_mesh_unificata), 1);
    C_mesh = repmat(color_tot, 1, 50);
    
    figure('Name', 'Analisi Strato Limite 3D', 'Color', 'w'); 
    hold on; 
    grid on; 
    
    set(gca, 'GridColor', [0.3 0.3 0.3], 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    
    surf(X_base, Z_base, Y_base, C_grigio, ...
        'EdgeColor', [0.4 0.4 0.4], ... 
        'EdgeAlpha', 0.4, ...           
        'FaceAlpha', 1, ...  
        'AmbientStrength', 0.5, ...
        'BackFaceLighting', 'reverselit',...  
        'FaceLighting', 'gouraud');     
    
    surf(X_mesh, Z_mesh, Y_mesh, C_mesh, ...
        'EdgeColor', [0.2 0.2 0.2], ...  
        'EdgeAlpha', 0, ...            
        'FaceColor', 'interp', ...       
        'FaceAlpha', 0.85, ...            
        'FaceLighting', 'gouraud', ...
        'BackFaceLighting', 'reverselit');
                 
    material dull;              
    colormap(jet); 
    cb = colorbar; 
    cb.Color = 'k'; 
    ylabel(cb, 'Spessore \delta^* [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k'); 
    
    max_val = max(color_tot(~isnan(color_tot)));
    if isempty(max_val) || max_val <= 0
        max_val = 1e-6; 
    end
    clim([0, max_val]);
    axis equal; 
    
    % Etichette degli assi 
    xlabel('Corda X [m]', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 15);
    ylabel('Apertura Alare Z [m]', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 15);
    zlabel('Spessore Y [m]', 'Color', 'k', 'FontWeight', 'bold');
    title('Strato Limite - Simulazione Viscosa', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
    
    view(-30, 25); 
    view(-30, 25); 
    delete(findall(gcf, 'Type', 'light')); 
    light('Position', [0, 0, 1], 'Style', 'infinite');
    light('Position', [0, 0, -1], 'Style', 'infinite');
    light('Position', [-1, -0.5, 0], 'Style', 'infinite');
    lighting gouraud;
    
    % =========================================================================
    % --- 8.5 PLOT 2D: PROFILO E STRATO LIMITE (STILE TESI)
    % =========================================================================
    
    fig_bl = figure('Name','Analisi Strato Limite 2D',...
                    'Color','w',...
                    'Position',[100 100 1200 500]);
    
    hold on;
    
    %% =========================================================================
    % 1. PROFILO ALARE
    % =========================================================================
    
    fill(coordinate(:,1),coordinate(:,2),...
        [0.7 0.7 0.7],...
        'EdgeColor','k',...
        'LineWidth',1.5,...
        'DisplayName','Profilo Alare');
    
    %% =========================================================================
    % 2. CALCOLO DELLO STRATO LIMITE
    % =========================================================================
    
    scala_2D = 30;
    
    delta_tot_2d = delta_tot;
    
    % Elimina punti non validi
    delta_tot_2d(isnan(delta_tot_2d)) = NaN;
    
    % Coordinate del bordo esterno dello strato limite
    x_bl = x_tot + nx_tot .* (delta_tot_2d*scala_2D);
    y_bl = y_tot + ny_tot .* (delta_tot_2d*scala_2D);
    
    %% =========================================================================
    % 3. COSTRUZIONE POLIGONO
    % =========================================================================
    
    idx = ~(isnan(x_tot) | isnan(y_tot) | ...
            isnan(x_bl)  | isnan(y_bl));
    
    x_prof = x_tot(idx);
    y_prof = y_tot(idx);
    
    x_ext = x_bl(idx);
    y_ext = y_bl(idx);
    
    % Disegna la fascia dello strato limite
    patch([x_prof; flipud(x_ext)],...
          [y_prof; flipud(y_ext)],...
          [0.35 0.65 1],...
          'FaceAlpha',0.4,...
          'EdgeColor','none',...
          'DisplayName','Spessore di Spostamento \delta^*');
    
    % Bordo esterno
    plot(x_ext,y_ext,'b--',...
        'LineWidth',2.5,...
        'DisplayName','Limite Esterno Strato Limite');
    
    %% =========================================================================
    % 4. ESTETICA
    % =========================================================================
    
    axis equal
    box on
    grid on
    
    set(gca,...
        'GridColor',[0.7 0.7 0.7],...
        'GridAlpha',0.5,...
        'Color','w',...
        'XColor','k',...
        'YColor','k',...
        'LineWidth',1.5,...
        'FontSize',14,...
        'FontWeight','bold');
    
    x_margin = 0.05;
    
    xlim([min(coordinate(:,1))-x_margin,...
          max(coordinate(:,1))+x_margin]);
    
    ylim([-0.2 0.2]);
    
    xlabel('x/c',...
           'FontWeight','bold',...
           'FontSize',20);
    
    ylabel('z/c',...
           'FontWeight','bold',...
           'FontSize',20);
    
    title('Visualizzazione dello Spessore di Spostamento \delta^*',...
          'FontSize',22,...
          'FontWeight','bold',...
          'Color','k');
    
    lgd = legend('Location','southwest',...
                 'FontSize',14);
    
    set(lgd,...
        'Color','w',...
        'EdgeColor','k',...
        'TextColor','k',...
        'LineWidth',1);
    
    %% =========================================================================
    % 6. ESPORTAZIONE
    % =========================================================================

    drawnow;

    exportgraphics(fig_bl,...
        'Strato_Limite_Tesi.png',...
        'Resolution',600,...
        'BackgroundColor','w');

    fprintf('>>> Grafico Strato Limite salvato correttamente a 600 DPI. <<<\n');
end

% if strcmp(caso_studio, 'BLASIUS')
%     % --- 8.1 Pulizia dei dati post-distacco ---
%     color_st_u = delta_st_u;
%     color_st_l = delta_st_l;
% 
%     if ~isempty(idx_sep_u)
%         delta_st_u(idx_sep_u:end) = NaN;
%         color_st_u(idx_sep_u:end) = NaN;
%     end
%     if ~isempty(idx_sep_l)
%         delta_st_l(idx_sep_l:end) = NaN;
%         color_st_l(idx_sep_l:end) = NaN;
%     end
% 
%     % --- 8.2 PROFILO BASE E NORMALI ---
%     % Calcoliamo le normali separate per le frecce della Sezione 8.5
%     dx_u = gradient(x2); dy_u = gradient(y2); ds_u = sqrt(dx_u.^2 + dy_u.^2);
%     nx_u = -dy_u ./ ds_u; tx_u = dx_u ./ ds_u; ny_u = dx_u ./ ds_u; ty_u = dy_u ./ ds_u;  
%     dx_l = gradient(x1); dy_l = gradient(y1); ds_l = sqrt(dx_l.^2 + dy_l.^2);
%     nx_l = dy_l ./ ds_l; tx_l = dx_l ./ ds_l; ny_l = -dx_l ./ ds_l; ty_l = dy_l ./ ds_l;
% 
%     % UNIFICAZIONE PER IL PLOT (Chiude il "buco" sul naso)
%     % Capovolgiamo il ventre e lo attacchiamo al dorso saltando il punto doppio
%     x_tot = [flipud(x1); x2(2:end)]; 
%     y_tot = [flipud(y1); y2(2:end)];
%     delta_tot = [flipud(delta_st_l); delta_st_u(2:end)];
%     color_tot = [flipud(color_st_l); color_st_u(2:end)];
% 
%     % Normali LISCIATE sull'intero profilo continuo per evitare spaccature
%     dx_tot = gradient(x_tot); dy_tot = gradient(y_tot); ds_tot = sqrt(dx_tot.^2 + dy_tot.^2);
%     nx_tot = -dy_tot ./ ds_tot; ny_tot = dx_tot ./ ds_tot; 
% 
%     % --- 8.3 MESH DELLO STRATO LIMITE UNIFICATA ---
%     z_mesh = linspace(0, 5, 50);
%     scala_spessore_plot = 80; % Fattore di ingrandimento visivo
% 
%     % Spingiamo in fuori i punti usando le normali continue
%     x_mesh_unificata = x_tot + (nx_tot .* delta_tot * scala_spessore_plot);
%     y_mesh_unificata = y_tot + (ny_tot .* delta_tot * scala_spessore_plot);
% 
%     % Matrici 3D Ala Base
%     X_base = repmat(x_tot, 1, 50); 
%     Y_base = repmat(y_tot, 1, 50); 
%     Z_base = repmat(linspace(-0.5, 5.5, 50), length(x_tot), 1);
% 
%     % Matrici 3D Strato Limite
%     X_mesh = repmat(x_mesh_unificata, 1, 50); 
%     Y_mesh = repmat(y_mesh_unificata, 1, 50);
%     Z_mesh = repmat(z_mesh, length(x_mesh_unificata), 1);
%     C_mesh = repmat(color_tot, 1, 50);
% 
%     % --- 8.4 Plotting ---
%     figure('Name', 'Analisi Strato Limite 3D Avanzata'); hold on; grid on;
% 
%     % Ingrandisce i numeri sugli assi X, Y, Z
%     set(gca, 'FontSize', 14); 
% 
%     % 1. Disegno ala grigia CONTINUA
%     surf(X_base, Z_base, Y_base, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 1);
% 
%     % 2. Disegno strato limite CONTINUO (texturemap aiuta MATLAB a gestire i NaN colorati)
%     surf(X_mesh, Z_mesh, Y_mesh, C_mesh, 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceColor', 'texturemap');
% 
%     colormap(jet); 
%     cb = colorbar; 
%     ylabel(cb, 'Spessore \delta^* [m]', 'FontSize', 14, 'FontWeight', 'bold'); % Colorbar label più grande
% 
%     % Troviamo il valore massimo da colorare ignorando i NaN
%     max_val = max(color_tot(~isnan(color_tot)));
% 
%     % SCUDO PER IL PLOT: Assicuriamoci che max_val sia STRETTAMENTE maggiore di 0
%     if isempty(max_val) || max_val <= 0
%         max_val = 1e-6; % Usiamo un micron di spessore simbolico per non far crashare clim
%     end
% 
%     % Ora applichiamo la mappa colori in totale sicurezza
%     clim([0, max_val]);
% 
%     axis equal; 
%     % Etichette e titolo ingranditi
%     xlabel('Corda X [m]', 'FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Apertura Alare Z [m]', 'FontSize', 15, 'FontWeight', 'bold');
%     zlabel('Spessore Y [m]', 'FontSize', 15, 'FontWeight', 'bold');
%     title('Strato Limite', 'FontSize', 17, 'FontWeight', 'bold');
%     view(-30, 30); 
%     camlight left; lighting gouraud;
% 
%     % =========================================================================
% 
%     % Calcolo della soluzione analitica esatta di Blasius
%     theta_blasius = 0.664 * sqrt((visc_c .* s2) ./ U_inf);
% 
%     figure('Name', 'Validazione Thwaites vs Blasius');
%     % Ingrandisce i numeri sugli assi per il grafico 2D
%     set(gca, 'FontSize', 14, 'LineWidth', 1.2); 
% 
%     plot(s2, theta_blasius, 'b-', 'LineWidth', 2); hold on; grid on;
%     % Plottiamo il risultato calcolato dal tuo codice (usiamo il dorso, theta_u)
%     plot(s2, theta_u, 'r--', 'LineWidth', 2);
% 
%     % Etichette 2D e titolo più grandi
%     xlabel('Distanza lungo la lastra s [m]', 'FontSize', 16, 'FontWeight', 'bold');
%     ylabel('Spessore Quantità di Moto \theta [mm]', 'FontSize', 16, 'FontWeight', 'bold');
%     title('Validazione Solutore Viscoso: Lastra Piana', 'FontSize', 18, 'FontWeight', 'bold');
% 
%     % Legenda ingrandita
%     lgd = legend('Soluzione Analitica Esatta (Blasius)', 'Metodo Integrale del Codice (Thwaites)', 'Location', 'northwest');
%     lgd.FontSize = 14; 
% 
%     % Calcolo e stampa dell'errore massimo percentuale
%     errore_max_perc = max(abs(theta_u(2:end) - theta_blasius(2:end)) ./ theta_blasius(2:end)) * 100;
%     fprintf('Test BLASIUS Completato. Errore massimo teorico di Thwaites: %.2f%%\n', errore_max_perc);
% end


if strcmp(caso_studio, 'BLASIUS')
    % --- 8.1 Pulizia dei dati post-distacco ---
    color_st_u = delta_st_u;
    color_st_l = delta_st_l;
    
    if ~isempty(idx_sep_u)
        delta_st_u(idx_sep_u:end) = NaN;
        color_st_u(idx_sep_u:end) = NaN;
    end
    if ~isempty(idx_sep_l)
        delta_st_l(idx_sep_l:end) = NaN;
        color_st_l(idx_sep_l:end) = NaN;
    end
    
    % --- 8.2 PROFILO BASE E NORMALI ---
    dx_u = gradient(x2); dy_u = gradient(y2); ds_u = sqrt(dx_u.^2 + dy_u.^2);
    nx_u = -dy_u ./ ds_u; ny_u = dx_u ./ ds_u;  
    dx_l = gradient(x1); dy_l = gradient(y1); ds_l = sqrt(dx_l.^2 + dy_l.^2);
    nx_l = dy_l ./ ds_l; ny_l = -dx_l ./ ds_l;
    
    x_tot = [flipud(x1); x2(2:end)]; 
    y_tot = [flipud(y1); y2(2:end)];
    delta_tot = [flipud(delta_st_l); delta_st_u(2:end)];
    color_tot = [flipud(color_st_l); color_st_u(2:end)];
    
    dx_tot = gradient(x_tot); dy_tot = gradient(y_tot); ds_tot = sqrt(dx_tot.^2 + dy_tot.^2);
    nx_tot = -dy_tot ./ ds_tot; ny_tot = dx_tot ./ ds_tot; 
    
    % --- 8.3 MESH DELLO STRATO LIMITE UNIFICATA ---
    z_mesh = linspace(0, 5, 50);
    scala_spessore_plot = 80;
    
    x_mesh_unificata = x_tot + (nx_tot .* delta_tot * scala_spessore_plot);
    y_mesh_unificata = y_tot + (ny_tot .* delta_tot * scala_spessore_plot);
    
    X_base = repmat(x_tot, 1, 50); Y_base = repmat(y_tot, 1, 50); Z_base = repmat(linspace(-0.5, 5.5, 50), length(x_tot), 1);
    X_mesh = repmat(x_mesh_unificata, 1, 50); Y_mesh = repmat(y_mesh_unificata, 1, 50); Z_mesh = repmat(z_mesh, length(x_mesh_unificata), 1);
    C_mesh = repmat(color_tot, 1, 50);
    
    % --- 8.4 PLOTTING OTTIMIZZATO PER TESI ---
    % Figura 3D
    fig3d = figure('Name', 'Analisi Strato Limite 3D', 'Color', 'w', 'Position', [100, 100, 1000, 600]); 
    surf(X_base, Z_base, Y_base, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8); hold on;
    surf(X_mesh, Z_mesh, Y_mesh, C_mesh, 'EdgeColor', 'none', 'FaceColor', 'texturemap');
    
    colormap(jet); cb = colorbar; ylabel(cb, 'Spessore \delta^* [m]', 'FontSize', 14);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    xlabel('Corda X [m]', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Apertura Alare Z [m]', 'FontSize', 16, 'FontWeight', 'bold');
    zlabel('Spessore Y [m]', 'FontSize', 16, 'FontWeight', 'bold');
    title('Sviluppo Strato Limite 3D', 'FontSize', 18, 'FontWeight', 'bold');
    view(-35, 25); axis equal; camlight left; lighting gouraud;
    
    exportgraphics(fig3d, 'Strato_Limite_3D.png', 'Resolution', 600);
    
    % Validazione Blasius
    theta_blasius = 0.664 * sqrt((visc_c .* s2) ./ U_inf);
    
   % =========================================================================
    % PLOT VALUTAZIONE BLASIUS - VERSIONE DEFINITIVA
    % =========================================================================
    
    % 1. Creazione figura con colore di sfondo FORZATO a bianco
    figVal = figure('Name', 'Validazione Thwaites vs Blasius', ...
                    'Color', [1 1 1], ...         % Sfondo bianco (RGB)
                    'InvertHardcopy', 'off', ...  % Impedisce a MATLAB di invertire i colori al salvataggio
                    'Position', [100, 100, 900, 500]); 
    
    % 2. Plot dei dati
    p1 = plot(s2, theta_blasius*1000, 'k-', 'LineWidth', 3); hold on; 
    p2 = plot(s2, theta_u*1000, 'r--', 'LineWidth', 3);
    
    % 3. Formattazione rigorosa
    grid on; box on; 
    
    % Imposta assi e testi a nero, sfondo area grafico a bianco
    set(gca, 'FontSize', 16, ...
             'LineWidth', 1.5, ...
             'XColor', 'k', ...
             'YColor', 'k', ...
             'Color', [1 1 1]); % Sfondo area grafico bianco
    
    % Etichette e Titolo forzati neri
    xlabel('Distanza lungo la lastra s [m]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('Spessore Quantità di Moto \theta [mm]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    title('Validazione Solutore Viscoso', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');
    
    % 3. Legenda FORZATA a bianco e nero
    lgd = legend([p1, p2], {'Soluzione Analitica (Blasius)', 'Metodo Integrale (Thwaites)'}, ...
                 'Location', 'northwest', 'FontSize', 14);
    
    % Forza esplicitamente lo sfondo della legenda a bianco e i bordi a nero
    set(lgd, 'Color', [1 1 1], ...      % Sfondo bianco (RGB)
             'EdgeColor', [0 0 0], ...  % Bordo nero (RGB)
             'TextColor', [0 0 0]);     % Testo nero (RGB)

    % 4. SALVATAGGIO (DPI elevati e senza barre strumenti)
    % Usiamo 'BackgroundColor' nel salvataggio per ulteriore sicurezza
    exportgraphics(figVal, 'Validazione_Blasius.png', 'Resolution', 600, 'BackgroundColor', 'white');
    
    errore_max_perc = max(abs(theta_u(2:end) - theta_blasius(2:end)) ./ theta_blasius(2:end)) * 100;
    fprintf('Test BLASIUS Completato. Errore massimo teorico di Thwaites: %.2f%%\n', errore_max_perc);
end
