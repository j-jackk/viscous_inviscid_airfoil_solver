function [Qt, Cp, mids, n_vec, panels, t_vec] = panels_method_v3(N, alpha, Q_inf,filename, profilo)
    
    Q_inf_a = [Q_inf * cos(alpha), Q_inf * sin(alpha)];

    if all(isnan(profilo))
        coordinate = import_airfoil(filename);
        [~, idx] = min(coordinate(:,1));

        x_u = flip(coordinate(1:idx, 1));
        z_u = flip(coordinate(1:idx, 2));
        x_l = coordinate(idx:end, 1);
        z_l = coordinate(idx:end, 2);

        [~, id_u] = unique([x_u, z_u], 'rows', 'stable'); x_u = x_u(id_u); z_u = z_u(id_u);
        [~, id_l] = unique([x_l, z_l], 'rows', 'stable'); x_l = x_l(id_l); z_l = z_l(id_l);

        s_u = [0; cumsum(sqrt(diff(x_u).^2 + diff(z_u).^2))];
        s_l = [0; cumsum(sqrt(diff(x_l).^2 + diff(z_l).^2))];

        c = max(coordinate(:,1));
        beta = linspace(0, pi, N + 1)';
        distribuzione = 0.5 * (1 - cos(beta)); 
        s_target_u = s_u(end) * distribuzione;
        s_target_l = s_l(end) * distribuzione;

        corner_x_u = pchip(s_u, x_u, s_target_u);
        corner_z_u = pchip(s_u, z_u, s_target_u);

        corner_x_l = pchip(s_l, x_l, s_target_l);
        corner_z_l = pchip(s_l, z_l, s_target_l);

        panels_u = [corner_x_u, corner_z_u];
        panels_l = [corner_x_l, corner_z_l];

        panels = [panels_l(N+1:-1:2, :) ; panels_u];
        panels(1,:) = [c, panels(end,2)]; 
        panels(end,1) = c;
    else
        panels = profilo;
        c = max(profilo(:,1)) - min(profilo(:,1));
    end
   
    N_tot = 2 * N;
    mids = zeros(N_tot, 2);
    n_vec = zeros(N_tot, 2);
    t_vec = zeros(N_tot, 2);
    sigma = zeros(N_tot, 1);
    coll_p = zeros(N_tot, 2);
    gamma = zeros(N_tot, 1);

    for k = 1: N_tot
        
        mids(k,1) = ((panels(k,1) + panels(k+1,1))/2); 
        mids(k,2) = ((panels(k,2) + panels(k+1,2))/2);
        vec = panels(k+1,:) - panels(k,:);
        gamma(k) = atan2(vec(2),vec(1));
        n_vec(k, :) = [sin(gamma(k)), -cos(gamma(k))];
        t_vec(k, :) = [cos(gamma(k)), sin(gamma(k))];
        sigma(k) = dot(n_vec(k,:), Q_inf_a);
        coll_p(k,:) =(n_vec(k,:) .* 0.05 * norm(vec)) + mids(k,:);
        
    end
    
    coll_p = mids;


    % figure;
    % plot(panels(:,1), panels(:,2), 'b', 'LineWidth', 2);
    % hold on
    % plot(panels(:,1), panels(:,2), 'rO', 'LineWidth', 4);
    % axis equal
    % hold on
    % plot(mids(:,1), mids(:,2), '^', 'Color', 'g', 'LineWidth', 1);
    % hold on
    % quiver(mids(:,1), mids(:,2), n_vec(:,1), n_vec(:,2), 0.05, 'k');
    % % hold on
    % % quiver(mids(:,1), mids(:,2), t_vec(:,1), t_vec(:,2), 0.05, 'm');
    % legend('Pannelli', 'Spigoli', 'Punti Medi', 'Normali');
    % axis equal
    % grid on

    % % =========================================================================
    % % --- SEZIONE PLOT GEOMETRIA PANNELLI (SFONDO BIANCO E EXPORT 600 DPI) ---
    % % =========================================================================
    % % Colori personalizzati ad alto contrasto
    % col_pan = [0 0.35 0.85];   % Blu elettrico per i pannelli
    % col_edge = [0.9 0 0];      % Rosso acceso per gli spigoli
    % col_mid = [0 0.7 0.2];     % Verde brillante per i punti medi
    % col_norm = [0.2 0.2 0.2];  % Grigio molto scuro per le normali
    % 
    % % --- FIGURA 1: PANORAMICA COMPLETA ---
    % fig_full = figure('Name', 'Discretizzazione Completa', 'Position', [100, 100, 1000, 600], 'Color', 'w');
    % hold on; axis equal; box on;
    % 
    % % FORZATURA SFONDO BIANCO DEGLI ASSI
    % set(gca, 'Color', 'w'); 
    % 
    % % 1. Pannelli (Linee continue)
    % plot(panels(:,1), panels(:,2), '-', 'Color', col_pan, 'LineWidth', 2.0, 'DisplayName', 'Pannelli');
    % 
    % % 2. Spigoli (Nodi rossi riempiti con bordo nero)
    % plot(panels(:,1), panels(:,2), 'o', 'Color', 'k', 'MarkerFaceColor', col_edge, ...
    %      'MarkerSize', 7, 'LineWidth', 1, 'DisplayName', 'Spigoli');
    % 
    % % 3. Punti Medi (Triangoli verdi riempiti con bordo nero)
    % plot(mids(:,1), mids(:,2), '^', 'Color', 'k', 'MarkerFaceColor', col_mid, ...
    %      'MarkerSize', 7, 'LineWidth', 1, 'DisplayName', 'Punti Medi');
    % 
    % % 4. Versori Normali (Lunghezza regolata dal parametro 0.04)
    % quiver(mids(:,1), mids(:,2), n_vec(:,1), n_vec(:,2), 0.04, ...
    %        'Color', col_norm, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'DisplayName', 'Normali');
    % 
    % % Estetica Assi
    % set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');
    % xlabel('x/c', 'FontWeight', 'bold', 'FontSize', 16);
    % ylabel('z/c', 'FontWeight', 'bold', 'FontSize', 16);
    % 
    % % Legenda: Sfondo bianco, bordo nero e TESTO NERO forzato
    % lgd1 = legend('Location', 'northeast', 'FontSize', 14);
    % lgd1.ItemTokenSize = [25, 18]; 
    % set(lgd1, 'Color', 'w', 'EdgeColor', 'k', 'TextColor', 'k'); 
    % 
    % 
    % % --- FIGURA 2: DETTAGLIO LEADING EDGE (ZOOM) ---
    % fig_zoom = figure('Name', 'Dettaglio Leading Edge', 'Position', [150, 150, 1000, 600], 'Color', 'w');
    % hold on; axis equal; box on;
    % 
    % % FORZATURA SFONDO BIANCO DEGLI ASSI
    % set(gca, 'Color', 'w'); 
    % 
    % % Ridisegniamo con spessori e dimensioni dei marker maggiorati per lo zoom
    % plot(panels(:,1), panels(:,2), '-', 'Color', col_pan, 'LineWidth', 2.5, 'DisplayName', 'Pannelli');
    % plot(panels(:,1), panels(:,2), 'o', 'Color', 'k', 'MarkerFaceColor', col_edge, ...
    %      'MarkerSize', 9, 'LineWidth', 1, 'DisplayName', 'Spigoli');
    % plot(mids(:,1), mids(:,2), '^', 'Color', 'k', 'MarkerFaceColor', col_mid, ...
    %      'MarkerSize', 9, 'LineWidth', 1, 'DisplayName', 'Punti Medi');
    % 
    % % Quiver senza autoscale (0) moltiplicato manualmente
    % scala_frecce = 0.015; 
    % quiver(mids(:,1), mids(:,2), n_vec(:,1)*scala_frecce, n_vec(:,2)*scala_frecce, 0, ...
    %        'Color', col_norm, 'LineWidth', 2.0, 'MaxHeadSize', 1.0, 'DisplayName', 'Normali');
    % 
    % % Inquadratura forzata sul bordo d'attacco
    % x_le = min(panels(:,1));
    % 
    % % ALLARGHIAMO LE X per creare un rettangolo orizzontale (es. ampiezza X = 0.35)
    % xlim([x_le - 0.02, x_le + 0.33]); 
    % 
    % % Manteniamo le Y uguali (ampiezza Y = 0.20)
    % ylim([-0.10, 0.10]);
    % 
    % % Estetica Assi Zoom
    % set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');
    % xlabel('x/c', 'FontWeight', 'bold', 'FontSize', 16);
    % ylabel('z/c', 'FontWeight', 'bold', 'FontSize', 16);
    % 
    % % Legenda Zoom: Sfondo bianco, bordo nero e TESTO NERO forzato
    % lgd2 = legend('Location', 'northeast', 'FontSize', 14);
    % set(lgd2, 'Color', 'w', 'EdgeColor', 'k', 'TextColor', 'k'); 
    % 
    % % =========================================================================
    % % --- ESPORTAZIONE DELLE IMMAGINI (PNG, 600 DPI) ---
    % % =========================================================================
    % % Assicurati che non vi siano barre degli strumenti visibili durante il salvataggio
    % ax_full = fig_full.CurrentAxes; ax_full.Toolbar.Visible = 'off';
    % ax_zoom = fig_zoom.CurrentAxes; ax_zoom.Toolbar.Visible = 'off';
    % drawnow;
    % 
    % exportgraphics(fig_full, 'discretizzazione_globale.png', 'Resolution', 600, 'BackgroundColor', 'w');
    % exportgraphics(fig_zoom, 'discretizzazione_zoom.png', 'Resolution', 600, 'BackgroundColor', 'w');
    % 
    % fprintf('\n ---> Immagini della discretizzazione salvate a 600 DPI!\n');

    B = zeros(N_tot, N_tot);
    C = zeros(N_tot, N_tot);

    for i = 1:N_tot
        xi = coll_p(i, 1);
        zi = coll_p(i, 2);
        
        for j = 1:N_tot
            xj1 = panels(j, 1);
            zj1 = panels(j, 2);
            xj2 = panels(j+1, 1);
            zj2 = panels(j+1, 2);
            
            B(i, j) = PHICS(1, xi, zi, xj1, zj1, xj2, zj2, gamma(j));
            
            if i == j
                C(i, j) = 0.5;
            else
                C(i, j) = PHICD(1, xi, zi, xj1, zj1, xj2, zj2, gamma(j));
            end
        end
    end

    x_TE = panels(1, 1);
    z_TE = panels(1, 2);
    
    L_wake = 100000 * c; 
    x_wake_end = x_TE + L_wake * cos(alpha);
    z_wake_end = z_TE + L_wake * sin(alpha);
    
    C_w = zeros(N_tot, 1);
    
    for i = 1:N_tot
        xi = mids(i, 1);
        zi = mids(i, 2);
        
        C_w(i) = PHICD(1, xi, zi, x_TE, z_TE, x_wake_end, z_wake_end, alpha);
    end
    
    A = C;
    A(:, 1)   = C(:, 1)   - C_w; 
    A(:, end) = C(:, end) + C_w; 
    
    RHS = -B * sigma;

    mu = A \ RHS;
    
    Qt = zeros(N_tot,1);
    Cp = zeros(N_tot,1);

    for t = 1:N_tot
        if t == N || t == N_tot
            num = mu(t-1) - mu(t);
            dl = norm(coll_p(t-1,:) - coll_p(t, :));
        elseif t == 1 || t == N + 1
            num = mu(t) - mu(t+1);
            dl = norm(coll_p(t, :) - coll_p(t+1,:));
        else
            num = mu(t-1) - mu(t+1);
            dl = norm(coll_p(t-1,:) - coll_p(t+1,:));
        end
       
        Qt_inf = dot(t_vec(t, :), Q_inf_a);
        Qt(t) = (num / dl) + Qt_inf;
        Cp(t) = 1 - ((Qt(t)^2) / (Q_inf^2));
        

    end

    % % --- PLOT STILE XFOIL SDOPPIATO ---
    % figure; hold on; grid on;
    % x_norm = coll_p(:, 1) / c;
    % 
    % % Dividiamo i vettori a metà (N pannelli sotto, N pannelli sopra)
    % x_lower = x_norm(1:N);
    % Cp_lower = Cp(1:N);
    % x_upper = x_norm(N+1:end);
    % Cp_upper = Cp(N+1:end);
    % 
    % % Plottiamo il Ventre in ROSSO e il Dorso in BLU
    % plot(x_lower, Cp_lower, 'rO', 'LineWidth', 4, 'DisplayName', 'Ventre');
    % plot(x_upper, Cp_upper, 'b*', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', 'Dorso', 'LineWidth', 4);
    % 
    % set(gca, 'YDir', 'reverse');
    % xlabel('x/c', 'FontSize', 15); ylabel('C_p', 'FontSize', 15);
    % legend('Location', 'best');
    % title(sprintf('Distribuzione del C_p (N = %d, \\alpha = %.1f deg)', N_tot, alpha * 180/pi), 'FontSize', 17);
    % ylim([-4, 1.5]);



    % % --- PLOT STILE XFOIL SDOPPIATO (BLINDATO PER LE SLIDE) ---
    % 
    % % 1. Creiamo la figura assegnandola alla variabile "fig_cp"
    % fig_cp = figure('Name', 'Distribuzione Cp', 'Color', 'w', 'Position', [100, 100, 1200, 700]); 
    % hold on; grid on; box on;
    % 
    % x_norm = coll_p(:, 1) / c;
    % x_lower = x_norm(1:N);
    % Cp_lower = Cp(1:N);
    % x_upper = x_norm(N+1:end);
    % Cp_upper = Cp(N+1:end);
    % 
    % % 2. Plot dei dati (Marcatori giganti)
    % plot(x_lower, Cp_lower, 'ro', 'MarkerSize', 12, 'LineWidth', 2.5, 'DisplayName', 'Ventre');
    % plot(x_upper, Cp_upper, 'b*', 'MarkerSize', 14, 'LineWidth', 2.5, 'DisplayName', 'Dorso');
    % 
    % % 3. FORZATURA ASSI E SFONDO (ANTI-DARK MODE)
    % set(gca, 'YDir', 'reverse', ...
    %          'FontSize', 22, ...
    %          'LineWidth', 2, ...
    %          'Color', 'w', ...          % Forza lo sfondo del grafico a bianco
    %          'XColor', 'k', ...         % Forza l'asse X e i numeri a nero
    %          'YColor', 'k', ...         % Forza l'asse Y e i numeri a nero
    %          'GridColor', 'k', ...      % Forza le linee della griglia a nero
    %          'GridAlpha', 0.2);         % Regola la visibilità della griglia
    % 
    % % 4. FORZATURA TESTI (Tutti esplicitamente neri)
    % xlabel('x/c', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k'); 
    % ylabel('Coefficiente di Pressione C_p', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k');
    % title(sprintf('Distribuzione del C_p (N = %d, \\alpha = %.1f^\\circ)', N_tot, alpha * 180/pi), 'FontSize', 28, 'FontWeight', 'bold', 'Color', 'k');
    % 
    % % 5. FORZATURA LEGENDA (Sfondo bianco, testo e bordi neri)
    % lgd = legend('Location', 'best', 'FontSize', 24, 'TextColor', 'k');
    % set(lgd, 'Color', 'w', 'EdgeColor', 'k', 'LineWidth', 2);
    % 
    % ylim([-4, 1.5]);
    % 
    % % =========================================================================
    % % --- ESPORTAZIONE DELLE IMMAGINI (PNG, 600 DPI) ---
    % % =========================================================================
    % % Nascondiamo la toolbar per evitare "sporco" visivo nell'immagine salvata
    % ax_cp = gca; 
    % ax_cp.Toolbar.Visible = 'off';
    % drawnow;
    % 
    % % Esportazione usando la variabile corretta della figura (fig_cp)
    % exportgraphics(fig_cp, 'cp_pannelli.png', 'Resolution', 600, 'BackgroundColor', 'w');
    % fprintf('>>> Grafico Cp salvato con successo a 600 DPI! <<<\n');
    % 
    % % Ripristiniamo la toolbar per permetterti di usare la figura a schermo
    % ax_cp.Toolbar.Visible = 'on';



    % % =========================================================================
    % % --- 1. LETTURA DATI XFOIL ---
    % % =========================================================================
    % % Inserisci qui il nome esatto del file di testo in cui hai salvato i dati
    % nome_file_xfoil = 'cp_0012_5.txt'; 
    % 
    % if isfile(nome_file_xfoil)
    %     % Leggiamo il file saltando le prime 3 righe di testo (Header)
    %     tabella_xfoil = readmatrix(nome_file_xfoil, 'NumHeaderLines', 3);
    %     x_xfoil  = tabella_xfoil(:, 1); % Estrae la colonna 1 (x)
    %     cp_xfoil = tabella_xfoil(:, 3); % Estrae la colonna 3 (Cp)
    %     dati_xfoil_caricati = true;
    % else
    %     warning('File %s non trovato! Il grafico verrà generato senza XFOIL.', nome_file_xfoil);
    %     dati_xfoil_caricati = false;
    % end
    % 
    % % =========================================================================
    % % --- 2. PLOT STILE XFOIL E CONFRONTO (OTTIMIZZATO PER CHIAREZZA) ---
    % % =========================================================================
    % fig_cp = figure('Name', 'Distribuzione Cp con XFOIL', 'Color', 'w', 'Position', [100, 100, 1200, 700]); 
    % hold on; grid on; box on;
    % 
    % % Dati del TUO solutore
    % x_norm = coll_p(:, 1) / c;
    % x_lower = x_norm(1:N);
    % Cp_lower = Cp(1:N);
    % x_upper = x_norm(N+1:end);
    % Cp_upper = Cp(N+1:end);
    % 
    % % 1. PLOT DI XFOIL (Sfondo/Riferimento)
    % % Lo disegniamo PER PRIMO in modo che finisca "sotto" i tuoi punti.
    % % Usiamo una linea continua nera, spessa ed elegante.
    % if dati_xfoil_caricati
    %     plot(x_xfoil, cp_xfoil, 'k-', 'LineWidth', 3, 'DisplayName', 'XFOIL');
    % end
    % 
    % % 2. PLOT DEL TUO SOLUTORE (In primo piano)
    % % Ventre: Cerchi rossi con contorno nero (MarkerFaceColor riempie il cerchio)
    % plot(x_lower, Cp_lower, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5, 'DisplayName', 'Ventre');
    % 
    % % Dorso: Triangoli blu con contorno nero
    % plot(x_upper, Cp_upper, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Dorso');
    % 
    % % FORZATURA ASSI E SFONDO (ANTI-DARK MODE)
    % set(gca, 'YDir', 'reverse', ...
    %          'FontSize', 22, ...
    %          'LineWidth', 2, ...
    %          'Color', 'w', ...          
    %          'XColor', 'k', ...         
    %          'YColor', 'k', ...         
    %          'GridColor', 'k', ...      
    %          'GridAlpha', 0.3);         
    % 
    % % FORZATURA TESTI 
    % xlabel('x/c', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k'); 
    % ylabel('Coefficiente di Pressione C_p', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k');
    % title(sprintf('Distribuzione del C_p (N = %d, \\alpha = %.1f^\\circ)', N_tot, alpha * 180/pi), 'FontSize', 28, 'FontWeight', 'bold', 'Color', 'k');
    % 
    % % FORZATURA LEGENDA 
    % lgd = legend('Location', 'northeast', 'FontSize', 24, 'TextColor', 'k');
    % set(lgd, 'Color', 'w', 'EdgeColor', 'k', 'LineWidth', 2);
    % 
    % ylim([-4, 1.5]);
    % % =========================================================================
    % % --- ESPORTAZIONE DELLE IMMAGINI (PNG, 600 DPI) ---
    % % =========================================================================
    % ax_cp = gca; 
    % ax_cp.Toolbar.Visible = 'off';
    % drawnow;
    % 
    % exportgraphics(fig_cp, 'cp_confronto_xfoil.png', 'Resolution', 600, 'BackgroundColor', 'w');
    % fprintf('>>> Grafico Cp con XFOIL salvato con successo a 600 DPI! <<<\n');
    % 
    % ax_cp.Toolbar.Visible = 'on';
end


function phi_s = PHICS(sigma, x_col, z_col, x1, z1, x2, z2, gamma)
    
    DCM = [[cos(gamma), sin(gamma)]; [-sin(gamma), cos(gamma)]];
    
    L_0 = [x2 - x1; z2 - z1];
    x_loc = x_col - x1;
    z_loc = z_col - z1;
    x_z = DCM * [x_loc; z_loc];
    L_vec = DCM * L_0;
    L = norm(L_vec);

    r1_sq = x_z(1)^2 + x_z(2)^2;
    r2_sq = (x_z(1) - L)^2 + x_z(2)^2;
    
    th1 = atan2(x_z(2), x_z(1));
    th2 = atan2(x_z(2), x_z(1) - L);
   
    term1 = x_z(1) * log(r1_sq);
    term2 = (x_z(1) - L) * log(r2_sq);
    term3 = 2 * x_z(2) * (th2 - th1);
    
    phi_s = (sigma / (4 * pi)) * (term1 - term2 + term3);


end

function phi_d = PHICD(mu, x_col, z_col, x1, z1, x2, z2, gamma)
    
    DCM = [[cos(gamma), sin(gamma)]; [-sin(gamma), cos(gamma)]];
    
    L_0 = [x2 - x1; z2 - z1];
    x_loc = x_col - x1;
    z_loc = z_col - z1;
    x_z = DCM * [x_loc; z_loc];
    L_vec = DCM * L_0;
    L = norm(L_vec);
    
    th1 = atan2(x_z(2), x_z(1));
    th2 = atan2(x_z(2), x_z(1) - L);
    
    phi_d = (-mu / (2 * pi)) * (th2 - th1);


end