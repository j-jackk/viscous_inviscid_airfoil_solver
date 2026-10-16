function coords = import_airfoil_v1(filename)
% IMPORT_AIRFOIL Legge, pulisce e formatta le coordinate di un profilo alare.
% Restituisce una matrice [N x 2] che parte dal Trailing Edge del Dorso,
% passa per il Leading Edge, e finisce al Trailing Edge del Ventre.

    % 1. Lettura Intelligente:
    rawData = importdata(filename);
    
    if isstruct(rawData)
        data = rawData.data; 
    else
        data = rawData;
    end
    
    % PULIZIA DATI "SPORCHI":
    data(any(isnan(data), 2), :) = []; % Cancella righe vuote
    data(data(:,1) > 2, :) = [];       % Cancella righe di intestazione numeriche
    
    % Forziamo l'estrazione in vettori colonna per sicurezza
    x = data(:, 1); x = x(:);
    z = data(:, 2); z = z(:);
    
    % 2. Rilevamento Automatico del Formato (Lednicer vs Selig)
    salti = diff(x);
    idx_salto = find(salti < -0.5, 1); 
    
    if ~isempty(idx_salto)
        % --- FORMATO LEDNICER (Curve separate) ---
        x_u = x(1 : idx_salto);
        z_u = z(1 : idx_salto);
        
        x_l = x(idx_salto+1 : end);
        z_l = z(idx_salto+1 : end);
    else
        % --- FORMATO SELIG (Giro continuo) ---
        [~, idx_LE] = min(x); 
        
        x_u = x(1 : idx_LE);
        z_u = z(1 : idx_LE);
        
        x_l = x(idx_LE : end);
        z_l = z(idx_LE : end);
        
        % Il formato Selig parte dalla coda e va al naso. 
        % Lo ribaltiamo così va dal naso alla coda.
        x_u = flip(x_u);
        z_u = flip(z_u);
    end
    
    % 3. Controllo Anti-Errore (Dorso vs Ventre)
    if mean(z_u) < mean(z_l)
        temp_x = x_u; temp_z = z_u;
        x_u = x_l; z_u = z_l;
        x_l = temp_x; z_l = temp_z;
    end
    
    % --- 4. IL FIX DEFINITIVO: INTERSEZIONE DELLE TANGENTI ---
    m_u = (z_u(end) - z_u(end-1)) / (x_u(end) - x_u(end-1));
    m_l = (z_l(end) - z_l(end-1)) / (x_l(end) - x_l(end-1));

    denom = m_u - m_l;

    % Valvola di sicurezza: se il profilo è già tronco e parallelo
    if abs(denom) < 1e-6
        x_int = max(x_u(end), x_l(end));
        z_int = (z_u(end) + z_l(end)) / 2;
    else
        x_int = (m_u * x_u(end) - m_l * x_l(end) + z_l(end) - z_u(end)) / denom;
        z_int = m_u * (x_int - x_u(end)) + z_u(end);
    end

    % Sostituiamo gli ultimi punti con il nuovo punto a punta
    x_u(end) = x_int; z_u(end) = z_int;
    x_l(end) = x_int; z_l(end) = z_int;
    
    % --- 5. ASSEMBLAGGIO DEL VETTORE FINALE ---
    % Al termine della Sezione 4, SIA x_u CHE x_l vanno dal Naso alla Coda.
    
    % Il Dorso deve andare dalla Coda al Naso, quindi lo ribaltiamo (flip)
    x_dorso_path = flip(x_u);
    z_dorso_path = flip(z_u);
    
    % Il Ventre va già dal Naso alla Coda.
    % Usiamo (2:end) per NON duplicare il punto esatto del Naso [0,0] 
    % che è già stato inserito come ultimo punto del dorso!
    x_ventre_path = x_l(2:end);
    z_ventre_path = z_l(2:end);
    
    % Creiamo la matrice finale [X, Z] unendo i due percorsi
    coords = [x_dorso_path, z_dorso_path; 
              x_ventre_path, z_ventre_path];

    % % --- VERIFICA VISIVA ---
    % figure;
    % hold on; grid on; axis equal;
    % plot(x_u, z_u, 'y.-', 'DisplayName', 'Dorso (Upper)');
    % plot(x_l, z_l, 'r.-', 'DisplayName', 'Ventre (Lower)');
    % plot(x_int, z_int, 'go', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Intersezione TE');
    % legend('Location', 'best');
    % title('Profilo Estrapolato (Chiusura Morbida)'); 

end