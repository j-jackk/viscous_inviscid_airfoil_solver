function [tangential_velocity, Cp, panel_midpoints, normal_vectors, panels, tangent_vectors] = panels_method_v3(N, alpha, freestream_speed, filename, airfoil_profile)

    freestream_velocity = [freestream_speed * cos(alpha), freestream_speed * sin(alpha)];

    if all(isnan(airfoil_profile))
        airfoil_coordinates = import_airfoil(filename);
        [~, leading_edge_idx] = min(airfoil_coordinates(:, 1));

        x_upper = flip(airfoil_coordinates(1:leading_edge_idx, 1));
        z_upper = flip(airfoil_coordinates(1:leading_edge_idx, 2));
        x_lower = airfoil_coordinates(leading_edge_idx:end, 1);
        z_lower = airfoil_coordinates(leading_edge_idx:end, 2);

        [~, unique_idx_upper] = unique([x_upper, z_upper], 'rows', 'stable');
        x_upper = x_upper(unique_idx_upper);
        z_upper = z_upper(unique_idx_upper);
        [~, unique_idx_lower] = unique([x_lower, z_lower], 'rows', 'stable');
        x_lower = x_lower(unique_idx_lower);
        z_lower = z_lower(unique_idx_lower);

        arc_length_upper = [0; cumsum(sqrt(diff(x_upper).^2 + diff(z_upper).^2))];
        arc_length_lower = [0; cumsum(sqrt(diff(x_lower).^2 + diff(z_lower).^2))];

        chord_length = max(airfoil_coordinates(:, 1));
        beta = linspace(0, pi, N + 1)';
        cosine_distribution = 0.5 * (1 - cos(beta));
        target_arc_length_upper = arc_length_upper(end) * cosine_distribution;
        target_arc_length_lower = arc_length_lower(end) * cosine_distribution;

        corner_x_upper = pchip(arc_length_upper, x_upper, target_arc_length_upper);
        corner_z_upper = pchip(arc_length_upper, z_upper, target_arc_length_upper);
        corner_x_lower = pchip(arc_length_lower, x_lower, target_arc_length_lower);
        corner_z_lower = pchip(arc_length_lower, z_lower, target_arc_length_lower);

        panels_upper = [corner_x_upper, corner_z_upper];
        panels_lower = [corner_x_lower, corner_z_lower];

        panels = [panels_lower(N+1:-1:2, :); panels_upper];
        panels(1, :) = [chord_length, panels(end, 2)];
        panels(end, 1) = chord_length;
    else
        panels = airfoil_profile;
        chord_length = max(airfoil_profile(:, 1)) - min(airfoil_profile(:, 1));
    end

    total_panels = 2 * N;
    panel_midpoints = zeros(total_panels, 2);
    normal_vectors = zeros(total_panels, 2);
    tangent_vectors = zeros(total_panels, 2);
    source_strength = zeros(total_panels, 1);
    control_points = zeros(total_panels, 2);
    panel_angles = zeros(total_panels, 1);

    for panel_idx = 1:total_panels
        panel_midpoints(panel_idx, 1) = (panels(panel_idx, 1) + panels(panel_idx + 1, 1)) / 2;
        panel_midpoints(panel_idx, 2) = (panels(panel_idx, 2) + panels(panel_idx + 1, 2)) / 2;
        panel_vector = panels(panel_idx + 1, :) - panels(panel_idx, :);
        panel_angles(panel_idx) = atan2(panel_vector(2), panel_vector(1));
        normal_vectors(panel_idx, :) = [sin(panel_angles(panel_idx)), -cos(panel_angles(panel_idx))];
        tangent_vectors(panel_idx, :) = [cos(panel_angles(panel_idx)), sin(panel_angles(panel_idx))];
        source_strength(panel_idx) = dot(normal_vectors(panel_idx, :), freestream_velocity);
        control_points(panel_idx, :) = (normal_vectors(panel_idx, :) .* 0.05 * norm(panel_vector)) + panel_midpoints(panel_idx, :);
    end

    control_points = panel_midpoints;

    B = zeros(total_panels, total_panels);
    C = zeros(total_panels, total_panels);

    for i = 1:total_panels
        x_control = control_points(i, 1);
        z_control = control_points(i, 2);

        for j = 1:total_panels
            x_panel_start = panels(j, 1);
            z_panel_start = panels(j, 2);
            x_panel_end = panels(j + 1, 1);
            z_panel_end = panels(j + 1, 2);

            B(i, j) = PHICS(1, x_control, z_control, x_panel_start, z_panel_start, x_panel_end, z_panel_end, panel_angles(j));

            if i == j
                C(i, j) = 0.5;
            else
                C(i, j) = PHICD(1, x_control, z_control, x_panel_start, z_panel_start, x_panel_end, z_panel_end, panel_angles(j));
            end
        end
    end

    x_trailing_edge = panels(1, 1);
    z_trailing_edge = panels(1, 2);

    wake_length = 100000 * chord_length;
    x_wake_end = x_trailing_edge + wake_length * cos(alpha);
    z_wake_end = z_trailing_edge + wake_length * sin(alpha);

    wake_influence = zeros(total_panels, 1);

    for i = 1:total_panels
        x_control = panel_midpoints(i, 1);
        z_control = panel_midpoints(i, 2);

        wake_influence(i) = PHICD(1, x_control, z_control, x_trailing_edge, z_trailing_edge, x_wake_end, z_wake_end, alpha);
    end

    A = C;
    A(:, 1) = C(:, 1) - wake_influence;
    A(:, end) = C(:, end) + wake_influence;

    right_hand_side = -B * source_strength;
    doublet_strength = A \ right_hand_side;

    tangential_velocity = zeros(total_panels, 1);
    Cp = zeros(total_panels, 1);

    for panel_idx = 1:total_panels
        if panel_idx == N || panel_idx == total_panels
            velocity_numerator = doublet_strength(panel_idx - 1) - doublet_strength(panel_idx);
            distance_span = norm(control_points(panel_idx - 1, :) - control_points(panel_idx, :));
        elseif panel_idx == 1 || panel_idx == N + 1
            velocity_numerator = doublet_strength(panel_idx) - doublet_strength(panel_idx + 1);
            distance_span = norm(control_points(panel_idx, :) - control_points(panel_idx + 1, :));
        else
            velocity_numerator = doublet_strength(panel_idx - 1) - doublet_strength(panel_idx + 1);
            distance_span = norm(control_points(panel_idx - 1, :) - control_points(panel_idx + 1, :));
        end

        freestream_tangential_velocity = dot(tangent_vectors(panel_idx, :), freestream_velocity);
        tangential_velocity(panel_idx) = (velocity_numerator / distance_span) + freestream_tangential_velocity;
        Cp(panel_idx) = 1 - (tangential_velocity(panel_idx)^2 / freestream_speed^2);
    end
end

function source_potential = PHICS(source_strength, x_control, z_control, x_start, z_start, x_end, z_end, panel_angle)

    rotation_matrix = [[cos(panel_angle), sin(panel_angle)]; [-sin(panel_angle), cos(panel_angle)]];

    panel_vector = [x_end - x_start; z_end - z_start];
    x_offset = x_control - x_start;
    z_offset = z_control - z_start;
    local_coordinates = rotation_matrix * [x_offset; z_offset];
    local_panel_vector = rotation_matrix * panel_vector;
    panel_length = norm(local_panel_vector);

    r_start_squared = local_coordinates(1)^2 + local_coordinates(2)^2;
    r_end_squared = (local_coordinates(1) - panel_length)^2 + local_coordinates(2)^2;

    theta_start = atan2(local_coordinates(2), local_coordinates(1));
    theta_end = atan2(local_coordinates(2), local_coordinates(1) - panel_length);

    term1 = local_coordinates(1) * log(r_start_squared);
    term2 = (local_coordinates(1) - panel_length) * log(r_end_squared);
    term3 = 2 * local_coordinates(2) * (theta_end - theta_start);

    source_potential = (source_strength / (4 * pi)) * (term1 - term2 + term3);
end

function doublet_potential = PHICD(doublet_strength, x_control, z_control, x_start, z_start, x_end, z_end, panel_angle)

    rotation_matrix = [[cos(panel_angle), sin(panel_angle)]; [-sin(panel_angle), cos(panel_angle)]];

    panel_vector = [x_end - x_start; z_end - z_start];
    x_offset = x_control - x_start;
    z_offset = z_control - z_start;
    local_coordinates = rotation_matrix * [x_offset; z_offset];
    local_panel_vector = rotation_matrix * panel_vector;
    panel_length = norm(local_panel_vector);

    theta_start = atan2(local_coordinates(2), local_coordinates(1));
    theta_end = atan2(local_coordinates(2), local_coordinates(1) - panel_length);

    doublet_potential = (-doublet_strength / (2 * pi)) * (theta_end - theta_start);
end
