clear all

% Simulation configuration
study_case = 'AIRFOIL';     % choose 'AIRFOIL' or 'BLASIUS'

reynolds_number = 3 * 10^6;
rho = 1.225;
dynamic_viscosity = 1.78 * 10^-5;
kinematic_viscosity = dynamic_viscosity / rho;
chord = 1.0;
U_inf = reynolds_number * dynamic_viscosity / (rho * chord);
alpha = deg2rad(0);
num_panels = 160;
coordinate_file = 'n0010_coords.dat';
load('Tab_Thwaites.mat');

% Select the requested simulation case
switch study_case
    case 'AIRFOIL'
        max_iterations = 40;
        displacement_tolerance = 1e-4;
        relaxation_factor = 0.2;
        current_file = coordinate_file;
        iteration = 0;
        convergence_error = 1;
        airfoil_profile = nan((2*num_panels)+1, 2);
        previous_displacement_thickness_full = zeros(2*num_panels, 1);

        % Iterate until the displacement thickness converges
        while convergence_error > displacement_tolerance && iteration < max_iterations
            iteration = iteration + 1;

            % Potential-flow panel solver
            [tangential_velocity, pressure_coefficient, panel_midpoints, panel_normals, ...
                previous_profile, panel_tangents] = ...
                panels_method_v3(num_panels, alpha, U_inf, current_file, airfoil_profile);

            if iteration == 1
                initial_chord = max(previous_profile(:,1));
                profile_coordinates = previous_profile .* (chord / initial_chord);
                mesh_coordinates = panel_midpoints .* (chord / initial_chord);
            end

            max_chordwise_x = max(profile_coordinates(:,1));
            flow_aligned_x = profile_coordinates(:,1) * cos(alpha) + profile_coordinates(:,2) * sin(alpha);
            [~, leading_edge_geometry_idx] = min(flow_aligned_x);
            total_points = length(profile_coordinates(:,1));
            search_window_half_width = floor(total_points * 0.2);
            search_indices = mod( ...
                (leading_edge_geometry_idx - search_window_half_width : ...
                leading_edge_geometry_idx + search_window_half_width) - 1, ...
                total_points) + 1;
            [~, local_max_idx] = max(pressure_coefficient(search_indices));
            split_idx = search_indices(local_max_idx);

            x_lower  = flip(panel_midpoints(1:split_idx, 1));
            y_lower  = flip(panel_midpoints(1:split_idx, 2));
            velocity_lower = flip(abs(tangential_velocity(1:split_idx)));

            x_upper  = panel_midpoints(split_idx+1:end, 1);
            y_upper  = panel_midpoints(split_idx+1:end, 2);
            velocity_upper = abs(tangential_velocity(split_idx+1:end));

            distance_lower = [0; cumsum(sqrt(diff(x_lower).^2 + diff(y_lower).^2))];
            distance_upper = [0; cumsum(sqrt(diff(x_upper).^2 + diff(y_upper).^2))];

            % Laminar boundary layer using the Thwaites method
            momentum_thickness_sq_lower = ...
                ((0.45 * kinematic_viscosity) ./ (velocity_lower .^ 6 + eps)) .* ...
                cumtrapz(distance_lower, velocity_lower .^ 5);
            momentum_thickness_sq_upper = ...
                ((0.45 * kinematic_viscosity) ./ (velocity_upper .^ 6 + eps)) .* ...
                cumtrapz(distance_upper, velocity_upper .^ 5);

            momentum_thickness_lower = sqrt(momentum_thickness_sq_lower);
            momentum_thickness_upper = sqrt(momentum_thickness_sq_upper);

            thwaites_parameter_lower = (momentum_thickness_sq_lower ./ kinematic_viscosity) .* ...
                gradient(velocity_lower, distance_lower);
            thwaites_parameter_upper = (momentum_thickness_sq_upper ./ kinematic_viscosity) .* ...
                gradient(velocity_upper, distance_upper);

            shape_factor_lower = interp1(lambda_table, H_table, thwaites_parameter_lower, 'linear', 'extrap');
            shape_factor_upper = interp1(lambda_table, H_table, thwaites_parameter_upper, 'linear', 'extrap');

            thwaites_l_parameter_lower = interp1(lambda_table, l_table, thwaites_parameter_lower, 'linear', 'extrap');
            thwaites_l_parameter_upper = interp1(lambda_table, l_table, thwaites_parameter_upper, 'linear', 'extrap');

            displacement_thickness_lower = shape_factor_lower .* momentum_thickness_lower;
            displacement_thickness_upper = shape_factor_upper .* momentum_thickness_upper;

            transition_idx_lower = find(thwaites_parameter_lower <= -0.09, 1);
            transition_idx_upper = find(thwaites_parameter_upper <= -0.09, 1);

            momentum_thickness_lower_transition = momentum_thickness_lower(transition_idx_lower);
            displacement_thickness_lower_transition = displacement_thickness_lower(transition_idx_lower);

            momentum_thickness_upper_transition = momentum_thickness_upper(transition_idx_upper);
            displacement_thickness_upper_transition = displacement_thickness_upper(transition_idx_upper);

            [turbulent_momentum_thickness_lower, turbulent_momentum_thickness_upper, ...
                turbulent_displacement_thickness_lower, turbulent_displacement_thickness_upper] = ...
                extension_thwaites_turbulent_v2(kinematic_viscosity, ...
                distance_lower(transition_idx_lower:end), distance_upper(transition_idx_upper:end), ...
                velocity_lower(transition_idx_lower:end), velocity_upper(transition_idx_upper:end), ...
                momentum_thickness_lower_transition, momentum_thickness_upper_transition, ...
                displacement_thickness_lower_transition, displacement_thickness_upper_transition);

            momentum_thickness_lower(transition_idx_lower:end) = turbulent_momentum_thickness_lower;
            displacement_thickness_lower(transition_idx_lower:end) = turbulent_displacement_thickness_lower;

            momentum_thickness_upper(transition_idx_upper:end) = turbulent_momentum_thickness_upper;
            displacement_thickness_upper(transition_idx_upper:end) = turbulent_displacement_thickness_upper;

            % Detect turbulent separation with the Alber parameter
            alber_parameter_lower = -(momentum_thickness_lower ./ (velocity_lower + eps)) .* ...
                gradient(velocity_lower, distance_lower);
            alber_parameter_upper = -(momentum_thickness_upper ./ (velocity_upper + eps)) .* ...
                gradient(velocity_upper, distance_upper);

            separation_threshold = 0.004;

            turbulent_separation_idx_lower = find(alber_parameter_lower(transition_idx_lower:end) >= separation_threshold, 1);
            if ~isempty(turbulent_separation_idx_lower)
                turbulent_separation_idx_lower = turbulent_separation_idx_lower + transition_idx_lower - 1;

                displacement_thickness_lower(turbulent_separation_idx_lower:end) = ...
                    displacement_thickness_lower(turbulent_separation_idx_lower);
            end

            turbulent_separation_idx_upper = find(alber_parameter_upper(transition_idx_upper:end) >= separation_threshold, 1);
            if ~isempty(turbulent_separation_idx_upper)
                turbulent_separation_idx_upper = turbulent_separation_idx_upper + transition_idx_upper - 1;

                displacement_thickness_upper(turbulent_separation_idx_upper:end) = ...
                    displacement_thickness_upper(turbulent_separation_idx_upper);
            end

            displacement_thickness_lower_final = displacement_thickness_lower;
            displacement_thickness_upper_final = displacement_thickness_upper;
            momentum_thickness_lower_final = momentum_thickness_lower;
            momentum_thickness_upper_final = momentum_thickness_upper;

            % Freeze the displacement thickness near the trailing edge
            geometry_freeze_x = 0.95 * max_chordwise_x;

            [~, lower_95_percent_idx] = min(abs(x_lower - geometry_freeze_x));
            displacement_thickness_lower(lower_95_percent_idx:end) = displacement_thickness_lower(lower_95_percent_idx);

            [~, upper_95_percent_idx] = min(abs(x_upper - geometry_freeze_x));
            displacement_thickness_upper(upper_95_percent_idx:end) = displacement_thickness_upper(upper_95_percent_idx);

            displacement_thickness_full = zeros(2*num_panels, 1);
            displacement_thickness_full(1:split_idx) = flip(displacement_thickness_lower);
            displacement_thickness_full(split_idx+1:end) = displacement_thickness_upper;

            % Bound and smooth the displacement-thickness distribution
            max_displacement_thickness = 0.03 * chord;
            displacement_thickness_full(~isfinite(displacement_thickness_full)) = max_displacement_thickness;
            displacement_thickness_full = min(max(displacement_thickness_full, 0), max_displacement_thickness);
            displacement_thickness_full = smoothdata(displacement_thickness_full, 'gaussian', 9);

            displacement_thickness_full = relaxation_factor * displacement_thickness_full + ...
                (1 - relaxation_factor) * previous_displacement_thickness_full;

            % Under-relaxation and convergence check
            convergence_error = max(abs(displacement_thickness_full - previous_displacement_thickness_full));
            fprintf('Iteration %d | Max Thickness variation: %.6e m\n', iteration, convergence_error);
            
            if convergence_error <= displacement_tolerance
                fprintf('>>> Viscous Convergence Reached in %d iterations! <<<\n', iteration);
                break;
            end

            previous_displacement_thickness_full = displacement_thickness_full;

            % Update the airfoil geometry
            corner_displacement_thickness = zeros((2*num_panels)+1,1);
            corner_normals = zeros((2*num_panels)+1,2);

            for i = 2 : (2*num_panels)
                corner_displacement_thickness(i) =  (displacement_thickness_full(i) + displacement_thickness_full(i-1)) / 2;
                normal_sum = panel_normals(i,:) + panel_normals(i-1,:);
                corner_normals(i,:) = normal_sum / norm(normal_sum);
            end

            corner_displacement_thickness(1) = displacement_thickness_full(1);
            corner_displacement_thickness(end) = displacement_thickness_full(end);
            corner_normals(1,:) = panel_normals(1,:);
            corner_normals(end,:) = panel_normals(end,:);

            corner_displacement_thickness = smoothdata(corner_displacement_thickness, 'movmean', 5);
            airfoil_profile = profile_coordinates + (corner_displacement_thickness .* (-corner_normals));

            trailing_edge_lower = airfoil_profile(1, :);
            trailing_edge_upper = airfoil_profile(end, :);

            trailing_edge_midpoint = (trailing_edge_lower + trailing_edge_upper) / 2;

            trailing_edge_shift_lower = trailing_edge_midpoint - trailing_edge_lower;
            trailing_edge_shift_upper  = trailing_edge_midpoint - trailing_edge_upper;

            % Smoothly close the trailing edge
            pinch_start_x = 0.90 * max_chordwise_x;

            for i = 1 : size(airfoil_profile, 1)
                current_x = profile_coordinates(i, 1);

                if current_x > pinch_start_x

                    linear_fraction = (current_x - pinch_start_x) / (max_chordwise_x - pinch_start_x);
                    linear_fraction = min(1, max(0, linear_fraction));

                    f = 0.5 * (1 - cos(pi * linear_fraction));

                    if i <= num_panels+1

                        airfoil_profile(i, :) = airfoil_profile(i, :) + f * trailing_edge_shift_lower;
                    else

                        airfoil_profile(i, :) = airfoil_profile(i, :) + f * trailing_edge_shift_upper;
                    end
                end
            end
        end

    case 'BLASIUS'
        num_points = 200;
        flat_plate_x = (linspace(0, 1, num_points)'.^2) * chord;
        flat_plate_y = zeros(num_points, 1);
        flat_plate_velocity = ones(num_points, 1) * U_inf;
        x_lower = flat_plate_x; y_lower = flat_plate_y; velocity_lower = flat_plate_velocity;
        x_upper = flat_plate_x; y_upper = flat_plate_y; velocity_upper = flat_plate_velocity;
        distance_lower = [0; cumsum(sqrt(diff(x_lower).^2 + diff(y_lower).^2))];
        distance_upper = [0; cumsum(sqrt(diff(x_upper).^2 + diff(y_upper).^2))];
        % Laminar boundary layer using the Thwaites method
        momentum_thickness_sq_lower = ...
                ((0.45 * kinematic_viscosity) ./ (velocity_lower .^ 6 + eps)) .* ...
                cumtrapz(distance_lower, velocity_lower .^ 5);
        momentum_thickness_sq_upper = ...
                ((0.45 * kinematic_viscosity) ./ (velocity_upper .^ 6 + eps)) .* ...
                cumtrapz(distance_upper, velocity_upper .^ 5);
        momentum_thickness_lower = sqrt(momentum_thickness_sq_lower);
        momentum_thickness_upper = sqrt(momentum_thickness_sq_upper);
        thwaites_parameter_lower = (momentum_thickness_sq_lower ./ kinematic_viscosity) .* ...
                gradient(velocity_lower, distance_lower);
        thwaites_parameter_upper = (momentum_thickness_sq_upper ./ kinematic_viscosity) .* ...
                gradient(velocity_upper, distance_upper);

        thwaites_l_parameter_lower = interp1(lambda_table, l_table, thwaites_parameter_lower, 'linear', 'extrap');
        shape_factor_lower = interp1(lambda_table, H_table, thwaites_parameter_lower, 'linear', 'extrap');

        thwaites_l_parameter_upper = interp1(lambda_table, l_table, thwaites_parameter_upper, 'linear', 'extrap');
        shape_factor_upper = interp1(lambda_table, H_table, thwaites_parameter_upper, 'linear', 'extrap');

        wall_shear_stress_lower = dynamic_viscosity .* ...
            (velocity_lower ./ momentum_thickness_lower) .* thwaites_l_parameter_lower;
        wall_shear_stress_upper = dynamic_viscosity .* ...
            (velocity_upper ./ momentum_thickness_upper) .* thwaites_l_parameter_upper;

        displacement_thickness_lower = shape_factor_lower .* momentum_thickness_lower;
        displacement_thickness_upper = shape_factor_upper .* momentum_thickness_upper;

        separation_idx_lower = find(thwaites_parameter_lower <= -0.09, 1);
        separation_idx_upper = find(thwaites_parameter_upper <= -0.09, 1);

end

% Compute aerodynamic coefficients and plot the airfoil boundary layer
if strcmp(study_case, 'AIRFOIL')
    displacement_thickness_lower_final(1:4) = displacement_thickness_lower_final(5);
    displacement_thickness_upper_final(1:4) = displacement_thickness_upper_final(5);

    % Pressure-based lift coefficient
    pressure_lift_coefficient = 0;
    airfoil_chord_length = max(airfoil_profile(:, 1));
    lift_direction = [-sin(alpha), cos(alpha)];

    panel_dx = diff(airfoil_profile(:, 1));
    panel_dz = diff(airfoil_profile(:, 2));
    ds = sqrt(panel_dx.^2 + panel_dz.^2);

    for j = 1 : 2*num_panels
        differential_force = -pressure_coefficient(j) * (ds(j) / airfoil_chord_length) * -panel_normals(j, :);

        pressure_lift_coefficient = pressure_lift_coefficient + dot(differential_force, lift_direction);
    end

    % Squire–Young drag estimate at 95% chord
    x_target = 0.95 * airfoil_chord_length;

    [~, upper_95_percent_idx] = min(abs(x_upper - x_target));
    momentum_thickness_sqy_upper = momentum_thickness_upper_final(upper_95_percent_idx);
    edge_velocity_sqy_upper = abs(velocity_upper(upper_95_percent_idx));
    shape_factor_sqy_upper = displacement_thickness_upper_final(upper_95_percent_idx) / (momentum_thickness_sqy_upper + eps);

    [~, lower_95_percent_idx] = min(abs(x_lower - x_target));
    momentum_thickness_sqy_lower = momentum_thickness_lower_final(lower_95_percent_idx);
    edge_velocity_sqy_lower = abs(velocity_lower(lower_95_percent_idx));
    shape_factor_sqy_lower = displacement_thickness_lower_final(lower_95_percent_idx) / (momentum_thickness_sqy_lower + eps);

    trailing_edge_momentum_thickness_total = momentum_thickness_sqy_upper + momentum_thickness_sqy_lower;
    trailing_edge_mean_velocity = (edge_velocity_sqy_upper + edge_velocity_sqy_lower) / 2;
    trailing_edge_mean_shape_factor = (shape_factor_sqy_upper + shape_factor_sqy_lower) / 2;

    drag_coefficient_sqy = 2 * (trailing_edge_momentum_thickness_total / airfoil_chord_length) * ...
        ((trailing_edge_mean_velocity / U_inf) ^ ((trailing_edge_mean_shape_factor + 5) / 2));

    fprintf('\nLift coefficient (Integral): %.4f\n', pressure_lift_coefficient);
    fprintf('Drag coefficient (Squire-Young): %.5f\n', drag_coefficient_sqy);

    % Remove boundary-layer data downstream of separation
    thickness_color_upper = displacement_thickness_upper_final;
    thickness_color_lower = displacement_thickness_lower_final;

    if exist('separation_idx_upper', 'var') && ~isempty(separation_idx_upper) && separation_idx_upper < transition_idx_upper
        thickness_color_upper(separation_idx_upper:end) = NaN;
    end
    if exist('separation_idx_lower', 'var') && ~isempty(separation_idx_lower) && separation_idx_lower < transition_idx_lower
        thickness_color_lower(separation_idx_lower:end) = NaN;
    end

    if exist('turbulent_separation_idx_upper', 'var') && ~isempty(turbulent_separation_idx_upper)
        thickness_color_upper(turbulent_separation_idx_upper:end) = NaN;
    end
    if exist('turbulent_separation_idx_lower', 'var') && ~isempty(turbulent_separation_idx_lower)
        thickness_color_lower(turbulent_separation_idx_lower:end) = NaN;
    end

    % Assemble the full airfoil surface and its normals
    x_surface = mesh_coordinates(:, 1);
    y_surface = mesh_coordinates(:, 2);
    displacement_thickness_total = [flipud(displacement_thickness_lower_final); displacement_thickness_upper_final(1:end)];
    thickness_color_total = [flipud(thickness_color_lower); thickness_color_upper(1:end)];

    normal_x_total = -panel_normals(:, 1);
    normal_y_total = -panel_normals(:, 2);

    if thickness_color_total(1) <= 1e-6
        thickness_color_total(1) = thickness_color_total(2);
        displacement_thickness_total(1) = displacement_thickness_total(2);
    end
    if thickness_color_total(end) <= 1e-6
        thickness_color_total(end) = thickness_color_total(end - 1);
        displacement_thickness_total(end) = displacement_thickness_total(end - 1);
    end

    max_chordwise_x = max(x_surface);
    x_target = 0.90 * max_chordwise_x;
    [~, leading_edge_idx] = min(x_surface);

    [~, lower_normal_freeze_idx] = min(abs(x_surface(1:leading_edge_idx) - x_target));
    for i = 1 : lower_normal_freeze_idx
        normal_x_total(i) = normal_x_total(lower_normal_freeze_idx);
        normal_y_total(i) = normal_y_total(lower_normal_freeze_idx);
    end

    [~, upper_normal_freeze_relative_idx] = min(abs(x_surface(leading_edge_idx:end) - x_target));
    upper_normal_freeze_idx = leading_edge_idx + upper_normal_freeze_relative_idx - 1;
    for i = upper_normal_freeze_idx : length(x_surface)
        normal_x_total(i) = normal_x_total(upper_normal_freeze_idx);
        normal_y_total(i) = normal_y_total(upper_normal_freeze_idx);
    end

    normal_magnitude = sqrt(normal_x_total.^2 + normal_y_total.^2);
    normal_x_total = normal_x_total ./ normal_magnitude;
    normal_y_total = normal_y_total ./ normal_magnitude;

    % Build the 3D boundary-layer mesh
    spanwise_mesh = linspace(0, 5, 50);
    thickness_plot_scale = 50;
    visual_offset = 0.02;

    x_boundary_layer_mesh = x_surface + ...
        (normal_x_total .* (displacement_thickness_total * thickness_plot_scale + visual_offset));
    y_boundary_layer_mesh = y_surface + ...
        (normal_y_total .* (displacement_thickness_total * thickness_plot_scale + visual_offset));

    X_base = repmat(x_surface, 1, 50);
    Y_base = repmat(y_surface, 1, 50);
    Z_base = repmat(linspace(-0.5, 5.5, 50), length(x_surface), 1);
    gray_color = repmat(reshape([0.4 0.4 0.4], 1, 1, 3), size(X_base, 1), size(X_base, 2));

    X_mesh = repmat(x_boundary_layer_mesh, 1, 50);
    Y_mesh = repmat(y_boundary_layer_mesh, 1, 50);
    Z_mesh = repmat(spanwise_mesh, length(x_boundary_layer_mesh), 1);
    C_mesh = repmat(thickness_color_total, 1, 50);

    % 3D boundary-layer visualization
    figure('Name', '3D boundary-layer visualization', 'Color', 'w');
    hold on;
    grid on;

    set(gca, 'GridColor', [0.3 0.3 0.3], 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

    surf(X_base, Z_base, Y_base, gray_color, ...
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
    ylabel(cb, 'Thickness \delta^* [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

    max_value = max(thickness_color_total(~isnan(thickness_color_total)));
    if isempty(max_value) || max_value <= 0
        max_value = 1e-6;
    end
    clim([0, max_value]);
    axis equal;

    xlabel('Chord X [m]', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 15);
    ylabel('Wingspan Z [m]', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 15);
    zlabel('Thickness Y [m]', 'Color', 'k', 'FontWeight', 'bold');
    title('Boundary Layer - Viscous Simulation', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');

    view(-30, 25);
    view(-30, 25);
    delete(findall(gcf, 'Type', 'light'));
    light('Position', [0, 0, 1], 'Style', 'infinite');
    light('Position', [0, 0, -1], 'Style', 'infinite');
    light('Position', [-1, -0.5, 0], 'Style', 'infinite');
    lighting gouraud;

    % 2D boundary-layer visualization
    boundary_layer_figure = figure('Name','2D boundary-layer visualization',...
                    'Color','w',...
                    'Position',[100 100 1200 500]);

    hold on;

    fill(profile_coordinates(:,1),profile_coordinates(:,2),...
        [0.7 0.7 0.7],...
        'EdgeColor','k',...
        'LineWidth',1.5,...
        'DisplayName','Airfoil');

    scale_2d = 30;

    displacement_thickness_2d = displacement_thickness_total;

    displacement_thickness_2d(isnan(displacement_thickness_2d)) = NaN;

    x_boundary_layer = x_surface + normal_x_total .* (displacement_thickness_2d*scale_2d);
    y_boundary_layer = y_surface + normal_y_total .* (displacement_thickness_2d*scale_2d);

    valid_points_mask = ~(isnan(x_surface) | isnan(y_surface) | ...
            isnan(x_boundary_layer)  | isnan(y_boundary_layer));

    x_profile = x_surface(valid_points_mask);
    y_profile = y_surface(valid_points_mask);

    x_outer = x_boundary_layer(valid_points_mask);
    y_outer = y_boundary_layer(valid_points_mask);

    patch([x_profile; flipud(x_outer)],...
          [y_profile; flipud(y_outer)],...
          [0.35 0.65 1],...
          'FaceAlpha',0.4,...
          'EdgeColor','none',...
          'DisplayName','Displacement Thickness \delta^*');

    plot(x_outer,y_outer,'b--',...
        'LineWidth',2.5,...
        'DisplayName','Outer Boundary Layer Limit');

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

    xlim([min(profile_coordinates(:,1))-x_margin,...
          max(profile_coordinates(:,1))+x_margin]);

    ylim([-0.2 0.2]);

    xlabel('x/c',...
           'FontWeight','bold',...
           'FontSize',20);

    ylabel('z/c',...
           'FontWeight','bold',...
           'FontSize',20);

    title('Visualization of Displacement Thickness \delta^*',...
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
end

% Flat-plate visualization and Blasius validation
if strcmp(study_case, 'BLASIUS')

    thickness_color_upper = displacement_thickness_upper;
    thickness_color_lower = displacement_thickness_lower;

    if ~isempty(separation_idx_upper)
        displacement_thickness_upper(separation_idx_upper:end) = NaN;
        thickness_color_upper(separation_idx_upper:end) = NaN;
    end
    if ~isempty(separation_idx_lower)
        displacement_thickness_lower(separation_idx_lower:end) = NaN;
        thickness_color_lower(separation_idx_lower:end) = NaN;
    end

    dx_upper = gradient(x_upper); dy_upper = gradient(y_upper); ds_upper = sqrt(dx_upper.^2 + dy_upper.^2);
    normal_x_upper = -dy_upper ./ ds_upper; normal_y_upper = dx_upper ./ ds_upper;
    dx_lower = gradient(x_lower); dy_lower = gradient(y_lower); ds_lower = sqrt(dx_lower.^2 + dy_lower.^2);
    normal_x_lower = dy_lower ./ ds_lower; normal_y_lower = -dx_lower ./ ds_lower;

    x_surface = [flipud(x_lower); x_upper(2:end)];
    y_surface = [flipud(y_lower); y_upper(2:end)];
    displacement_thickness_total = [flipud(displacement_thickness_lower); displacement_thickness_upper(2:end)];
    thickness_color_total = [flipud(thickness_color_lower); thickness_color_upper(2:end)];

    dx_surface = gradient(x_surface); dy_surface = gradient(y_surface); ds_surface = sqrt(dx_surface.^2 + dy_surface.^2);
    normal_x_total = -dy_surface ./ ds_surface; normal_y_total = dx_surface ./ ds_surface;

    % Build the 3D boundary-layer mesh
    spanwise_mesh = linspace(0, 5, 50);
    thickness_plot_scale = 80; 
    visual_offset = 0.02;      
    
    x_boundary_layer_mesh = x_surface + ...
        (normal_x_total .* (displacement_thickness_total * thickness_plot_scale + visual_offset));
    y_boundary_layer_mesh = y_surface + ...
        (normal_y_total .* (displacement_thickness_total * thickness_plot_scale + visual_offset));
        
    X_base = repmat(x_surface, 1, 50);
    Y_base = repmat(y_surface, 1, 50);
    Z_base = repmat(linspace(-0.5, 5.5, 50), length(x_surface), 1);
    gray_color = repmat(reshape([0.4 0.4 0.4], 1, 1, 3), size(X_base, 1), size(X_base, 2));
    
    X_mesh = repmat(x_boundary_layer_mesh, 1, 50);
    Y_mesh = repmat(y_boundary_layer_mesh, 1, 50);
    Z_mesh = repmat(spanwise_mesh, length(x_boundary_layer_mesh), 1);
    C_mesh = repmat(thickness_color_total, 1, 50);
    
    % 3D boundary-layer visualization
    figure_3d = figure('Name', '3D boundary-layer mesh', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
    hold on;
    grid on;
    set(gca, 'GridColor', [0.3 0.3 0.3], 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

    % Base airfoil 
    surf(X_base, Z_base, Y_base, gray_color, ...
        'EdgeColor', [0.4 0.4 0.4], ...
        'EdgeAlpha', 0.4, ...
        'FaceAlpha', 1, ...
        'AmbientStrength', 0.5, ...
        'BackFaceLighting', 'reverselit',...
        'FaceLighting', 'gouraud');
        
    % Boundary layer mesh 
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
    ylabel(cb, 'Thickness \delta^* [m]', 'FontSize', 14);
    
    max_value = max(thickness_color_total(~isnan(thickness_color_total)));
    if isempty(max_value) || max_value <= 0
        max_value = 1e-6;
    end
    clim([0, max_value]);
    axis equal;
    
    xlabel('Chord X [m]', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Wingspan Z [m]', 'FontSize', 16, 'FontWeight', 'bold');
    zlabel('Thickness Y [m]', 'FontSize', 16, 'FontWeight', 'bold');
    title('3D Boundary Layer Development', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    
    view(-35, 25);
    
    delete(findall(gcf, 'Type', 'light'));
    light('Position', [0, 0, 1], 'Style', 'infinite');
    light('Position', [0, 0, -1], 'Style', 'infinite');
    light('Position', [-1, -0.5, 0], 'Style', 'infinite');
    lighting gouraud;
    
    % Compare the numerical result with the Blasius solution 
    theta_blasius = 0.664 * sqrt((kinematic_viscosity .* distance_upper) ./ U_inf);
    validation_figure = figure('Name', 'Validation Thwaites vs Blasius', ...
                    'Color', [1 1 1], ...
                    'InvertHardcopy', 'off', ...
                    'Position', [100, 100, 900, 500]);
    p1 = plot(distance_upper, theta_blasius*1000, 'k-', 'LineWidth', 3); hold on;
    p2 = plot(distance_upper, momentum_thickness_upper*1000, 'r--', 'LineWidth', 3);
    grid on; box on;
    
    set(gca, 'FontSize', 16, ...
             'LineWidth', 1.5, ...
             'XColor', 'k', ...
             'YColor', 'k', ...
             'Color', [1 1 1]);
             
    xlabel('Distance along the plate s [m]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('Momentum Thickness \theta [mm]', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    title('Viscous Solver Validation', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');
    
    lgd = legend([p1, p2], {'Analytical Soultion (Blasius)', 'Integral Method (Thwaites)'}, ...
                 'Location', 'northwest', 'FontSize', 14);
    set(lgd, 'Color', [1 1 1], ...
             'EdgeColor', [0 0 0], ...
             'TextColor', [0 0 0]);
end
