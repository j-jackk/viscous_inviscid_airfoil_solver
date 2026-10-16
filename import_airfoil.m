function coordinates = import_airfoil(filename)

    % Import coordinate data
    raw_data = importdata(filename);

    if isstruct(raw_data)
        data = raw_data.data;
    else
        data = raw_data;
    end

    % Remove invalid rows and possible header values
    data(any(isnan(data), 2), :) = [];
    data(data(:, 1) > 2, :) = [];

    % Ensure column vectors
    x = data(:, 1);
    x = x(:);
    z = data(:, 2);
    z = z(:);

    % Detect Lednicer or Selig coordinate format
    x_differences = diff(x);
    split_idx = find(x_differences < -0.5, 1);

    if ~isempty(split_idx)
        % Lednicer format: upper and lower surfaces are separate
        x_upper = x(1:split_idx);
        z_upper = z(1:split_idx);

        x_lower = x(split_idx + 1:end);
        z_lower = z(split_idx + 1:end);
    else
        % Selig format: the airfoil contour is continuous
        [~, leading_edge_idx] = min(x);

        x_upper = x(1:leading_edge_idx);
        z_upper = z(1:leading_edge_idx);

        x_lower = x(leading_edge_idx:end);
        z_lower = z(leading_edge_idx:end);

        % Reverse upper-surface ordering
        x_upper = flip(x_upper);
        z_upper = flip(z_upper);
    end

    % Ensure the upper and lower surfaces are correctly identified
    if mean(z_upper) < mean(z_lower)
        temp_x = x_upper;
        temp_z = z_upper;
        x_upper = x_lower;
        z_upper = z_lower;
        x_lower = temp_x;
        z_lower = temp_z;
    end

    % Close the trailing edge by extrapolating surface tangents
    slope_upper = (z_upper(end) - z_upper(end - 1)) / ...
                  (x_upper(end) - x_upper(end - 1));
    slope_lower = (z_lower(end) - z_lower(end - 1)) / ...
                  (x_lower(end) - x_lower(end - 1));
    slope_difference = slope_upper - slope_lower;

    % Handle nearly parallel tangents
    if abs(slope_difference) < 1e-6
        x_intersection = max(x_upper(end), x_lower(end));
        z_intersection = (z_upper(end) + z_lower(end)) / 2;
    else
        x_intersection = (slope_upper * x_upper(end) - ...
                          slope_lower * x_lower(end) + ...
                          z_lower(end) - z_upper(end)) / slope_difference;
        z_intersection = slope_upper * (x_intersection - x_upper(end)) + z_upper(end);
    end

    % Set the theoretical trailing-edge intersection
    x_upper(end) = x_intersection;
    z_upper(end) = z_intersection;
    x_lower(end) = x_intersection;
    z_lower(end) = z_intersection;

    % Assemble the final coordinate array
    x_upper_path = flip(x_upper);
    z_upper_path = flip(z_upper);

    % Append the lower surface without duplicating the leading edge
    x_lower_path = x_lower(2:end);
    z_lower_path = z_lower(2:end);

    coordinates = [x_upper_path, z_upper_path;
                   x_lower_path, z_lower_path];
end
