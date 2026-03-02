% This experiment investigates the effect of various geometric and flow parameters on the zero-lift wave drag of a delta wing with a diamond-shaped airfoil. 
% Mach numbers
M = [1.5 1.75 2.0 2.25 2.5 3.0];
% Sweep angles
sweeps = [45, 60, 75];
% Define CD data for m = 0.5
CD_m05 = struct( ...
   't004', [ ...
       0.0068 0.0049 0.0040 0.0033 0.0026 0.0013; ...
       0.0043 0.0038 0.0035 0.0037 0.0032 0.0025; ...
       0.0026 0.0024 0.0023 0.0022 0.0021 0.0018], ...
   't006', [ ...
       0.0153 0.0111 0.0089 0.0074 0.0058 0.0029; ...
       0.0097 0.0085 0.0079 0.0084 0.0071 0.0055; ...
       0.0058 0.0055 0.0052 0.0049 0.0046 0.0042], ...
   't010', [ ...
       0.0424 0.0309 0.0248 0.0205 0.0162 0.0080; ...
       0.0270 0.0235 0.0219 0.0233 0.0198 0.0153; ...
       0.0160 0.0153 0.0144 0.0136 0.0128 0.0116]);
% Define CD data for m = 0.2
CD_m02 = struct( ...
   't004', [ ...
       0.0135 0.0085 0.0066 0.0053 0.0041 0.0018; ...
       0.0047 0.0051 0.0055 0.0071 0.0056 0.0041; ...
       0.0017 0.0017 0.0018 0.0019 0.0021 0.0023], ...
   't006', [ ...
       0.0303 0.0191 0.0148 0.0120 0.0093 0.0039; ...
       0.0106 0.0114 0.0123 0.0160 0.0126 0.0092; ...
       0.0039 0.0039 0.0041 0.0044 0.0047 0.0052], ...
   't010', [ ...
       0.0841 0.0531 0.0412 0.0334 0.0258 0.0109; ...
       0.0295 0.0316 0.0342 0.0444 0.0350 0.0255; ...
       0.0107 0.0109 0.0114 0.0121 0.0132 0.0145]);
% Helper function to plot
function plot_CD(M, CD_data, sweeps, titleText)
   colors = {'-or','-ob','-ok'};
   fields = fieldnames(CD_data);
   for i = 1:numel(fields)
       figure;
       hold on;
       data = CD_data.(fields{i});
       for j = 1:size(data,1)
           colorIdx = mod(j-1, numel(colors)) + 1; % wrap around if too many
           plot(M, data(j,:), colors{colorIdx}, 'DisplayName', sprintf('\\Lambda = %d°', sweeps(j)));
      end
       xlabel('Mach Number');
       ylabel('Drag Coefficient (C_D)');
       title([titleText ', \tau = ' fields{i}(2:end)]);
       legend show;
       grid on;
   end
end
% Plot for m = 0.5
plot_CD(M, CD_m05, sweeps, 'm = 0.5');
% Plot for m = 0.2
plot_CD(M, CD_m02, sweeps, 'm = 0.2');
function plot_K_vs_BcotLambda_combined(CD_struct, M, Lambda_deg, m_val)
   % Define color map for sweep angles
   colorList = {'r','b','k','g','m','c'};
   markerList = {'o','s','^','d','x','v'}; % different markers for tau
   Lambda_rad = deg2rad(Lambda_deg);
   cotLambda = cot(Lambda_rad);
   tanLambda = tan(Lambda_rad);
   B = sqrt(M.^2 - 1);
   fields = fieldnames(CD_struct); % 't004', 't006', etc.
   figure; hold on;
   legend_entries = {};
  
   for i = 1:numel(fields) % each tau
       % Extract tau from field name
       tau_str = fields{i}(2:end);
       tau = str2double(tau_str) / 100;
      
       data = CD_struct.(fields{i}); % size = numSweeps x numMach
      
       for j = 1:length(Lambda_deg) % each sweep angle
           CD_vals = data(j, :);
           K_vals = CD_vals .* tanLambda(j) / (4 * tau^2);
           x_vals = B .* cotLambda(j);
           colorIdx = mod(j-1, numel(colorList)) + 1;
           markerIdx = mod(i-1, numel(markerList)) + 1;
           % Combine color and marker for plotting
           style = [colorList{colorIdx} markerList{markerIdx} '-'];
           plot(x_vals, K_vals, style, 'LineWidth', 1.5);
           % Add to legend with both tau and sweep info
           legend_entries{end+1} = sprintf('\\tau = %.2f, \\Lambda = %d°', tau, Lambda_deg(j));
       end
   end
   xlabel('B cot(\Lambda)');
   ylabel('K');
   title(sprintf('K vs. B cot(\\Lambda), m = %.1f', m_val));
   legend(legend_entries, 'Location', 'bestoutside');
   grid on;
end
Lambda_deg = [45, 60, 75];
plot_K_vs_BcotLambda_combined(CD_m05, M, Lambda_deg, 0.5);
plot_K_vs_BcotLambda_combined(CD_m02, M, Lambda_deg, 0.2);
