%% Pareto Front Modeling Script
% In Workspace have results of running 'GPR_FOM_script.m'
close all;

%% Calculate Pareto front of regression model
    % Define the two-objective function for the Pareto search
    func = @(x)[-predict(mdl_T,x),-predict(mdl_G,x)];
        % paretosearch is min finder, since we want max T and G 
        % we use negative of models
    % Define lower and upper bounds on the domain X
    lb = [0 0 0]; ub = [1 1 1];

    % Calculate (scaled) Pareto front of regression model
    options = optimoptions(@paretosearch,'Display','off'); % Suppresses messages

    [X_pareto_calc,Y_pareto_calc] = ...
        paretosearch(func,3,[],[],[],[],lb,ub,[],options);

    % Unscale values back to experimental units
    X_pareto_calc = X_pareto_calc .* (maxX - minX) + minX;
    Y_pareto_calc(:,1) = -Y_pareto_calc(:,1) * (max(transmit_percent)-min(transmit_percent))...
        + min(transmit_percent); % Lead minus undoes minus in func
    Gsheet = 1./Rsheet;
    Y_pareto_calc(:,2) = -Y_pareto_calc(:,2) * (max(Gsheet)-min(Gsheet)) + min(Gsheet);

    Y_pareto_calc(:,2) = 1./Y_pareto_calc(:,2);  % Convert back to resistance
    
    % Remove rows with non-positive values in modeled conductance
    idx = find(Y_pareto_calc(:,2) <= 0);
    Y_pareto_calc = removerows(Y_pareto_calc,idx);
    X_pareto_calc = removerows(X_pareto_calc,idx);

%% Plot transmission vs Rsheet data and those points on Pareto front
    fig1 = figure(1); fig1.Position = [10 550 550 400];
    semilogy(transmit_percent,Rsheet,'ko');
    hold on
    % Plot regression model Pareto front
    scatter(Y_pareto_calc(:,1),Y_pareto_calc(:,2),'blue','square');
    xlabel('Transmission (%)'); ylabel('Rsheet (ohm/sq)');
    legend('Data','Regression Model Pareto Front',...
        'Location','northwest');
    title('Pareto Front');
    hold off

%% Plot X-values used for data and those corresponding to points on 
        % Regression Pareto front
    fig2 = figure(2); fig2.Position = [600 550 550 400];
    scatter3(X_raw(:,1),X_raw(:,2),X_raw(:,3),'ko');
    hold on
    set(gca,'YDir','reverse') % Reverses direction of y-axis
    scatter3(X_pareto_calc(:,1),X_pareto_calc(:,2),...
        X_pareto_calc(:,3),'blue','square');
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)');
        zlabel('Volume (µL)');
    legend('Data Inputs','Modeled Pareto Front Inputs',...
        'Location','northwest');
    title('Experimental Input Parameters');
    hold off
%