%% Gaussian Process Regression Modeling for AgNW Experiment
    load('Round_4_data.mat')

%% Choose FOM (uncomment one of the following lines)
    %FOM_type = 'unitless';
    FOM_type = 'T10';

%% Choose kernel function (uncomment one of the following lines)
    %kernel_type = 'ardsquaredexponential';
    kernel_type = 'ardmatern52';
    %kernel_type = 'ardmatern32';

%% Choose hyperparameter optimization method (uncomment one of the following lines)
    hyperopt = 'exact';    % Maximize log likelihood of improvement (fitrgp default)
    %hyperopt = 'kfold'; k = 5;    % k-fold cross-validation, default k = 5
    %hyperopt = 'leaveout'; % leave-one-out cross-validation

%% Customize min and max of X domain if desired (uncomment one of the following lines)
    minX = min(X_raw); maxX = max(X_raw); % Use as default if no custom values are desired
    
%% Scale raw data to [0 1]
    Gsheet = 1./Rsheet;

    X_data = (X_raw - minX)./(maxX - minX);
    T_data = (transmit_percent - min(transmit_percent))./(max(transmit_percent) - min(transmit_percent));
    G_data = (Gsheet - min(Gsheet))./(max(Gsheet) - min(Gsheet));

%% Set up domain axes, n points per axis (n^3 total points in grid)
    n = 101;  % Number of grid points per axis
    x1 = linspace(0,1,n)';
    x2 = linspace(0,1,n)';
    x3 = linspace(0,1,n)';

    [X1,X2,X3] = ndgrid(x1,x2,x3);
    X = [X1(:), X2(:), X3(:)];

%% Gaussian Process Regression
    sigma0 = 0.05;  % Uncertainty in measured data, as fraction of full-scale

% Set up kernel hyperparameters
    switch hyperopt
        case 'exact'
            fitmethod = 'exact';
            kernparamT = []; % Enter best guess for initial kernelparam inside [],
            kernparamG = []; % Leave empty if using fitrgp default guess
        case 'kfold'
            fitmethod = 'none';
            kernparamT = cvhyperopt_3D(X_data,T_data,hyperopt,k,kernel_type,sigma0);
            kernparamG = cvhyperopt_3D(X_data,G_data,hyperopt,k,kernel_type,sigmaG0);
        case 'leaveout'
            fitmethod = 'none';
            kernparamT = cvhyperopt_3D(X_data,T_data,hyperopt,1,kernel_type,sigma0);
            kernparamG = cvhyperopt_3D(X_data,G_data,hyperopt,1,kernel_type,sigmaG0);
    end

% Gaussian process regression model for transmission data
    mdl_T = fitrgp(X_data,T_data,...
                 KernelFunction=kernel_type,...
                 Sigma=sigma0,ConstantSigma=true,...
                 KernelParameters=kernparamT,...
                 FitMethod=fitmethod);
    
    [mu_T,sd_T] = predict(mdl_T,X);

T_model_data = predict(mdl_T,X_data);

% Gaussian process regression model for conductance data
    mdl_G = fitrgp(X_data,G_data,...
                 KernelFunction=kernel_type,...
                 Sigma=sigma0,ConstantSigma=true,...
                 KernelParameters=kernparamG,...
                 FitMethod=fitmethod);
    
    [mu_G,sd_G] = predict(mdl_G,X);

G_model_data = predict(mdl_G,X_data);

%% Restore scaled values to original values
    X1 = X1 * (maxX(:,1) - minX(:,1)) + minX(:,1);
    X2 = X2 * (maxX(:,2) - minX(:,2)) + minX(:,2);
    X3 = X3 * (maxX(:,3) - minX(:,3)) + minX(:,3);
    X = [X1(:), X2(:), X3(:)];

    Xdomain = [minX(:,1),maxX(:,1);...
               minX(:,2),maxX(:,2);...
               minX(:,3),maxX(:,3)];

    x1 = x1 * (maxX(:,1) - minX(:,1)) + minX(:,1);
    x2 = x2 * (maxX(:,2) - minX(:,2)) + minX(:,2);
    x3 = x3 * (maxX(:,3) - minX(:,3)) + minX(:,3);

    mu_T = mu_T * (max(transmit_percent) - min(transmit_percent)) + min(transmit_percent);
    T_model_data = T_model_data * (max(transmit_percent) - min(transmit_percent)) + min(transmit_percent);
    sd_T = sd_T * (max(transmit_percent) - min(transmit_percent));

    mu_G = mu_G * (max(Gsheet) - min(Gsheet)) + min(Gsheet);
    G_model_data = G_model_data * (max(Gsheet) - min(Gsheet)) + min(Gsheet);
    sd_G = sd_G * (max(Gsheet) - min(Gsheet));

%% Calculate FOM from T and G models
switch FOM_type
    case 'unitless'
        mu_FOM = 188.5 * mu_G ./ ((mu_T/100).^-0.5 - 1);
        sd_FOM = FOM1stdev(mu_T,sd_T,mu_G,sd_G); % Std dev of FOM
    case 'T10'
        mu_FOM = (mu_T/100).^10.* mu_G;
        sd_FOM = FOM2stdev(mu_T,sd_T,mu_G,sd_G); % Std dev of FOM
end

%% Calculate quantitative regression quality parameters
    [slope_T,intercept_T,r2_T] = model_to_data_fit(transmit_percent,T_model_data);

    [slope_G,intercept_G,r2_G] = model_to_data_fit(Gsheet,G_model_data);

 %% Find extrema of regression models
    [mu_T_max,posT] = max(mu_T);
    X_max_T = X(posT,:);

    [mu_G_max,posG] = max(mu_G);
    X_max_G = X(posG,:);

    [mu_FOM_max,pos] = max(mu_FOM);
    X_max_FOM = X(pos,:);

%% Print regression parameters
    fprintf('Transmission Measured v. Modeled slope: %.5f\n', slope_T);
    fprintf('Transmission Measured v. Modeled intercept: %.3f\n', intercept_T);
    fprintf('Transmission Model r-squared: %.5f\n', r2_T);
    fprintf('Modeled transmittance max: µ_T_max(%.1f, %.0f, %.1f) = %.1f ± %.1f %%\n',...
            X_max_T(1), X_max_T(2), X_max_T(3), mu_T_max, sd_T(posT));
    fprintf('Modeled Resistance @ Transmittance max: %.2f ohm/sq\n', 1/mu_G(posT));

    fprintf('\nG Sheet Conductance Measured v. Modeled slope: %.5f\n', slope_G);
    fprintf('G Measured v. Modeled intercept: %.4f\n', intercept_G);
    fprintf('G Model r-squared: %.5f\n', r2_G);
    fprintf('Modeled G sheet max: µ_G_max(%.1f, %.0f, %.1f) = %.4f ± %.4f sq/ohm\n',...
            X_max_G(1), X_max_G(2), X_max_G(3), mu_G_max, sd_G(posG));
    fprintf('   (Modeled R Sheet Min = %.2f ± %.2f ohm/sq)\n',1/mu_G_max,sd_G(posG)/mu_G_max^2);
    fprintf('Modeled Transmittance @ Gsheet max: %.2f ohm/sq\n', mu_T(posG));
%%
    fprintf('\nModeled FOM max: µ_FOM_max(%.1f, %.0f, %.1f) = %.4f ± %.4f\n',...
            X_max_FOM(1), X_max_FOM(2), X_max_FOM(3), mu_FOM_max, sd_FOM(pos));
    fprintf('   Transmission at FOM max: = %.1f ± %.1f %%\n',mu_T(pos),sd_T(pos));
    fprintf('   Gsheet at FOM max: = %.4f ± %.4f sq/ohms\n',mu_G(pos),sd_G(pos));
    fprintf('     (Rsheet at FOM max: = %.2f ± %.2f ohms/sq)\n',1/mu_G(pos),sd_G(pos)/mu_G(pos)^2);
    

%% Plot regression results
    % Generate meshgrid from ndgrid coords; required for slice plotting
        X1mesh = permute(X1,[2,1,3]);
        X2mesh = permute(X2,[2,1,3]);
        X3mesh = permute(X3,[2,1,3]);
    % Generate slice planes for 3D rendering of posterior mean function
        x1slice = Xdomain(1,1):(Xdomain(1,2)-Xdomain(1,1))/20:Xdomain(1,2);
        x2slice = Xdomain(2,1):(Xdomain(2,2)-Xdomain(2,1))/20:Xdomain(2,2);
        x3slice = Xdomain(3,1):(Xdomain(3,2)-Xdomain(3,1))/20:Xdomain(3,2);

%% Plot regression results for transmittance
% Plot transmittance GPR posterior mean function on 3D volume heat plot
fig1 = figure(1); fig1.Position = [10 550 1000 400];
    subplot(121)
    mu_T_cube = reshape(mu_T,[n,n,n]);       % reshape uses ndgrid coords
    mu_T_cube_mesh = permute(mu_T_cube,[2,1,3]);  % Convert to meshgrid coords for slice plotting
    h1 = slice(X1mesh,X2mesh,X3mesh,mu_T_cube_mesh,x1slice,x2slice,x3slice);
    set(h1,'edgecolor','none','FaceAlpha',0.2)
    ax = gca; ax.YDir = 'reverse';
    xlim("tight"); ylim("tight"); zlim("tight");
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)'); zlabel('Volume (µL)');
    title('GPR Posterior Mean Function for Transmittance');
    cb = colorbar; colormap parula; cb.Label.String = 'Modeled Transmittance (%)';

% Plot measured vs modeled transmittance
    subplot(122);
    scatter(transmit_percent,T_model_data,25,'g','filled','MarkerEdgeColor','k');
    lsline;
    hold on
    title('Predicted vs. Modeled for Transmittance');
    xlabel('Measured Transmittance (%)'); ylabel('Modeled Transmittance (%)');
    legend('Data Points','Linear Fit','Location','northwest');
    
    annotation('textbox',[0.75, 0.2, 0.1, 0.1],...
        'String',...
        "Slope = " + slope_T + newline + ...
        "Intercept = " + intercept_T + newline + ...
        "Fit r^2 = " + r2_T);
    hold off;

%% Plot projected transmittance maxima on 2D planes
fig2 = figure(2); fig2.Position = [10 40 1400 300];

    M12 = max(mu_T_cube,[],3);
    M13 = max(mu_T_cube,[],2); M13 = squeeze(M13);
    M23 = max(mu_T_cube,[],1); M23 = squeeze(M23);
    Mmin = min([M12(:);M13(:);M23(:)]);

subplot(131);
[X1,X2] = ndgrid(x1,x2);
contourf(X1,X2,M12,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)');
    cb = colorbar; colormap summer; cb.Label.String = 'Modeled Transmittance Max';
    clim([Mmin mu_T_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,1),X_raw(:,2),20,'green','filled','MarkerEdgeColor','k');
scatter(X_max_T(:,1),X_max_T(:,2),50,'*','red');
legend('','Data Pts','Modeled T Max','Location','northeast');
hold off;
%
subplot(132);
[X1,X3] = ndgrid(x1,x3);
contourf(X1,X3,M13,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Volume (µL)');
    cb = colorbar; colormap summer; cb.Label.String = 'Modeled Transmittance Max';
    clim([Mmin mu_T_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,1),X_raw(:,3),20,'green','filled','MarkerEdgeColor','k');
scatter(X_max_T(:,1),X_max_T(:,3),50,'*','red');
legend('','Data Pts','Modeled T Max','Location','northeast');
hold off;
%
subplot(133);
[X2,X3] = ndgrid(x2,x3);
contourf(X2,X3,M23,50,'FaceColor','flat','LineStyle','none');
    xlabel('Spin Speed (rpm)'); ylabel('Volume (µL)');
    cb = colorbar; colormap summer; cb.Label.String = 'Modeled Transmittance Max';
    clim([Mmin mu_T_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,2),X_raw(:,3),20,'green','filled','MarkerEdgeColor','k');
scatter(X_max_T(:,2),X_max_T(:,3),50,'*','red');
legend('','Data Pts','Modeled T Max','Location','northeast');
hold off;

sgtitle('Projected Transmittance Maxima')
%
%% Plot regression results for sheet resistance
% Plot GPR posterior mean function on 3D volume heat plot
fig3 = figure(3); fig3.Position = [100 550 1000 400];    
    subplot(121)
    mu_G_cube = reshape(mu_G,[n,n,n]);
    mu_G_cube_mesh = permute(mu_G_cube,[2,1,3]); % Convert to meshgrid coords for plotting
    h2 = slice(X1mesh,X2mesh,X3mesh,mu_G_cube_mesh,x1slice,x2slice,x3slice);
    set(h2,'edgecolor','none','FaceAlpha',0.2)
    ax = gca; ax.YDir = 'reverse';
    xlim("tight"); ylim("tight"); zlim("tight");
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)'); zlabel('Volume (µL)');
    title('GPR Posterior Mean Function for Gsheet');
    cb = colorbar; colormap turbo; cb.Label.String = 'Modeled G (sq/ohm)';

% Plot modeled vs measured sheet conductance
    subplot(122);
    scatter(Gsheet,G_model_data,25,'g','filled','MarkerEdgeColor','k');
    lsline;
    xlim("tight"); ylim("tight"); zlim("tight");
    hold on;
    title('Predicted vs. Modeled for Sheet Conductance');
    xlabel('Measured Conductance (sq/ohm)'); ylabel('Modeled Conductance (sq/ohm)');
    legend('Training Pts','Linear Fit','Location','northwest');
    
    annotation('textbox',[0.75, 0.2, 0.1, 0.1],...
        'String',...
        "Slope = " + slope_G + newline + ...
        "Intercept = " + intercept_G + newline + ...
        "Fit r^2 = " + r2_G);
    hold off;

%% Plot projected G maxima on 2D planes
fig4 = figure(4); fig4.Position = [100 40 1400 300];

    M12 = max(mu_G_cube,[],3);
    M13 = max(mu_G_cube,[],2); M13 = squeeze(M13);
    M23 = max(mu_G_cube,[],1); M23 = squeeze(M23);
    Mmin = min([M12(:);M13(:);M23(:)]);

subplot(131);
[X1,X2] = ndgrid(x1,x2);
contourf(X1,X2,M12,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)');
    cb = colorbar; colormap cool; cb.Label.String = 'Modeled G Max (sq/ohm)';
    clim([Mmin mu_G_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,1),X_raw(:,2),20,'green','filled','MarkerEdgeColor','k');
scatter(X_max_G(:,1),X_max_G(:,2),50,'*','red');
legend('','Data Pts','Modeled G Max','Location','southeast');
hold off;
%
subplot(132);
[X1,X3] = ndgrid(x1,x3);
contourf(X1,X3,M13,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Volume (µL)');
    cb = colorbar; colormap cool; cb.Label.String = 'Modeled G Max (sq/ohms)';
    clim([Mmin mu_G_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,1),X_raw(:,3),20,'green','filled','MarkerEdgeColor','k');
scatter(X_max_G(:,1),X_max_G(:,3),50,'*','red');
legend('','Data Pts','Modeled G Max','Location','southeast');
hold off;
%
subplot(133);
[X2,X3] = ndgrid(x2,x3);
contourf(X2,X3,M23,50,'FaceColor','flat','LineStyle','none');
    xlabel('Spin Speed (rpm)'); ylabel('Volume (µL)');
    cb = colorbar; colormap cool; cb.Label.String = 'Modeled R min (ohms/sq)';
    clim([Mmin mu_G_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_raw(:,2),X_raw(:,3),10,'green','filled','MarkerEdgeColor','k');
scatter(X_max_G(:,2),X_max_G(:,3),50,'*','red');
legend('','Data Pts','Modeled G Max','Location','southeast');
hold off;

sgtitle('Projected Gsheet Maxima');
%
%% Plot projected unitless FOM maximums on 2D planes
fig5 = figure(5); fig5.Position = [10 40 1400 300];

    mu_FOM_cube = reshape(mu_FOM,[n,n,n]);  % reshape uses ndgrid coords
    mu_FOM_cube(mu_FOM_cube > 200) = 200;
    M12 = max(mu_FOM_cube,[],3);
    M13 = max(mu_FOM_cube,[],2); M13 = squeeze(M13);
    M23 = max(mu_FOM_cube,[],1); M23 = squeeze(M23);
    Mmin = min([M12(:);M13(:);M23(:)]);

subplot(131);
[X1,X2] = ndgrid(x1,x2);
contourf(X1,X2,M12,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Spin Speed (rpm)');
    cb = colorbar; colormap bone; cb.Label.String = 'Modeled FOM';
    clim([Mmin mu_FOM_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_max_FOM(:,1),X_max_FOM(:,2),80,'*','red');
legend('','Modeled FOM Max','Location','northeast');
hold off;
%
subplot(132);
[X1,X3] = ndgrid(x1,x3);
contourf(X1,X3,M13,50,'FaceColor','flat','LineStyle','none');
    xlabel('Concentration (mg/mL)'); ylabel('Volume (µL)');
    cb = colorbar; colormap bone; cb.Label.String = 'Modeled FOM';
    clim([Mmin mu_FOM_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_max_FOM(:,1),X_max_FOM(:,3),80,'*','red');
legend('','Modeled FOM Max','Location','northeast');
hold off;
%
subplot(133);
[X2,X3] = ndgrid(x2,x3);
contourf(X2,X3,M23,50,'FaceColor','flat','LineStyle','none');
    xlabel('Spin Speed (rpm)'); ylabel('Volume (µL)');
    cb = colorbar; colormap bone; cb.Label.String = 'Modeled FOM';
    clim([Mmin mu_FOM_max]);  % Use same colorbar min and max for all 3 plots
hold on;
scatter(X_max_FOM(:,2),X_max_FOM(:,3),80,'*','red');
legend('','Modeled FOM Max','Location','northeast');
hold off;

sgtitle('Projected FOM Maxima');
%
%
%% Function to find linear fit parameters of model to data
function [slope,intercept,r2] = model_to_data_fit(real_data,model_data)
    linfit = LinearModel.fit(real_data,model_data); % Fit real data vs modeled data to a line
    Coefficients = table2array(linfit.Coefficients); % Extract fit parameters
    slope = Coefficients(2,1);
    intercept = Coefficients(1,1);
    r2 = linfit.Rsquared.ordinary;
end
%
%% Function to calculate standard deviation of unitless FOM
function sd_FOM = FOM1stdev(mu_T,sd_T,mu_G,sd_G)
    T = mu_T/100;
    sd_T = sd_T/100; % stdev of T in decimal given stdev of T in %

    fT = (T.^-0.5 - 1).^-1; % Expectation value of 1/(T^-0.5 - 1) given T
    fT_sd = 0.5*((T.^0.5 -1).^2.*T.^0.5).^-1.*sd_T; % stdev of 1/(T^-0.5 - 1) given T, stdev of T

    sd_FOM = 188.5*sqrt((fT.*sd_G).^2 + (mu_G.*fT_sd).^2 + (sd_G.*fT_sd).^2);
end
%
%% Function to calculate standard deviation of T^10 FOM
function sd_FOM = FOM2stdev(mu_T,sd_T,mu_G,sd_G)
    T = mu_T/100;
    sd_T = sd_T/100; % stdev of T in decimal given stdev of T in %

    fT = T.^10; % Expectation value of T^10 given T
    fT_sd = 10*(T.^9).*sd_T; % stdev of T^10 given T, stdev of T

    sd_FOM = sqrt((fT.*sd_G).^2 + (mu_G.*fT_sd).^2 + (sd_G.*fT_sd).^2);
end
%