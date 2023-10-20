classdef MLEs
    %   Maximum Likelihood Estimator functions
    %   Maximum Likelihood Estinators for monoexponential and translational
    %   diffusion.
    %   written by Joseph S. Beckwith
    %       Yusuf Hamied Department of Chemistry, Lensfield Road, 
    %       University of Cambridge, Cambridge, CB2 1EW, UK
        
    methods
        function [TrackDataPoints,RRTrackX,DtrueX,SigmaX,SigmaXError,RRTrackY,DtrueY,...
                SigmaY,SigmaYError,RRTrackZ,DtrueZ,SigmaZ,SigmaZError,min_varx_index,...
                min_vary_index,min_varz_index,DxTrans,DxTransError,DyTrans,DyTransError,...
                DzTrans,DzTransError,DxyzTrans,DxyzTransError,DcorrX,DcorrY,DcorrZ,N] = ...
                transMLE(~,X,Y,Z,TimeStep,nlags)
            % MLE of diffusion coefficient from Montiel, D.; Cang, H.; Yang, H.
            % Quantitative Characterization of Changes in Dynamical Behavior 
            % for Single-Particle Tracking Studies.
            % J. Phys. Chem. B 2006, 110 (40), 19763–19770.
    %%%%%%%%%fprintf('%s','Translational MLE... ')

            % Save the number of data points in the segment
            TrackDataPoints = length(X);

            % Number of overlapping time lags
            lags = length(X)-1;
            if lags > nlags
                lags = nlags;
            end

            % Preallocate arrays for apparent and corrected diffusion
            % coefficient MLEs, as well as number of time lags of a
            % given length
            DcorrX = zeros(lags,1);
            DcorrY = zeros(lags,1);
            DcorrZ = zeros(lags,1);
            N = zeros(lags,1);

            % Calculate apparent and corrected diffusion coefficients
            % using *overlapping* time lags
            parfor dn = 1:lags
                dX = X(1+dn:end) - X(1:end-dn);
                dY = Y(1+dn:end) - Y(1:end-dn);
                dZ = Z(1+dn:end) - Z(1:end-dn);
                N(dn) = length(dX);
                DcorrX(dn) = (mean(dX.^2)./(2*dn*TimeStep));
                DcorrY(dn) = (mean(dY.^2)./(2*dn*TimeStep));
                DcorrZ(dn) = (mean(dZ.^2)./(2*dn*TimeStep));
            end

            fitMin = 1;
            fitMax = 5000;
            if fitMax > nlags
                fitMax = nlags;
            end

            delta = (fitMin:fitMax)'*TimeStep;
            DcorrX = DcorrX(fitMin:fitMax);
            DcorrY = DcorrY(fitMin:fitMax);
            DcorrZ = DcorrZ(fitMin:fitMax);
            N = N(fitMin:fitMax);

            % Fits are performed for each axis, again according to the
            % procedure given by the reference: Xu, C.
            % et al. J. Phys. Chem. C 111, 32-35 (2007)
            s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[Inf, Inf],'Startpoint',[mean(DcorrX), 0.015],'Weights',N(fitMin:fitMax));
            f = fittype('a+((c^2)/x)','independent', 'x', 'dependent', 'y','options',s);
            [fitMLEx,gof,output] = fit(delta,DcorrX,f);
            RRTrackX = gof.rsquare; % R^2 of fit
            DtrueX = fitMLEx.a; % Diffusion coefficient MLE (um^2/s)
            SigmaX = fitMLEx.c; % Mobile localization precision (um)
            [~,Error] = nlparci([fitMLEx.a fitMLEx.c],output.residuals,'jacobian',output.Jacobian); % Standard errors of fitting parameters
            SigmaXError = Error(2); % (microns)

            s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[Inf, Inf],'Startpoint',[mean(DcorrY), 0.015],'Weights',N(fitMin:fitMax));
            f = fittype('a+((c^2)/x)','independent', 'x', 'dependent', 'y','options',s);
            [fitMLEy,gof,output] = fit(delta,DcorrY,f);
            RRTrackY = gof.rsquare;
            DtrueY = fitMLEy.a;
            SigmaY = fitMLEy.c;
            [~,Error] = nlparci([fitMLEy.a fitMLEy.c],output.residuals,'jacobian',output.Jacobian);
            SigmaYError = Error(2);

            s = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[Inf, Inf],'Startpoint',[mean(DcorrZ), 0.015],'Weights',N(fitMin:fitMax));
            f = fittype('a+((c^2)/x)','independent', 'x', 'dependent', 'y','options',s);
            [fitMLEz,gof,output] = fit(delta,DcorrZ,f);
            RRTrackZ = gof.rsquare;
            DtrueZ = fitMLEz.a;
            SigmaZ = fitMLEz.c;
            [~,Error] = nlparci([fitMLEz.a fitMLEz.c],output.residuals,'jacobian',output.Jacobian);
            SigmaZError = Error(2);

            % Corrected Variance of diffusion coefficient MLE
            % Variance of diffusion coefficient MLE
            var_DtrueX = 2*(SigmaX^2 + delta*DtrueX).^2./(N.*delta.^2);
            var_DtrueY = 2*(SigmaY^2 + delta*DtrueY).^2./(N.*delta.^2);
            var_DtrueZ = 2*(SigmaZ^2 + delta*DtrueZ).^2./(N.*delta.^2);

            % Find time lags which minimize the variance in the
            % diffusion coefficient MLEs
            [~,min_varx_index] = min(var_DtrueX);
            [~,min_vary_index] = min(var_DtrueY);
            [~,min_varz_index] = min(var_DtrueZ);

            % Calculate translational diffusion parameters and their
            % errors using the minimal variance time lags
            DxTrans = (DcorrX(min_varx_index) - SigmaX^2/delta(min_varx_index))*10^-12; % m^2/s
            DxTransError = (sqrt(var_DtrueX(min_varx_index)))*10^-12; % m^2/s
            DyTrans = (DcorrY(min_vary_index) - SigmaY^2/delta(min_vary_index))*10^-12; % m^2/s
            DyTransError = (sqrt(var_DtrueY(min_vary_index)))*10^-12; % m^2/s
            DzTrans = (DcorrZ(min_varz_index) - SigmaZ^2/delta(min_varz_index))*10^-12; % m^2/s
            DzTransError = (sqrt(var_DtrueZ(min_varz_index)))*10^-12; % m^2/s
            DxyzTrans = (DxTrans+DyTrans+DzTrans)/3; % m^2/s
            DxyzTransError = (1/3)*sqrt(DxTransError^2 + DyTransError^2 + DzTransError^2); % m^2/s

%%%%%%%%%%%%%%%%%%%fprintf('%s\n','done.')
        end  
        
    end
end

