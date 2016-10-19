%--------------------------------------------------------------------------
%                   Majorized Minimization (MoM)                          |
%-------------------------------------------------------------------------+
%   This function is designed to solve the optimization task of the       |
% Assisted Dictionary Learning (ADL) via Majorization Method also known   |
% as the Majorized Minimization algorithm.                                |
%                                                                         |
%-------------------------------------------------------------------------+
%     19 Oct 2016 - Vs 5.5                                         -mMm-  |
%-------------------------------------------------------------------------+
%        PARAMETERS:                                                      |
%   param.data  -> double T x N matrix with the data                      |
%   param.K     -> Number of components                                   |
%   param.iter  -> Total number of iterations                             |
%   param.Verb  -> (optional) Display low/off             (low by default)|
%                                                                         |
%     -Spatial Maps                                                       |
%   param.Ts    -> Number of iteration to compute the spatial maps        |
%   param.lambda-> The λ of the problem                                   |
%   param.Es    -> (optional) If it is selected, the stop criteria of     |
%          the spatial maps is changed to the error Es between iterations |
%          given by ||S^[n]-S^[n-1]||_{F} < Es                            | 
%          by default the maximum of iteration is 100                     |
%                                                                         |
%     -Dictionary                                                         |
%   param.Td    -> Number of iteration to compute the dictionary          |
%   param.ccl   -> Value of the normalization of each atom                |
%   param.Ed    -> (optional) The stop criteria is changed to use an      |
%          error check step ||D^[n]-D^[n-1]||_{F} < Ed                    |
%          by default the total number of iterations is 100               |
%                                                                         |
%    - Canonical Dictionary                                               |
%   param.cdl   -> Parameter of the proximity of the canonical atom on    |
%                     respect the current atom of the dictionary          |
%   param.Del   -> double T x M matrix which contains the canonical       |
%                     dictionary with the constrained atoms               |
%-------------------------------------------------------------------------+
%        RETURN:                                                          |
%   D  -> double T x K matrix  with the dictionary atoms                  |
%   S  -> double K x N matrix  with the coefficient                       |
%   E  -> A vector with the error per iteration (optional)                |
%--------------------------------------------------------------------------
function [D,S,E] = MoM(param) % Main Function

%%%%%%% PARAMETERS %%%%%%%%
	Y = param.data;         % Data
	[T,N] = size(Y);        % Number of time component and voxels
	K = param.K;            % Number of components
	lamb = param.lambda;    % Lambda parameter
    Tt = param.iter;        % Number of iterations


	% Normalization of the parameter λ
	lamb = lamb*sqrt(norm(Y,'fro')/(T*N));
	param.lambda = lamb;



	% Check Verbose
    if(~isfield(param,'Verb'))
		Verb=true;
    else 
        if(strcmp(param.Verb,'low'))
			Verb=true;
		else
			Verb=false;
        end
    end
    

%%%%% INTIALIZATION %%%%%
	D = 1-2*rand(T,K);
	S = zeros(K,N);
    E  = zeros(Tt,1);
    
    % Canonical Dictionary
    if(isfield(param,'Del'))
		[~,M]=size(param.Del);
        for i=1:M;
            D(:,i)=param.Del(:,i);
        end
    end
    
    

if(Verb);fprintf('Param [Ok] --- Prgs >> [                    ] 00.0 %%');end



%%=== Main Loop ===========================================================
	for t=1:Tt

        
		S = NewCoef(S,D,param);  % Update Coefficients
        
		D = NewDict(S,D,param);  % Update Dictionary

    
   		%  Error
    	E(t) = sqrt(norm(Y-D*S,'fro')/(T*N));


		% VERBOSE STUFF ---------------------------------------------------
        if (Verb)
            if(mod(t,ceil(Tt/100.00))==0)
            	vprct = t*100.00/Tt;
            	if(vprct>100);vprct=100;end
            	vrat  = ceil((t-1)*20.00/Tt);
                for i=1:28 ; fprintf('\b');end
                for i=1:vrat ; fprintf('=');end
                for i=vrat:19 ; fprintf(' ');end
                fprintf('] %4.1f %%',vprct);
            end
        end
    	%------------------------------------------------------------------

	end

	if(Verb);fprintf('\n   Completed! \\(^^) \n\n');end
end


%--------------------------------------------------------------------------
%        INTERNAL FUNCTIONS                                               |
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%     function NewCoef (So,D,param)                                       |
% calculates the new coefficients according to the MM algorithm           |
%--------------------------------------------------------------------------
function [S] = NewCoef(So,D,param)

	% Parameters needed
	Y = param.data;         % Data
	lam = param.lambda;     % Lambda parameter

	% Check the stop criteria
	if(isfield(param,'Es'))
		Ts = 500;
		Es = param.Es;
	else
		Ts = param.Ts;
		Es = 0;
	end

	% Initializations
	I = diag(ones(1,length(D'*D)));
	S = So;

    Dux = D'*D;
    cS  = max(sqrt(eig(Dux.'*Dux))); % Faster Spectral Norm
    
    % Constant
    DY = D'*Y/cS;
    Aq = I-D'*D/cS;


%   ----- Main Internal loop ----------------------------------
	Err = 1;
	t = 1;
    while (t<=Ts && Err>Es)

		A = DY+Aq*S;
        
        
        A = wthresh(A,'s',0.5*lam);
        

		% Control step
		if(isfield(param,'Es'))
			Err = norm(S-A,'fro');
		end

		% Actualization
		S = A;
		t = t+1;
    end
end
%--------------------------------------------------------------------------
%     function NewDict (S,Do,param)                                       |
% calculates the new Dictionary atoms according to the MM algorithm       |
%--------------------------------------------------------------------------
function [D] = NewDict(S,Do,param)

	% Parameters needed
	Y = param.data;            % Data
	ccl = param.ccl;           % Value of the normalization of each atom

	% Canonical Dictionary Selection (If it exist)
    if (isfield(param,'Del'))
        Del = param.Del;
        [~,M] = size(Del);
        cdl = param.cdl;
    else
        M = 0;
    end

    % Check stop criteria
    if (isfield(param,'Ed'))
    	Ed = param.Ed;
    	Td = 100;
    else
    	Td = param.Td;
    	Ed = 0;
    end


    % Initializations 
    I = diag(ones(1,length(S*S')));
    D = Do;

    Sux = S*S'; 
    cD  = max(sqrt(eig(Sux.'*Sux))); % Faster Spectral Norm

    
    % Constants
    YS = Y*S'/cD;
    Bq = I-S*S'/cD;


%   ----- Main Internal Loop ----------------------------------
	Err = 1;
	t = 1;
    while (t<=Td && Err > Ed)
        
        B = YS + D*Bq;

		% Projection
        
        % Closeness
        if (M>0)       % If there exist the Canonical Dictionary
            for j=1:M

                b = B(:,j);
                d = Del(:,j);
            
                if (norm(b-d)^2 > cdl)
                    B(:,j) = d+cdl*(b-d)/norm(b-d)^2;
                end
            end
        end
        
        % Normalization
        Kv = sqrt(sum(B.^2));
        
        Kv(Kv < ccl) = 1;         % Condition of normalization
        Kv(1:M) = 1;              % Avoid change the Canonical Atoms

        
        B = bsxfun(@rdivide,B,Kv); % Normalize the rest
             

		% Control step
		if(isfield(param,'Ed'))
			Err = norm(B-D,'fro');
		end

		% Actualizations
		D = B;
		t = t+1;
    end
end