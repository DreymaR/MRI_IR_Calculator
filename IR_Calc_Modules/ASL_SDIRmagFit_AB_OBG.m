function ASL_SDIRmagFit_AB_OBG()                         % Mt = DIRmagFitASL(nInv)
%% From Atle B: Calculate and plot optimal dual TI for background suppression (ASL)
% Uses ga(fitnessfunction,nvars,LB,UB) genetic algorithm to find multiple TI.
% This is the main script - run this one. It calls magFunc for fitting and magPlot for plotting.
% NOTE: This script may end at slightly different TI values each run.
% GadOE: Edits and fiddlings (quite a lot of them)
% 
% DONE: - Reformulate times from TI1+TAG and TI2-PLD to generic DIR values. Use the former only in reporting.
%       - Move everything into this file as nested functions (avoids code repetition in the CalcMagZ function)
% 
% TODO: - Button to include readout, showing how suppression is lost slice by slice. Specify slice time/number?
%       - Vectorize the fitness function for speed? Then use 'UseVectorized',true in the gaoptimset?

%% INIT
% clear all                                           % Start without any variables predefined (good for debug).
global NInv TI0 DUR T1s NT1s                        % Inversion recovery and tissue/T1 parameters
global TIn MagZ Mtot                                % Results
saveFigureImage = false;                            % Save figure as image in Figures dir?

sSz = get(0,'ScreenSize');                          % Find the screen size
fSz = uint16(sSz(4)/10.00);                         % Scaling factor for figures (was 100)
figName = 'ASL Background Sat. TI time calculator'; % A (hopefully unique) name for this figure
myFig = findobj( 'Type', 'Figure',              ... % Find this script's figure (possibly among others)
        'Name', figName                         );
if isempty(myFig)                                   % Iff my figure does not yet exist...
    myFig = figure(                             ... % ...make a new figure!
        'Name', figName,                        ...
        'Position', [7*fSz 1*fSz 8*fSz 6*fSz] );    % At pos. (70%,10%) of ScreenSize, etc.
else
    figure(myFig);                                  % If my figure exists, make it the CurrentFigure
end

NInv = 2;                                           % Number of inversion pulses. NB: Keep this to 2 for now!
TAG = 1650; PLD = 1500;                             % Standard values for LD (label duration) and PLD, in ms
% inp = inputdlg('Enter Label Duration',   'TAG', 1, {'1650'});   % TODO: Make GUI UI instead!
% TAG = str2double(inp{1});
% inp = inputdlg('Enter Post Label Delay', 'PLD', 1, {'1500'});
% PLD = str2double(inp{1});                         % Pulse labeling delay in ms; change as needed
TI0 = TAG; DUR = TAG + PLD;                         % Lowest and highest possible inversion time
T1s   = [791, 1445, 4163];                          % T1 times to be suppressed. Only 3 values allowed for now?
T1str = [ "WM ", "GM ", "CSF" ];                    % Names of tissues corresponding to T1s above
% T1c = { "WM ",   791;           ...                 % Cell array of T1 labels and times. Looks cool, but gets complex.
%         "GM ",  1445;           ...                 % E.g., T1str(i) = string(T1c(i,1))
%         "CSF",  4163            };
NT1s = length(T1s);

%% MAIN
% Display presets, quick-n-dirty:
disp([ 'TAG(TImin) = ', num2str(TAG), ' ms; ',  ...
       'PLD(delTI) = ', num2str(PLD), ' ms; ',  ...
       'DUR(TImax) = ', num2str(DUR), ' ms'     ]);
fprintf( '\n' );

TIn = DIRmagFitFunc();                              % The optimization proper. Needs access to [ NInv, TI0, DUR, T1s ].

PlotMagZ();                                         % Plot the Z magnetization under inversions for all TI/T1
ASLPlotDeco();                                      % Decorate the plot with ASL relevant stuff

MagZ = CalcMagZ( TIn ); Mtot = mean(abs( MagZ ));   % Calculate the end Z magnetization. Not used: Mtot=Mz(:,DUR);
MoreOutput( saveFigureImage );                      % Output results

%% FUNCTIONS
% Fitting function for minimizing the total (weighted) magnetization Mt after NInv # of inversions at TI ms
function TInv = DIRmagFitFunc()     % global NInv TI0 DUR; cost function also needs T1s
    % Set initial TI-values and constraints for genetic algorithm (ga) fitting
    Tmin        = TI0 + 10;         % Is minimum TI 10 ms after end of label period too short? Maleki needed 20 ms.
    TIini(1)    = Tmin + 30;        % The starting value for ga() optimization
    TImin(1)    = Tmin;
    TImax(1)    = Tmin + 200;       % This number is arbitrary; could that be a problem? (Was 300)
    TIini(NInv) = DUR - 200;        % Was PLD - 200
    TImin(NInv) = TImax(1)+100;     % Minimum gap between inversion pulses = 100 ms..
    TImax(NInv) = DUR - 10;         % Was PLD - 10

    MzFn = @(TI)DIRmagCostFunc(TI); % Anon. handle to the magnetization function Mz(TI) (for passing the function to ga)

    % Enable these lines for constrained minimization using genetic algorithm:
    %-----------------------------------------------
    fprintf('Run 1... ');
    options = gaoptimset('InitialPopulation',TIini,'TolFun',.1e-8,'Generations',1000);
    [TInv,~,~,output] = ga(MzFn,NInv,[],[],[],[],TImin,TImax,[],options);     % [TI,fval,exitflag,output]
    TIini = TInv;
    fprintf('Run 2... ');
    options = gaoptimset('InitialPopulation',TIini,'TolFun',.1e-8,'Generations',2000);
    [TInv,~,~,output2] = ga(MzFn,NInv,[],[],[],[],TImin,TImax,[],options);
    %-----------------------------------------------

    % Enable these lines for unconstrained minimization (no TImin/TImax; won't work for ASL if TI1 < TAG):
    %-----------------------------------------------
    % options=optimset('TolX',.0005,'TolFun',1e-10,'MaxIter',1000);
    % [TInv,~,~,output] =fminsearch(Mz,TI0,options);
    %-----------------------------------------------

%     disp( output );                                 % Output info from the optimization (ga or fminsearch)
    TInv = round(TInv);                             % The array of found inversion times (now from t = 0)
end

% Fitness cost function for the total (weighted) magnetization Mt after inversions at TI ms
function Mt = DIRmagCostFunc( TI )
    Mz = CalcMagZ( TI );                            % ( TI, DUR, T1s )
    Mt = 0; pFact = 10;                             % Add penalty weighted by pFact on negative Mz-values
    for i = 1:NT1s
        if Mz(i)<0
            Mt = Mt + pFact * abs(Mz(i));
        else
            Mt = Mt + Mz(i);
        end
    end
end

% Calculate the relative Z magnetization Mz after [TI] inversions
function MagZ = CalcMagZ( TI )                      % ( TI, DUR, T1s )
%     NInv = length(TI);                              % # of inversion time points
    TID = [ TI, DUR ];
    dT = zeros(1,NInv+1);                           % Vector of the differences between TIs
    dT(1) = TID(1);                                 % Time point for the first inversion (was TAG + TI1)
    for i = 1:NInv
         dT(i+1) = TID(i+1) - TID(i);               % dT(2) = TI(2) - TI(1); dT(3) = DUR - TI(2)
    end

%     NT1s = length(T1s);                             % # of T1 times for which to calculate Mz
    MagZ = zeros(1,NT1s);                           % End Z magnetization for each T1
    for i = 1:NInv+1
%         MagZ = -MagZ.*exp(-dT(i)./T1s) + 1 - exp(-dT(i)./T1s);
        MagZ = 1 - (1 + MagZ).*exp(-dT(i)./T1s);
    end
end

% Plot the magnetization Mz of several T1s under multiple inversions
function PlotMagZ()                                 % ( TIn, DUR, T1s )
    dT = 1; Tp = 0:dT:DUR-1;                        % dT is the sampling interval (ms), Tp all time points
    Mz = zeros( NT1s+1, DUR );                      % MagZ for all tissues (plus a zero line), for all time points
    Mz(:,1) = 0;                                    % Starting Mz for each T1. Zero here.
    for j = 2:DUR
        t = (j-1)*dT;
        if( any( TIn(:) == t ))                     % Inversion time, so invert the Mz (was TI+TAG)
            Mz(:,j) = -Mz(:,j-1);
        else                                        % Signal evolvement ( if(t~=TI(:)+TAG) )
            for i = 1:NT1s
                Mz(i,j) = Mz(i,j-1)*exp(-dT/T1s(i)) + (1-exp(-dT/T1s(i)));
            end
        end
    end
    plot( Tp, Mz,                                   ...
            'LineWidth', 1.5                        )
end

% Decorate the plot with ASL relevant stuff
function ASLPlotDeco()
    figtxt = [  'ASL BgSuppr. optimization for ',   ... % Figure title
                'LD = ',  int2str(TAG), ' ms, ',    ... % Labeling Duration
                'PLD = ', int2str(PLD), ' ms:'      ];  % Post Label Delay
    title( figtxt,                                  ...
        'FontSize', 14                              );
    
    for i = 0.005:0.010:0.040                           % Draw this dotted line several times for an "RF comb" look
        line( [ 1 TAG ], [ i i ],                   ... % Horizontal line to represent ASL tagging
            'Color', [ 0.9 0.3 0.3 ],               ... % 'black'
            'LineWidth', 3,                         ... % Tip: Line width affects dot distance
            'DisplayName', ' ASL labeling',         ... % Tip: If declared after legend, gets in it
            'LineStyle', ':'                        );
    end
    
    LegS = strings(1,NT1s);                             % Figure legends
    for i = 1:NT1s
        LegS(i) = sprintf('%4s (T1 = %4i ms)', T1str(i), T1s(i) );  %legend('WM', 'GM', 'CSF'...
    end
    legend( LegS(1), LegS(2), LegS(3), 'Location', 'SouthWest');
    
    %TI1 = num2str(TI(1)); TI2 = num2str(TI(2));
    %TI_values = text(100,-0.5, sprintf('{\\itTI-values: %d %d }', TI1, TI2));
    % text(100,-0.2, '\\itInversion Pulses:')
    % text(100,-0.3, sprintf('\\bfTI_1 = %4.d ms',TI(1)))
    % text(100,-0.4, sprintf('\\bfTI_2 = %4.d ms',TI(2)))
    tStr = ['\\itInversion pulses:\\rm\n'           ... % Figure info text box
            '\\bfTI_1 = %4.d ms\n'                  ...
            '\\bfTI_2 = %4.d ms\\rm'                ];  % \n'              ...
%             '\\itEnd Z magn.:\\rm\n'                ...
%             '\\bfMtot = %4.d (a.u.)'                ];
    tVar = [ TIn(1), TIn(2) ];                          % , Mtot ];
    figtxt = sprintf( tStr, tVar );                     % Formatted string "print"
    xlm = xlim; xlm = xlm(2);                           % Max value of axis; (x,y) is a plot point!
    ylm = ylim; ylm = ylm(2);                           % --"--
    text( 0.800*xlm, -0.800*ylm, figtxt,            ...
        'EdgeColor', 'black',                       ...
        'Margin', 4,                                ...
        'LineWidth', 0.5,                           ...
        'LineStyle', '-'                            );
end

line( [ 1 TAG ], [ 0.005 0.005 ],               ... % Horizontal line to represent ASL tagging
    'Color', [ 0.9 0.3 0.3 ],                   ... % 'black'
    'LineWidth', 3,                             ... % Tip: Line width affects dot distance
    'DisplayName', ' ASL labeling',             ... % Tip: If declared after legend, gets in it
    'LineStyle', ':'                            );

% Quick-n-dirty results output
function MoreOutput( saveFigImg )                       % ( TAG, PLD ), ( MagZ, Mtot, TI )
    % Save the figure as an image (this can be done from its window menu too)
    if saveFigImg
        if exist('figurer','dir')~=7,mkdir('figurer'),end
        cd figurer
        print('-djpeg','-r300',['PLD_' int2str(PLD) '_and_LD_'  int2str(TAG) '-InvTider'])
        cd '..'
    end

    % Display answers, quick-n-dirty:
    fprintf(  '\nOptimization info:\n' );
    disp([ 'Mz(T1s) = ',    num2str(MagZ)           ]);
    disp([ 'MtotW   = ',    num2str(Mtot)           ]);
    disp([ 'TI1,TI2 = ',    num2str(TIn(1)),' ms, ', num2str(TIn(2)),' ms' ]);
end

end
