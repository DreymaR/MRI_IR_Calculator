function DIR_TI_Calculator_v20171221()
%% DIR MRI TI/T1 calculation
%   - In mode 1 & 2, find TI time(s) that null given T1(s) in DIR & 1IR (FLAIR/STIR/etc).
%   - In mode 3, calculate/plot relative signal strength for WM & WML as a function of TI2 in a DIR sequence,
%       - find TI2 and TI1(TI2) such that there is no T1 weighting between WM & WML while nulling CSF,
%       - and find a T1 that's nulled at the resulting TIs, corresponding to a fictive tissue.
% 
%  NOTE:
%   - TI times are defined vendor-style here! That is, TI1 is the full time and TI2 only the time to readout.
%       - So TI1 = time_TI1 & TI2 = time_TI1 - time_TI2. This makes the signal formulae simpler.
%       - Used to call the time between inversions TI1. Now it's TI1_new = TI1_old + TI2.
%   - Don't make signal equations generally symbolic rather than functions, as symexpr are slow (esp. in r2014a?)!
%   - Anonymous function handle for CalcTI1, for easier passing to other functions such as vpasolve()? Likely not.
%   - There's a good explanation of inversion mathematics at http://xrayphysics.com/contrast.html
%   - About (single) incomplete inversion (OliG: 95% is realistic), see Miho Kita et al, MRI 31(9) 2013 p 1631–1639.
%       - http://www.sciencedirect.com/science/article/pii/S0730725X13002312
% 
%  TODO:
%   - Include AtleB's ASL BS TI optimalization? Use a press button for it, as it takes some time to run.
%   - Solve Mz(0) in TIR/DIR by formula [1 - E_(TR-TI1)] instead of iteration? But it works well as it is.
%   - Correct the generic S0 formula CalcS0(). It gives very wrong results! (See the debug for DIR.)
%   - Implement incomplete spin lock during readout! CPMG inversion pulses may be 120-150° instead of 180°.
%       - So readout CPMG inversion efficiency is then -cos(FA)=50-87%
% 
%  DONE:
%   - Find DIR TI1 & TI2 such that two given T1 times are nulled.
%       - WN: Then we can see any discrepancies compared to vendor (GE/Sie.) calculations!
%   - Switch between T1W-nulling mode and T1 tissue nulling mode (and figs?).
%   - GUI controls for tissue T1s. B0 settings update them, but also input.
%   - Single IR (STIR/FLAIR/etc) and DIR tissue nulling.
%   - TIR/DIR Mz(t) curves w/ spin locked readout starting at TI.
%       - Starting Mz(0) of these solved by iterating over a number of TR (~5–8).
%       - Realized that readout should start at TI where the nulling is.
%       - BUT: What if there's a wait before the echo train? Not accounted for. Haven't seen the sequence diagrams.
%   - Solve optimal TI1/TI2 symbolically instead of numerically (taking the minimum of the difference)
%   - Implement excitation 
%       - Without rewinding? What is the case in our sequences?
%       - Add a GUI for effective FA (90° for full excitation, 0° if refocused). Maybe also add rewinder?
%   - Implement incomplete inversion using the Inversion Efficiency = -cos(FA_inv)

%% INIT
saveImg = false;                                    % Save figure as image in Figures dir? (Can use the menu instead)
debugInfo  = true; %false;                          % Show extra info in the command window output

% sSz = get(0,'ScreenSize');                          % Get the screen size (actually, in R2015b+ pixels are 1/96")
fSz = 10.80; % sSz(4)/100.00;                       % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])    % Default fig. position & size
wh = [ 1 1 1 ]; figClr = 0.94*wh; myGray = 0.50*wh; % Colors of fig. background and grayed-out UI elements
% fig1 = get(groot, 'CurrentFigure');                 % Get handle of the current fig. ('gcf' creates one)
figName = 'DIR inversion time calculator';
myFig = findobj( 'Type', 'Figure',              ... % Find this script's figure (possibly among others)
        'Name', figName                         );
if isempty(myFig)                                   % Iff my figure does not yet exist... (was ~ishghandle)
    myFig = figure(                             ... % ...make a new figure!
        'Name', figName,                        ...
        'Color', figClr,                        ...
        'Position', [ 80 20 80 60 ]*fSz         );  % At pos. (80,20)% of ScreenSize, etc.
%         'NumberTitle', 'off',                   ... % Removes 'Figure #: ' from the window title
else
    figure(myFig);                                  % If my figure exists, make it the CurrentFigure
end % if

global RunMode pStr pVar                                % Program control
global UIon uiBox uiIES uiFAsS                          % UI elements
global uiTn1T uiTn2T uiTn2S uiSTnT uiSTnS uiB0T uiB0S   % --"--
global sTR sT_Ro sTR_Eff IEf IE IE2 FA REW              % System/sequence variables
global T1sel T1str T1tag T1set T1z T1plt T1txt          % T1 times/labels
global B0sel B0str B0tag                                % B0 settings/labels (System B0 and reference values)
global T1_CSF T1_WM T1_MS T1_GM T1_Fat                  % Tissue relaxation times (based on "B0")
global T2_CSF T2_WM T2_MS %T2_GM isn't known?           % --"--
global TI2_max T1_n1 T1_n2 TI1_n TI2_n SRo MzTI MzIt    % TI calculation, T1 tissue nulling, magnetization

t = sym('t');                                       % Set the symbolic explicitly as it's used in nested functions
% TI2s = sym('TI2s'); TI1s(TI2s) = TI2s;              % Temporary assignment globalizes the sym expression
% DIR_S0(t) = ( 1 - 2*exp(-TI2/t) + 2*exp(-TI1/t) - exp(-TReff/t) ); DIR S0(T1) as a sym expression
% TI1s(TI2s) = T1*log(2/(2*exp(-TI2s/T1) + exp(-sTR_Eff/T1) - 1)); % TI1(TI2) nulling T1_n1, using syms
% % Note: Do not replace CalcTI1(TI2,T1) with TI1s(t) above: It's much slower.

% Parameter presets for the first run (most can be changed in the UI)
if isempty(sTR)    ; sTR     = 6000.0 ; end         % Repetition time (ms); as float
if isempty(sT_Ro)  ; sT_Ro   =  780.0 ; end         % Readout time (ms) ETL * ESP (GE: Check CVs after Download)
if isempty(IEf)    ; IEf     =   95   ; end         % Inversion efficiency in %
FA =  90;           %MEX = cos(FA*pi/180);          % MEX is MagZ after excitation with angle FA
REW = 0.00;                                         % Rewinder after readout, regaining some of the pre-readout Mz?

UIon = true;                                        % Show the UI control panel?
if isempty(RunMode); RunMode =    1   ; end         % Start in mode 1, then keep the run mode over reruns
if isempty(T1sel)  ; T1sel   =    3   ; end         % Start with a value (CSF), then keep the selection over reruns
setRelaxTs( 3 );                                    % Select 3T setup, with relaxation times from literature (GE sett.)
SetMode( RunMode );                                 % Note: Need to set RunMode -> T1sel -> RelaxTs -> SetMode.

%% MAIN
main();
function main()
    sTR_Eff = sTR - sT_Ro;                          % Recalculate based on current parameter settings
    IE = IEf/100; IE2 = 1 + IE;                     % MagZ to recover after incomplete inv. (was log(2/...) etc)
    T1 = T1_n1;                                     % Find parameters that null the longest T1 (usually CSF)
    TI2_max = log(IE2/(1 + IE*exp(-sTR_Eff/T1)))*T1; % Max TI2 nulling T1_n1 (at max TI1 = TReff)
%     TI2_max = T1*log(2/(1 + 2*exp(-sTR/T1) - exp(-sTR_Eff/T1))); % Old TI2_max calculation without InvEff
    
    clf                                             % Clear the figure so we can change it
    switch RunMode
        case 1                                      % Mode 1 (1IR tissue nulling)
            mainSIRNul();
        case 2                                      % Mode 2 (DIR tissue nulling)
            mainDIRNul();
        case 3                                      % Mode 3 (T1W nulling)
            mainT1WNul();
    end
    createUI();
    if saveImg, saveFigAsImg(),end                  % Automatically save the figure as an image (or use the fig. menu)?
end % main

% 1) Calculate TI in a 1IR/STIR/FLAIR sequence, nulling a specified T1 time
function mainSIRNul()
    T1z = T1_n1;                                    % What T1 are we concerned with nulling?
    TI1_n = log(IE2/(1 + IE*exp(-sTR_Eff/T1z)))*T1z; % Calc. 1IR TI nulling T1
%     TI1_n = log(2)*T1z;                             % The simple case, TR>>T1, gives TIn = ln(2)*T1 ? 0.69*T1
    TI1_n = uint16( TI1_n );                        % Make value int for figure display
    
    TInts = [ 0, TI1_n, sTR ];                      % SIR: Inversion is at t = 0, then readout at TI
    Mzt = PlotMagZ( TInts );                        % Plot Z magnetization over time, return end MagZ vector
    
    figtxt = ['Single IR tissue nulling at ' B0tag ', for TR_{eff} = ' num2str(sTR_Eff) ' ms'];
    title(  figtxt,                             ...
            'FontSize', 16                      );
    
    fStr = ['\\itTissue T1 time:\n\\rm'         ...
            '  T1_{n}   = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling IR:\n\\rm\\bf'  ...
            '  TI_{n}    = %4.i ms'             ];
    fVar = [ uint16(T1_n1), uint16(TI1_n) ];
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    pStr = ['\n'                                ...
            '*******************************\n' ...
            '***    IR tissue nulling:   ***\n' ...
            '*******************************\n' ...
            '*  Tissue T1       = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  Nulling TI      = %4.i ms  *\n' ...
            '*******************************\n' ];
    pVar = [ T1_n1, TI1_n ];

end % mainSIRNul

% 2) Calculate TI1+TI2 in a DIR sequence, nulling two specified T1 times (T1_n1 = T1_CSF by default)
function mainDIRNul()
    T1z = T1_n2;                                    % What T1 are we concerned with nulling?
    [ TI1_n, TI2_n ] = Find2TIn2T1(T1_n1,T1_n2);    % Inversion times that null 2 T1
    
%     S0DIR_CSF_debug = CalcS0D(TI1_n, TI2_n, T1_CSF) % Debug DIR S0 calc.
%     S0nIR_CSF_debug = CalcS0([TI1_n,TI2_n], [T1_CSF]) % Debug generic S0 calc.
    
    TInts = [ 0, (TI1_n-TI2_n), TI1_n, sTR ];       % DIR: Vector of all time points (inversion, readout and end)
    Mzt = PlotMagZ( TInts );                        % Plot Z magnetization over time, return end MagZ vector
    
    figtxt = ['DIR dual tissue nulling at ' B0tag ', for TR_{eff} = ' num2str(sTR_Eff) ' ms'];
    title(  figtxt,                             ...
            'FontSize', 16                      );
    
    fStr = ['\\itTissue T1 times:\n\\rm'        ...
            '  T1_{n1}  = %4.i ms\n'            ...
            '  T1_{n2}  = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling DIR:\n\\rm\\bf' ...
            '  TI_{1}    = %4.i ms\n'           ...
            '  TI_{2}    = %4.i ms'             ];
    fVar = [ T1_n1, T1_n2, TI1_n, TI2_n ];
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    pStr = ['\n'                                ...
            '*******************************\n' ...
            '***   DIR tissue nulling:   ***\n' ...
            '*******************************\n' ...
            '*  T1 tissue 1     = %4.i ms  *\n' ...
            '*  T1 tissue 2     = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  TI1             = %4.i ms  *\n' ...
            '*  TI2             = %4.i ms  *\n' ...
            '*  TI1 - TI2       = %4.i ms  *\n' ...
            '*******************************\n' ];
    pVar = [ T1_n1, T1_n2, TI1_n, TI2_n, (TI1_n-TI2_n) ];

end % mainDIRNul

% 3) Calculate & plot MRI S0 of WM vs MS/WML in a DIR sequence, and TI1+TI2 nulling their T1 weighting
function mainT1WNul()
    T1z = T1set(T1sel);                             % What T1 are we concerned with nulling?
    TI2 = 1:1:TI2_max;                              % Vector of TI2 values 1 ms apart (could've used linspace() here)
    LEN = length(TI2);                              % Length = Array size in the biggest dimension
    TI1 = zeros(LEN,1);                             % Vector of TI1 values for all TI2s; replace by TI1s(TI2s)?
    S0_DIR = zeros(LEN,2);                          % Array of two DIR signals@TE=0 for TI2s
    
%     TI1s(t) = CalcTI1(t,T1_n1);                     % TI1(TI2) calc.
    for i = 1:LEN
        TI1(i)      = CalcTI1D(TI2(i),T1_n1);       % TI1(TI2) calc.
        S0_DIR(i,1) = CalcS0D(TI1(i),TI2(i),T1_WM); % S0(TI2) calc.
        S0_DIR(i,2) = CalcS0D(TI1(i),TI2(i),T1_MS); % --"--
    end
    
%     S0_DIFF = abs(S0_DIR(:,1)-S0_DIR(:,2));         % Calculate signal difference between two T1s
%     [~,minind] = min(S0_DIFF);                      % Find the T1Wnull crossing graphically (as [minval,minind])
%     TI2_n = TI2(minind);
    TI2_n = uint16(FindT1WnTI2(T1_WM, T1_MS));      % Find [TI1,TI2] to null the DIR T1 weighting between two T1s
    TI1_n = uint16(CalcTI1D(TI2_n,T1_n1));          % TI1(TI2) calc. by function
%     TI2_n = uint16(TI2_n);                        % Make figure text nicer w/o e+ notation
%     TI1_n = uint16(TI1s( TI2_n ));                % TI1(TI2) calc. by symbolic expr. (slower?!)
    
    T1_n2 = uint16(FindT1nulled(TI1_n,TI2_n));      % Second tissue nulled by DIR. UInt for display.
    
    hold on
    plot(TI2,S0_DIR(:,1), 'r--')                    % Plot the signals calculated above
    plot(TI2,S0_DIR(:,2), 'b--')
%     plot(TI2,TI1);                                  % Debug: Plot the TI1 calculated for each TI2
%     fplot( t, TI1s(t) );                            % Debug: Plot the symbolic TI1(TI2) - WIP
%     plot(TI2,S0_DIFF);                              % Debug: Plot the WM/MS signal difference
    hold off
    
%     ax = gca;                                       % Get the handle of the current axes
    whitebg( 'white' );
%     whitebg( [ 0.70 1.00 0.70 ] );                  % Debug: See white plot/GUI elements
    figtxt =                                    ...
        ['T1-DIR (T1W nulling) at ' B0tag ', for TR_{eff} = ' num2str(sTR_Eff) ' ms'];
    title(  figtxt,                             ...
            'FontSize', 16                      );
    ylabel('Rel. MR signal (TE = 0 ms)', 'FontSize', 12 )
    xlabel('TI_2 inversion time (ms)',   'FontSize', 12 )
    xlm = [ 0 200*ceil( TI2_max/200 ) ]; xlim(xlm); % Round up to 200 on the x axis
    ylm = [ -1 1 ];                      ylim(ylm); % Rel. S0 is in [-1,1]
    
    lxlm = [ TI2_n TI2_n ];                         % S0 crossing x (=TI2) value
    lylm = [ ylm(1) S0_DIR(lxlm(1),1) ];            % Lower/upper y of line (=S0 curve)
    line( lxlm, lylm,                           ... % Vertical line @ crossing
        'Color', 'black',                       ... % (Tip: If declared after legend, the line appears in it)
        'DisplayName', 'TI_{2} nullT1W',        ...
        'LineStyle', ':'                        );
    
    Tis = [ 1 4 ];                                  % Make a legend for tissue 1(WM) & 4(MS)... [TODO: Use T1plt/T1txt?]
    NT1s = length(Tis); LegS = strings(1,NT1s);
    for i = 1:NT1s, j = Tis(i);                     % ...showing the tissue short names and T1 times
        LegS(i) = sprintf('%s (T1 = %i ms)', T1tag(j), T1set(j) );
    end
    legend( LegS(1), LegS(2),                   ... % Was: legend( 'WM', 'MS', ...
            'Location', 'NorthWest'             );  % 'Best' may conflict with UI and TextBox
    
%     fStr = ['\\itTissue T1 times:\n\\rm'        ...
%             '  T1_{WM} = %4.i ms\n'             ...
%             '  T1_{MS} = %4.i ms\n\n'           ...
    fStr = ['\\itT1W-0 DIR times:\n\\rm\\bf'    ... % Tex \bf\it\rm = bold/italic/normal; escape \ for sprintf!
            '  TI_{1}    = %4.i ms\n'           ...
            '  TI_{2}    = %4.i ms\n'           ...
            '  T1_{nt}  = %4.i ms'              ];
    fVar = [ TI1_n, TI2_n, T1_n2 ];                 % Was: T1_WM, T1_MS, etc
    figTextBox( fStr, fVar );                       % Display an info text box on the figure

    pStr = ['\n'                                ...
            '*******************************\n' ...
            '***    DIR T1-W nulling:    ***\n' ...
            '*******************************\n' ...
            '*  T1 nulled (CSF) = %4.i ms  *\n' ...
            '*  T1 also nulled  = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  TI1             = %4.i ms  *\n' ...
            '*  TI2             = %4.i ms  *\n' ...
            '*  TI1 - TI2       = %4.i ms  *\n' ...
            '*******************************\n' ];
    pVar = [ T1_n1, T1_n2, TI1_n, TI2_n, (TI1_n-TI2_n) ];
    
end % mainT1WNul

%% FUNCTIONS
function FltOrSym = FoS( InNr )                     % Cast as double for, e.g., exp() unless it's used as a sym
    if      isa( InNr, 'sym' ), FltOrSym = InNr;
    else,                       FltOrSym = double(InNr);
    end % if
end % fn

function S0_IR = CalcS0( Ti, T1 )                   % Calculate generic IR relative S0 (def: S0 = S/(E2*M0) at TE=0)
    Ti = FoS(Ti); T1 = FoS(T1);                     % Cast input as float unless it's a sym
    NInv = length( Ti ); S0_IR = ones(1,length(T1));
    for i = 0:NInv-1
        S0_IR = S0_IR - (-IE^i)*IE2*exp(-Ti(NInv-i)./T1); % Note: TIv() is defined as time(s) to readout!
    end
    S0_IR = S0_IR - (-IE^NInv)*exp(-sTR_Eff./T1);
end % fn

function S0_IR = CalcS0S( Ti, T1 )                  % Calculate single IR relative S0 (def: S0 = S/(E2*M0) at TE=0)
    Ti = FoS(Ti); T1 = FoS(T1);                     % Cast input as float unless it's a sym
    S0_IR = 1 - IE2*exp(-Ti/T1) + IE*exp(-sTR_Eff/T1);
end % fn

function S0_DIR = CalcS0D( Ti1, Ti2, T1 )           % Calculate DIR relative S0 (def: S0 = S/(E2*M0) at TE=0)
    Ti1 = FoS(Ti1); Ti2 = FoS(Ti2); T1 = FoS(T1);   % Cast input as float unless it's a sym
    S0_DIR = 1 - IE2*exp(-Ti2/T1) + IE*IE2*exp(-Ti1/T1) - IE^2*exp(-sTR_Eff/T1); % (New def.: TI1_new=TI1_old+TI2_old)
end % fn

function TI1_DIR = CalcTI1D( Ti2, T1 )              % Calc. TI1 (new TI1, the old one was - TI2) nulling T1 given TI2
    Ti2 = FoS(Ti2); T1 = FoS(T1);                   % Cast input as float unless it's a sym
    TI1_DIR = log(IE*IE2/(IE2*exp(-Ti2/T1) + IE^2*exp(-sTR_Eff/T1) - 1))*T1; % This may be a sym, so don't cast to int
end % fn

function T1WnTI2 = FindT1WnTI2( T1w1, T1w2 )        % Find [TI1,TI2] to null the DIR T1 weighting between two T1s
    T1w1 = FoS(T1w1); T1w2 = FoS(T1w2);             % Cast input as float unless it's a sym
    T1WnTI2 = vpasolve((CalcS0D(CalcTI1D(t,T1_n1),t,T1w1) ...   % Solve S0DIR == 0 for T1
                      - CalcS0D(CalcTI1D(t,T1_n1),t,T1w2)), t); % Return TI2 (then calculate TI1(TI2) later)
end % fn

function T1n_DIR = FindT1nulled( Ti1o, Ti2o )       % Find T1 time nulled by DIR
    Ti1o = FoS(Ti1o); Ti2o = FoS(Ti2o);             % Cast input as float unless it's a sym
    T1n_DIR = vpasolve(CalcS0D(Ti1o,Ti2o,t), t);    % Solve S0DIR == 0 for T1
end % fn

function [ Ti1n, Ti2n ] = Find2TIn2T1( T1n1, T1n2 ) % Find TI1/TI2 that null T1n1/2
    T1n1 = FoS(T1n1); T1n2 = FoS(T1n2);             % Cast input as float unless it's a sym
    Ti2n = vpasolve(CalcS0D(CalcTI1D(t,T1n1),t,T1n2), t); % (Used VPAsolve for numeric output but still got syms?)
    Ti1n = uint16(CalcTI1D(Ti2n,T1n1));             % Make value int for figure display; round() doesn't change the sym!
    Ti2n = uint16( Ti2n );                          % --"--
end % fn

% % Calculate the relative Z magnetization Mz after inversions at time points [TI]
% function MagZ = CalcMagZ( TI, DUR, T1v )                 % ( TI, DUR, T1set )
%     NInv = length(TI);                              % # of inversion time points
%     
%     dTI = CalcDelTI( [ TI, DUR ] );
% %     NT1s = length(T1v);                             % # of T1 times for which to calculate Mz
%     MagZ = zeros(1,length(T1v));                    % End Z magnetization for each T1
%     for i = 1:NInv+1
%         MagZ = 1 - (1 + IE*MagZ).*exp(-dTI(i)./T1v);
%     end
% end % fn

% % Calculate a vector of dTI() from TIT[] times (absolute times)
% function dTI = CalcDelTI( TIT )
%     NIT = length(TIT);                              % # of time points (TI# and DUR)
%     dTI = zeros(1,NIT);                             % Vector of the differences between TIs
%     dTI(1) = TIT(1);                                % Time point for the first inversion
%     for i = 1:NIT-1
%          dTI(i+1) = TIT(i+1) - TIT(i);              % dT(2) = TI(2) - TI(1); dT(3) = DUR - TI(2)
%     end
% end % fn

function SetMode( mode )                            % Sets the run mode and resets the Nulling selection (T1sel)
    RunMode = mode;                                 % Make the mode setting global
    
    T1_n1 = T1_CSF;                                 % Default for all modes: Null out CSF as the longest T1
    switch RunMode
        case 1                                      % 1IR default: Null out CSF (=FLAIR)
            SetT1n( 3 );
        case 2                                      % DIR default: Null out WM as tissue 2
            SetT1n( 1 )
        case 3                                      % T1W-DIR default: (Null out CSF as the longest T1)
    end % switch
end % fn

function SetT1n( tissue )                           % Sets the Nulling selection and T1_n#
    T1sel = tissue;
    switch RunMode
        case 1                                      % 1IR (only one tissue nulled)
            T1_n1 = T1set(T1sel);
        otherwise                                   % DIR (set second T1 to null; the first is CSF by default)
            if T1set(T1sel) >= T1_n1, T1sel=1; end  % Don't null T1_n1/CSF twice; use WM as tissue 2 instead then.
            T1_n2 = T1set(T1sel);
    end % switch
    SetT1plt()                                      % Reset which T1 times to plot (based on RunMode)
end % fn

function SetT1plt()                                 % Select which T1 times to plot [TODO: When to call this?]
    switch RunMode                                  
        case 1                                      % 1IR: Plot the signal for 3 brain tissues (WM, GM, CSF)
            i = 1:3;
        case 2                                      % DIR: Plot the signal for 3 brain tissues (WM, GM, CSF)
            i = 1:3;
        case 3                                      % T1W-DIR: Plot the signal for WM and WML
            i = [ 1 4 ];                            % TODO: This is not implemented yet
    end % switch
    T1plt = T1set(i); T1txt = T1tag(i);             % Update the T1 times to plot and their label texts
end % fn

function setRelaxTs( B0opt )
    B0sel = B0opt;                                  % Ensures right selection in the UI
    B0str = [   "1.5 T (Alsop)", "3.0 T (Alsop)", ...                   % Array of B0/sources for relax. times
                "3.0 T (GE)", "1.5 T (MRIQue)", "3.0 T (Lalande)" ];    %   (used in the UI)
    B0tag = char(B0str(B0sel)); B0tag = B0tag(1:5); % Char array of only the field strength, for figure titles
    T1str = [   " T1 WM     ", " T1 GM     ",   ...                     % Array of tissue names
                " T1 CSF    ", " T1 WML/MS ", " T1 Fat    " ];          %   (for formatted output)
    T1 = split(T1str); T1 = char(T1(:,:,3));        % Array of only the tissue type (3rd word) for legends,...
    T1tag = string(T1(:,1:3,:));                    % ...trimmed to length 3 (e.g., "CSF" or "WML")
    
    switch B0opt
        case 1                  % 1.5 T (Alsop); Madhuranthakam et al, Mag Res Med 2012 67(1):81-88
            T1_WM  =  650.0;            % 1.5 T Alsop
            T1_GM  = 1300.0;            % 1.5 T --"--
            T1_CSF = 4200.0;            % 1.5 T --"--
            T1_MS  = T1_GM;             % Estimate T1_WML ? T1_GM
            T2_WM  =   70.0;            % 1.5 T Alsop
            %T2_GM  =    ???
            T2_CSF = 2000.0;            % 1.5 T Alsop
            T2_MS  =  100.0;            % 1.5 T --"--
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            T1str(4) = "(T1 MS/WML)";   % This value was derived from elsewhere
            T1str(5) = "(T1 Fat)   ";   % This value was derived from elsewhere
        case 2                  % 3.0 T (Alsop)
            T1_WM  =  750.0;            % 3.0 T Alsop
            T1_GM  = 1400.0;            % 3.0 T litteratur
            T1_CSF = 4200.0;            % 3.0 T litteratur
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean)
            T1str(2) = "(T1 GM)    ";   % This value was derived from elsewhere
            T1str(5) = "(T1 Fat)   ";   % This value was derived from elsewhere
        case 3                  % 3.0 T (GE)
            T1_WM  =  825.0;            % 3.0 T litteratur (GE har 825 ms, andre 850 ms?)
            T1_GM  = 1400.0;            % 3.0 T --"--
            T1_CSF = 4200.0;            % 3.0 T --"--
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean)
            T1str(4) = "(T1 MS/WML)";   % This value was derived from elsewhere
            T1str(5) = "(T1 Fat)   ";   % This value was derived from elsewhere
        case 4                  % 1.5 T (MRIQ); MRIquestions.com and other online resources
            T1_WM  =  580.0;            % 1.5 T MRIQ
            T1_GM  =  940.0;            % 1.5 T --"--
            T1_CSF = 3600.0;            % 1.5 T --"--
            T1_MS  = T1_GM;             % Estimate T1_WML ? T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQ
            T1str(4) = "(T1 MS/WML)";   % This value was derived from elsewhere
        case 5                  % 3.0 T (Lalande); Bojorquez et al, MRI 35(2017) 69–80
            T1_WM  =  842.0;            % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR)
            T1_GM  = 1425.0;            % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR)
            T1_CSF = 4300.0;            % 3.0 T Lalande (trimmed mean of 4(6) reports)
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean; 405.0 in Barral 2010 w/ SE-IR)
            T1str(4) = "(T1 MS/WML)";   % This value was derived from elsewhere
        otherwise
            error('ERROR: Undefined B0/source!');
    end % switch B0opt
    T1set = [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ];
    
    T1_n1 = T1_CSF;                     % Default for all modes: Null out CSF as the longest T1
    switch RunMode                      % When a scheme is selected, adjust the nulling selection accordingly.
        case 1                          % 1IR T1 nulling
            SetT1n( T1sel );
        case 2                          % DIR T1 nulling
            SetT1n( T1sel );            % (Check and set T1_n2)
        case 3                          % T1W nulling
    end % switch
end % setRelaxTs

%% PLOT

% Plot the magnetization Mz of several T1s under multiple inversion time points (NB: Absolute TI times are used here!)
function MagZ = PlotMagZ( TIabs )                   % ( TIabs, T1ts, T1cs ); TIabs = [ dTI1 dTI1+dTI2 Tacq DUR ]
    NIT = length(TIabs)-2;                          % # of inversion time points, without readout time and duration
    TIn = TIabs(1:NIT);                             % The actual inversion time points
    DUR = TIabs(NIT+2);                             % Total duration for the plot
    TIR = TIabs(NIT+1); xRo = double([ TIR TIR ]);  % Readout time (=TI)
%     sZRo = 0.02*sT_Ro;                              % NOTE: sZRo is irrelevant as spin lock dictates readout from TIR?
%     TRo = xRo + [ -sZRo, (sT_Ro-sZRo) ];            % Readout time range. Linear(sZRo=0.5) or center-out (sZRo small)?
    TRo = xRo + [ 0, sT_Ro ];                       % Readout time range. Due to spin lock, this should start at TIR.
    dT = 1; Tp = (0:dT:DUR-1);                      % Sampling interval dT (ms, Int), Tp is the time point vector
    
    NT1s = length( T1plt );                         % # of T1 times for which to calculate Mz
    if ~ismember( T1z, T1plt )                      % Add T1z if not in the preselected T1 ensemble
        NT1s  = NT1s + 1;
        T1plt = [ T1plt, T1z ];
        T1txt = [ T1txt, sprintf('T1#%s', num2str(NT1s)) ];
    end % if
    Mz = ones( NT1s, DUR );                         % MagZ for all tissues and time points (was zeros(NT1s+1...)
    
    % WIP: Starting Z magnetization calculation.
%     Mz0 = ones(1,length(T1plt));                    % Test/debug: Start with relaxed signal (TR >> T1)
%     Mz0 = -IE*Mz0;                                  % Test/debug: Start with inversion
%     Mz0 = 0*Mz0;                                    % Test/debug: Start with saturation
%     Mz0 = -CalcS0( [ sTR_Eff, 0 ], T1plt )          % --"-- CalcS0 uses TI differences.
%     Mz0 = CalcMagZ( [ 0 ], sTR_Eff, T1plt )         % --"-- [1 - exp(-sTR_Eff./T1plt)]
    % TODO: Use derivation instead of iteration to arrive at SS? Maybe not so critical.
    its=8; MzIt = ones(NT1s,its);                   % Iterate to get to a steady state
    for i = 1:its                                   % Iteration loop (Note: Needs only one TR if Mz(TRo)=0.)
        Mz(:,1) = Mz(:,DUR);                            % The starting magnetization (before first inv.) is MagZ(TR)
        for j = 2:DUR                                   % Step through all time points
            Ts = (j-2)*dT;                              % Was (j-1)*dT, but I want the time point 0 to be defined
            Mz1 = Mz(:,j-1);                            % MagZ at the previous time point
            if ( any( Ts == TIn(:) ) )                  % Inversion time, so invert the Mz (accounting for Inv.Eff.)
                Mz(:,j) = -IE*Mz1;
            elseif ( Ts == TRo(1) )                     % Readout time
                MzTI    = Mz1;                          % MagZ at readout
                SRo     = sin(FA*pi/180)*Mz1;           % Signal strength at readout
                Mz(:,j) = cos(FA*pi/180)*Mz1;
                MzRo    = MzTI - Mz(:,j);
            elseif ( Ts > TRo(1) ) && ( Ts < TRo(2) )   % Readout with CPMG train, so don't evolve signal
                Mz(:,j) = Mz1;                          % TODO: Evolve slowly with adjusted T1s instead?
            elseif ( Ts == TRo(2) )                     % Readout end
                Mz(:,j) = Mz1 + REW*MzRo;               % Rewinder reclaims the signal lost in excitation? Or does it?
            else                                        % Signal evolvement ( if(t~=TI(:)) )
                for k = 1:NT1s                          % Can't assign all tissues at once with exp, so use a loop
    %                 Mz(k,j) = Mz(k,j-1)*exp(-dT/T1v(k)) + (1-exp(-dT/T1v(k)));
                    Mz(k,j) = 1 - (1 - Mz(k,j-1))*exp(-dT/T1plt(k));
                end % for k
            end % if Ts
        end % for j
        MzIt(:,i) = Mz(:,DUR);                      % The end Z magnetization for all T1s, for this iteration
    end % for i
    MagZ = Mz(:,DUR)';                              % Return the end Z magnetization for all T1s as a vector
    
    Mz = [ Mz; zeros(1,DUR) ];                      % Add a zero line to the plot
    plot( Tp, Mz )                                  % Plot the (steady state) Z magnetization
    
    xlabel('t (ms)', 'FontSize', 12 )
    ylabel('Rel. S0', 'FontSize', 12 )              % ''
    xlm = [  0  DUR ]; xlim(xlm);
    ylm = [ -1   1  ]; ylim(ylm);
%     set( gca, 'XTick', [] );                        % Remove x axis ticks
%     set( gca, 'YTick', [] );                        % Remove y axis ticks

    line( TRo, [ 0 0 ],                         ... % Horizontal line representing readout
        'Color', [ 0.3 0.7 0.3 ],               ... % 'black'
        'LineWidth', 6,                         ...
        'DisplayName', 'Echo train',            ...
        'LineStyle', '-'                        );
    
    line( xRo, ylm,                             ... % Vertical line at TEeff (=TI)
        'Color', 'black',                       ... % (Tip: If declared after legend, the line appears in it)
        'DisplayName', 'TI_{1}',                ... % 'TE_{eff}'
        'LineStyle', ':'                        );
    
    LegS = strings(1,NT1s);                             % Figure legends
    for i = 1:NT1s
        LegS(i) = sprintf('%4s (T1: %4i ms)', T1txt(i), T1plt(i) );  %legend('WM', 'GM', 'CSF'...
    end % for i
    legend( LegS, 'Location', 'NorthWest');
end % PlotMagZ

function figTextBox( fStr, fVar )
    figtxt = sprintf( fStr, fVar );                 % Formatted string "print"
    xlm = xlim; xt = double(xlm(2));                % Max value of axis; (x,y) is a plot point!
    ylm = ylim; yt = double(ylm(2));                % --"--
    text( 0.794*xt, -0.620*yt, figtxt,          ... % Note: annotation('textbox' etc doesn't scale with text
        'EdgeColor', 'black',                   ...
        'Margin', 4,                            ...
        'LineWidth', 1,                         ...
        'LineStyle', '-'                        );
end % fn

function saveFigAsImg()
    imgName = ['TReff_' int2str(sTR_Eff) '_TI1_' int2str(TI1_n) '_TI2_'  int2str(TI2_n) ];
    if exist('Figures','dir')~=7, mkdir('Figures'), end
    cd Figures
    % TODO: What's the difference between print to image and the imwrite command?
    print('-djpeg', '-r300', imgName )
    cd '..'
end % fn

%% USER INTERFACE

function createUI()                                 % Display a button to show/hide the UI controls
    uiX = 74.0*fSz; uiY = 50.0*fSz;                 % Default start position for masterUI controls
    eW = 50       ; eH = 40;                        % Button width/height for masterUI controls
    
    uicontrol( 'Style', 'pushbutton',           ... % uiMode
        'String', 'Mode',                       ... % 'togglebutton'
        'FontWeight', 'bold',                   ...
        'ToolTipString', 'Switch mode',         ...
        'Position', [ uiX uiY-00 eW eH ],       ...
        'Callback', @cModes_callback            );  % UI to set run mode
    
    uicontrol( 'Style', 'pushbutton',           ... % uiCtrl
        'String', 'Ctrl',                       ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', 'Show/hide UI',        ...
        'Position', [ uiX uiY-60 eW eH ],       ...
        'Callback', @masterUI_callback          );  % UI to show/hide the UI control panel
    
    createUIPanel();
    switch UIon
        case true ; uiBox.Visible = 'on';
        case false; uiBox.Visible = 'off';
    end % switch UIon
    switch RunMode
        case 1                                      % 1IR T1 nulling
        uiTn1T.String  = 'T1_n:                         ms'; % There is only one T1 to null in 1IR
        uiSTnT.Position = uiSTnT.Position + [ 0 30 0 0 ];   % Move the tissue selector one line up
        uiSTnS.Position = uiSTnS.Position + [ 0 30 0 0 ];   % --"--
%         uiTn2S.Visible = 'off';                         % Hide the tissue 2 edit text/box (ui###T/S, if created)
        case 2                                      % DIR T1 nulling
%         uiB0S.ForegroundColor = myGray;                 % Gray out the B0 select text/box (ui###T/S, if created)
        case 3                                      % T1W nulling
        uiTn2S.ForegroundColor = myGray;                % Gray out the tissue 2 edit box (it's set automatically)
        uiFAsS.ForegroundColor = myGray;                % Gray out the FA edit box (FA is irrelevant in this mode)
    end % switch RunMode
    if T1z ~= T1set(T1sel)
        uiSTnS.ForegroundColor = myGray;        % Gray out the preselection box if entering values manually
    end % if T1z
end % createUI

function createUIPanel()
    % Make UI control buttons for saving images instead of the variable, like with printVars?
    %      Probably not necessary, since you can save the image from the File menu.
    
    boxClr = figClr;                                % UI panel background color (same as fig background)
    uiX = 58.3*fSz; uiY = 54.5*fSz;                 % Default start position for UI controls
    mX = 10; mY =  6;                               % x/y margins in UI panel
    eN = 10; eW = 120; eH = 30;                     % Number of UI elements; default element width/height
    uiH = eN*eH + 2*mY ; uiT = uiH - mY +  3;       % UI panel height (pixels); top position for elements
    
    uiBox = uipanel( myFig,                     ... % uiBox
        'BackgroundColor', boxClr,              ... %  [ 0.8 0.8 0.4 ]
        'BorderType', 'etchedin',               ... % default 'etchedin'
        'BorderWidth', 1.0,                     ... % default 1
        'Clipping', 'off',                      ... % default 'on'; 'off' is good for debugging
        'Units', 'pixels',                      ... % With normalized units (def.), got Pos. [ .726 .420 .172 .500 ]
        'Position', [ uiX uiY-uiH eW+2*mX uiH ] );  % UI panel container for the controls (pixel units)
%         'Title', 'Controls',                    ... % If specifying a title, it'll show up on the panel
%         'FontSize', 10,                         ...
%     uistack( uiBox, 'top' );                        % Bring the panel to the top of the graphic stack (not necessary)
    
    for i = 1:eN; eY = uiT - i*eH;                    % Generate each UI element on a new row in the panel
        switch i
        case 1
        uicontrol( 'Parent', uiBox,                 ... % uipVar
            'Position', [ mX eY+0 eW 24 ],          ...
            'Style', 'pushbutton',                  ...
            'String', 'Print results',              ...
            'Callback', @printVar_callback          );  % UI to print result variables in the command window
        
        case  2
        uicontrol( 'Parent', uiBox,                 ... % uiTR[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Repetition time',     ...
            'String', 'TR:                             ms');
        uicontrol( 'Parent', uiBox,                 ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'String',   sTR,                        ...
            'Value',    sTR,                        ...
            'Callback', @tr_set_callback            );  % UI to set TR

        case  3
        uicontrol( 'Parent', uiBox,                 ... % uiTRo[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Readout time',        ...
            'String', 'Tro:                            ms');
        uicontrol( 'Parent', uiBox,                 ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'String',   sT_Ro,                      ...
            'Value',    sT_Ro,                      ...
            'Callback', @tro_set_callback           );  % UI to set T_Ro

        case  4
        uicontrol( 'Parent', uiBox,                 ... % uiFAs[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Readout flip angle',  ...
            'String', 'FAex:                         °'); % Note: UIControl text can't support Tex/LaTeX
        uiFAsS = uicontrol( 'Parent', uiBox,        ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 10000, 'Max', 10000,             ...
            'String',   FA,                         ...
            'Value',    FA,                         ...
            'Callback', @fa_set_callback            );  % UI to set FA (in degrees)

        case  5
        uicontrol( 'Parent', uiBox,                 ... % uiIE[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Inversion efficiency',...
            'String', 'Inv.eff.:                      %');
        uiIES = uicontrol( 'Parent', uiBox,         ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'String',   IEf,                        ...
            'Value',    IEf,                        ...
            'ToolTipString', '-cos(FA_inv)',        ...
            'Callback', @ie_set_callback            );  % UI to set inversion efficiency

        case  6
        uiB0T = uicontrol( 'Parent', uiBox,         ... % uiB0[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Ref. values for T1s', ...
            'String', 'B0:'                         );
        uiB0S = uicontrol( 'Parent', uiBox,         ...
            'Position', [ mX+22 eY+0 eW-22 23 ],    ... % Note: Height setting won't change the actual height here
            'Style', 'popup',                       ...
            'String',   B0str,                      ...
            'Value',    B0sel,                      ...
            'Callback', @b0_sel_callback            );  % UI to set B0 (relaxation times)

        case  7
        uicontrol( 'Parent', uiBox,                 ... % uipT1s
            'Position', [ mX eY+0 eW 24 ],          ...
            'Style', 'pushbutton',                  ...
            'String', 'Print T1 times',             ...
            'Callback', @printT1s_callback          );  % UI to print T1s in the command window

        case  8
        uiTn1T = uicontrol( 'Parent', uiBox,        ... % uiTn1[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'T1 time to null',     ...
            'String', 'T1_n1:                       ms'); % Note: UIControl text can't support Tex/LaTeX
        uicontrol( 'Parent', uiBox,                 ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 10000, 'Max', 10000,             ...
            'String',   T1_n1,                      ...
            'Value',    T1_n1,                      ...
            'Callback', @t1n1_set_callback          );  % UI to set T1_n1 (longer T1, typically CSF)

        case  9
            if ismember( RunMode, [ 2 3 ] )                 % DIR T1 nulling and T1-DIR null two T1
            uiTn2T = uicontrol( 'Parent', uiBox,        ... % uiTn2[T|S]
                'Position', [ mX eY+1 eW 16 ],          ...
                'Style', 'text',                        ...
                'BackgroundColor', boxClr,              ...
                'HorizontalAlignment', 'Left',          ...
                'ToolTipString', 'Shorter T1 to null',  ...
                'String', 'T1_n2:                       ms');
            uiTn2S = uicontrol( 'Parent', uiBox,        ...
                'Position', [ mX+40 eY+0 60 20 ],       ...
                'Style', 'edit',                        ...
                'Min', 10000, 'Max', 10000,             ...
                'String',   T1_n2,                      ...
                'Value',    T1_n2,                      ...
                'Callback', @t1n2_set_callback          );  % UI to set T1_n2 (shorter T1; GM/WM/etc)
            end % if

        case 10
            if ismember( RunMode, [ 1 2 ] )                 % 1IR and DIR T1 nulling use tissue selector
            uiSTnT = uicontrol( 'Parent', uiBox,        ... % uiSTn[T|S]
                'Position', [ mX eY+1 eW 16 ],          ...
                'Style', 'text',                        ...
                'BackgroundColor', boxClr,              ...
                'HorizontalAlignment', 'Left',          ...
                'ToolTipString', 'Ref. tissue to null', ...
                'String', 'Nulling:'                    );
            uiSTnS = uicontrol( 'Parent', uiBox,        ...
                'Position', [ mX+40 eY+0 eW-40 22 ],    ...
                'Style', 'popup',                       ...
                'String',   T1str,                      ...
                'Value',    T1sel,                      ...
                'Callback', @t1_sel_callback            );  % UI to set T1_n(2) by tissue
%                 'String', {'WM','GM','CSF','Fat'},      ...
            end % if
        end % switch i
    end % for i
end % createUIPanel

function masterUI_callback(~,~)                     % UI that shows/hides the UI control panel
    UIon = ~UIon;
    switch UIon
        case true ; uiBox.Visible = 'on';
        case false; uiBox.Visible = 'off';
    end
end % fn

function cModes_callback(~,~)                       % Switch between calculation modes
    switch RunMode
        case 1                                      % 1IR T1 nulling    ->|
            SetMode(2);
        case 2                                      % DIR T1 nulling    ->|
            SetMode(3);
        case 3                                      % T1W nulling       <<-
            SetMode(1);
    end                                             % Note: The figure title shows the active mode
    main();
end % fn

function b0_sel_callback(src,~)                     % UI to set B0 (rel. times)
    setRelaxTs( src.Value );                        % Was get(src,'Value') pre-R2014b
    main();
end % fn

function tr_set_callback(src,~)                     % UI to set TR
    sTR = str2double( src.String );
    Lim = sT_Ro + 5; if sTR < Lim, sTR = Lim; end   % TR > T_Ro
    main();
end % fn

function tro_set_callback(src,~)                    % UI to set T_Ro
    sT_Ro = str2double( src.String );
    Lim = sTR - 5; if sT_Ro > Lim, sT_Ro = Lim; end % TR > T_Ro
    main();
end % fn

function t1n1_set_callback(src,~)                   % UI to set T1_n1 (longer T1)
    T1_n1 = str2double( src.String );
    if RunMode > 1                                  % Using multiple inversions...
        Lim = T1_n2 + 1;
        if T1_n1 < Lim, T1_n1 = Lim; end            % T1_n1 is the longest T1 time to null
    end
    main();
end % fn

function t1n2_set_callback(src,~)                   % UI to set T1_n2 (shorter T1)
    T1_n2 = str2double( src.String );
    Lim = T1_n1 - 1;
    if T1_n2 > Lim, T1_n2 = Lim; end                % T1_n1 is the longest T1 time to null
    main();
end % fn

function t1_sel_callback(src,~)                     % UI to select a tissue T1 to null
    SetT1n( src.Value );
    main();
end % fn

function fa_set_callback(src,~)                     % UI to select the flip angle (in degrees)
    FA = str2double( src.String );
    main();
end % fn

function ie_set_callback(src,~)                     % UI to select the inversion efficiency (in percent)
    IEf = str2double( src.String );
    main();
end % fn

function printVar_callback(~,~)                     % UI that prints variables to the command window
    fprintf( pStr, pVar );
    if debugInfo && ismember( RunMode, [ 1 2 ] )    % Debug info for plot
        len = length(MzIt);
        fprintf('\nMagZ(TR) over the first %i repetitions (starting at 1):\n', len);
        disp( MzIt );                               % Show the iterative Z magnetization over the first TRs
%         disp( MzIt(:,2:its+1) );
        fStr = [ 'MagZ and Signal(0) by T1 at the'  ...
            ' last readout with flip angle %3.1f°'  ...
            ' [=%1.3f*MagZ(TI)]:\n'                  ];
        fprintf(fStr, round(FA,1,'significant'), sin(FA*pi/180));
        fprintf("\t "); fprintf("%-10s", T1txt); fprintf("\n");
        fStr = [ MzTI'; SRo' ]; disp ( fStr );
    end
end % fn

function printT1s_callback(~,~)                     % UI that prints variables to the command window
    fStr = ['\n'                                ...
            '*******************************\n' ...
            '***    Relaxation times     ***\n' ...
            '***    for %-15s',       '  ***\n' ...
            '*******************************\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ... % Printing T1s directly, %4.i was used:
            '*******************************\n' ];  % fprintf( fStr, B0str(B0sel), T1set );
%     fVar = [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ];
    T1info = [ T1str; num2cell(T1set) ];            % Use this trick to weave strings and numbers into a string array
    fprintf( fStr, B0str(B0sel), T1info{:} );
end % fn

end % script fn