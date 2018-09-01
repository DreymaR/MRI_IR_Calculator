function IR_Calculator()
%% (D)IR MRI TI/T1 calculation
%   - In mode 1 & 2, find TI time(s) that null given T1(s) in DIR & 1IR (FLAIR/STIR/etc).
%   - In mode 3, calculate/plot relative signal strength for WM & WML as a function of TI2 in a DIR sequence,
%       - find TI2 and TI1(TI2) such that there is no T1 weighting between WM & WML while nulling CSF,
%       - and find a T1 that's nulled at the resulting TIs, corresponding to a fictive tissue.
%   - Also plot the T2 decay at readout (in a separate file), to determine the optimal TE.
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
%   - Rewrite with handles structure ( use guihandles(), guidata() etc ) like GUIDE; check GUIDE first.
%   - Use less globals! Maybe also a struct/cell array for a set of properties/settings?
%   - Include "T1W WM/WML" in DIR nulling options!?
%       - Instead of the Mode button, make it show "IR"/"DIR", depending.
%       - Then, with DIR, add another button for "T1W nulling"
%       - When leaving T1W nulling mode, keep the T1_n2 so DIR shows it.
%   - Include AtleB's ASL BS TI optimalization? Use a press button for it, as it takes some time to run.
%   - Implement incomplete spin lock during readout? CPMG inversion pulses may be 120-150° instead of 180°.
%       - So readout CPMG inversion efficiency is then -cos(FA)=50-87%
%       - Merely simulating incomplete inversions doesn't lead to a change in behaviour, as spin lock is still strong.
%       - Likely we'll need a more thorough Bloch simulation to account for all the spin components.
%       - Or is it true that the T1 recovery is completely negligible during a FSE train, even with low FAs?
%   - Simulate "T2 prep"? In which shorter T2 times are dephased so their signal becomes nearly saturation recovery.
%       - Replace the unused(?) excitation FA UI with a T2 prep (ms) one?!?
%       - Allow T2 entry? Crowds the UI. Alternatively, only emulate T2prep if times aren't set manually.
%       - Could force T2prep = 0 when setting times manually?
%   - Make the T2 plot a subplot instead of a new figure!? Send its script our figure handle.
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
%   - Solve Mz(0) in TIR/DIR by formula [1 - E_(TR-TI1)] instead of iteration? No, it's good to simulate it as a check.
%   - Implement excitation 
%       - Without rewinding? What is the case in our sequences?
%       - Add a GUI for effective FA (90° for full excitation, 0° if refocused). Maybe also add rewinder?
%   - Implement incomplete inversion using the Inversion Efficiency = -cos(FA_inv)
%   - Added a generic S0 formula CalcS0() instead of CalcS0S() and CalcS0D().
%   - Made a T2-TE figure that can be updated from the GUI panel for the current signal/plot settings
%   - Correct CNR for sB0

%% INIT
% saveImg = false;                                    % Save figure as image in Figures dir? (Can use the menu instead)
debugInfo  = true; %false;                          % Show extra info in the command window output

global Ps UIs Tis                                       % Program/UI/tissue data structs
global pInit UIsBox                                     % Some globals that didn't thrive in the structs
global sTR sETD sTR_ef sIEf IE IE2 sFA sETL sESP sIET sRew sT2p    % System/sequence variables
global TisSel TsLeg TsPlt                               % Tissue labels/tags/legends
global sB0 B0sel B0str B0tag           oldB oldS newS   % B0 settings/labels (System B0 and reference values)
global Tset Ttag T1z                                    % T1/T2 tissue relaxation times and tags (based on "B0")
% global T1_WM T1_GM T1_CSF T1_MS T1_Fat T2_WM T2_GM T2_CSF T2_MS T2_Fat
global TI2_max T1_n1 T1_n2 TI1_n TI2_n SRo MzTI MzIt    % TI calculation, T1 tissue nulling, magnetization

% sSz = get(0,'ScreenSize');                          % Get the screen size (actually, in R2015b+ pixels are 1/96")
fSz = 10.80; % sSz(4)/100.00;                       % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])    % Default fig. position & size
wh = [ 1 1 1 ]; figClr = 0.94*wh; myGray = 0.50*wh; % Colors of fig. background and grayed-out UI elements
% fig1 = get(groot, 'CurrentFigure');                 % Get handle of the current fig. ('gcf' creates one)
figName = 'MR-IR inversion time calculator';
Ps.Fig1 = findobj( 'Type', 'Figure',            ... % Find this script's figure (possibly among others)
        'Name', figName                         );
if isempty(Ps.Fig1)                                 % Iff my figure does not yet exist... (was ~ishghandle)
    Ps.Fig1 = figure(                           ... % ...make a new figure!
        'Name', figName,                        ...
        'Color', figClr,                        ...
        'KeyPressFcn', @keyPress_callback,      ...
        'Position', [ 90 25 80 60 ]*fSz         );  % At pos. (90,25)% of ScreenSize, etc.
%         'NumberTitle', 'off',                   ... % Removes 'Figure #: ' from the window title
else
    figure(Ps.Fig1);                                % If my figure exists, make it the CurrentFigure
end % if

t = sym('t');                                       % Set the symbolic explicitly as it's used in nested functions
% TI2s = sym('TI2s'); TI1s(TI2s) = TI2s;              % Temporary assignment globalizes the sym expression
% DIR_S0(t) = ( 1 - 2*exp(-TI2/t) + 2*exp(-TI1/t) - exp(-TReff/t) ); DIR S0(T1) as a sym expression
% TI1s(TI2s) = T1*log(2/(2*exp(-TI2s/T1) + exp(-sTR_ef/T1) - 1)); % TI1(TI2) nulling T1_n1, using syms
% % Note: Do not replace CalcTI1(TI2,T1) with TI1s(t) above: It's much slower.

% Parameter presets for the first run (most can be changed in the UI)
if isempty( pInit )
    pInit = true;                                   % The program has been initialized

sTR     = 6000.0    ;                               % Repetition time (ms); as float (if isempty(sTR); sTR = 6000.0; end)
sETD    =  780.0    ;                               % Readout time (ms) ETL * ESP (GE: Check CVs after Download)
sIEf    =   95      ;                               % Inversion efficiency in percent
sT2p    =    0.0    ;                               % T2 preparation time (ms)
sFA     =   90      ;   %MEX = cos(sFA*pi/180);     % MEX is MagZ after excitation with angle FA
sRew    =    0.00   ;                               % Rewinder after readout, regaining some of the pre-readout Mz?
sESP = 3.6; sIET = 0.5; % sETL = 122;               % WIP: Echo Train Length, Echo Spacing and Inv.Eff. for TSE readout
                                                    % WIP: Just assuming incomplete inversion doesn't lead to a T1rho!

UIs.Show    = true  ;                               % Show the UI control panel
Ps.Mode     =    1  ;                               % Start in mode 1, then keep the run mode over reruns
TisSel      =    3  ;                               % Start with a value (CSF), then keep the selection over reruns
setRelaxTimes( 3 )  ;                               % Select 3T setup, with relaxation times from literature (GE sett.)
SetMode( Ps.Mode )  ;                               % Note: Need to set Ps.Mode -> TisSel -> RelaxTimes -> SetMode.
% TODO: Want instead to have if isempty(Ps.Mode)  ; setMode( 1 );
end % if pInit

%% MAIN
main();
function main()
    sTR_ef = sTR - sETD - sT2p;                     % Recalculate based on current parameter settings
    IE = sIEf/100; % IE2 = 1 + IE;                  % MagZ to recover after incomplete inv. (was log(2/...) etc)
    sETL = ceil(sETD/sESP);                         % Echo train length (# pulses) for TSE readout in plot
    
    clf                                             % Clear the figure so we can change it
    switch Ps.Mode
        case 1                                      % Mode 1 (1IR tissue nulling)
            mainSIRNul( T1_n1 );
        case 2                                      % Mode 2 (DIR tissue nulling)
            mainDIRNul( T1_n1, T1_n2 );
        case 3                                      % Mode 3 (T1W nulling)
            testT1_n2 = mainT1WNul( Tset(1,3), Tset(1,1), Tset(1,4) ); % ( T1_CSF, T1_WM, T1_MS )
    end
    createUI();
%     figure( gcf );                                  % Ensure focus on this figure? Doesn't work?
%     if saveImg, saveFigAsImg(),end                  % Automatically save the figure as an image (or use the fig. menu)?
end % main

% 1) Calculate TI in a 1IR/STIR/FLAIR sequence, nulling a specified T1 time
function mainSIRNul( T1z1 )
    T1z = T1z1;                                     % What T1 are we concerned with nulling?
    if ( sT2p == 0 )                                % Only use the T2prep formula if T2 is known
        E2 = 1;
    else
        E2 = exp( -sT2p/Tset(2,TisSel) );
    end % if sT2p
    TI1_n = T1z*log( (1+IE*E2)/(1 + IE*E2*exp(-sTR_ef/T1z))); % Calc. 1IR TI nulling T1
%     TI1_n = log(2)*T1z;                             % The simple case, TR>>T1, gives TIn = ln(2)*T1 ? 0.69*T1
    TI1_n = uint16( TI1_n );                        % Make value int for figure display
    
    TInts = [ 0, TI1_n+sT2p, sTR ];                 % SIR: Inversion is at t = 0, then readout at TI WIP: T2prep @ t=0?
    Mzt = PlotMagZ( TInts );                        % Plot Z magnetization over time, return end MagZ vector
    
    txt = ['Single IR tissue nulling at ' B0tag ', for TR_{eff} = ' num2str(sTR_ef) ' ms'];
    title(  txt,                                ...
            'FontSize', 16                      );
    
    fStr = ['\\itTissue T1 time:\n\\rm'         ...
            '  T1_{n}   = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling IR:\n\\rm\\bf'  ...
            '  TI_{n}    = %4.i ms\n\\rm'       ...
            '  TI_{n,c}  = %4.i ms'             ];
    fVar = [ uint16(T1z), uint16(TI1_n), uint16(TI1_n+sT2p) ];  % WIP: T2 prep
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    Ps.pStr = ['\n'                             ...
            '*******************************\n' ...
            '***    IR tissue nulling:   ***\n' ...
            '*******************************\n' ...
            '*  Tissue T1       = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  Nulling TI      = %4.i ms  *\n' ...
            '*******************************\n' ];
    Ps.pVar = [ T1z, TI1_n ];

end % mainSIRNul

% 2) Calculate TI1+TI2 in a DIR sequence, nulling two specified T1 times (T1_n1 = T1_CSF by default)
function mainDIRNul( T1z1, T1z2 )
    T1z = T1z2;                                      % What T1 are we concerned with nulling?
    [ TI1_n, TI2_n ] = Find2TIn2T1(T1z1,T1z2);    % Inversion times that null 2 T1
    
    TInts = [ 0, (TI1_n-TI2_n), TI1_n, sTR ];       % DIR: Vector of all time points (inversion, readout and end)
    Mzt = PlotMagZ( TInts );                        % Plot Z magnetization over time, return end MagZ vector

    txt = ['DIR dual tissue nulling at ' B0tag ', for TR_{eff} = ' num2str(sTR_ef) ' ms'];
    title(  txt,                                ...
            'FontSize', 16                      );
    
    fStr = ['\\itTissue T1 times:\n\\rm'        ...
            '  T1_{n1}  = %4.i ms\n'            ...
            '  T1_{n2}  = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling DIR:\n\\rm\\bf' ...
            '  TI_{1}    = %4.i ms\n'           ...
            '  TI_{2}    = %4.i ms'             ];
    fVar = [ T1z1, T1z2, TI1_n, TI2_n ];
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    Ps.pStr = ['\n'                             ...
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
    Ps.pVar = [ T1z1, T1z2, TI1_n, TI2_n, (TI1_n-TI2_n) ];

end % mainDIRNul

% 3) Calculate & plot MRI S0 of WM vs MS/WML in a DIR sequence, and TI1+TI2 nulling their T1 weighting
function T1z2 = mainT1WNul( T1z1, T1_T1, T1_T2 )    % T1_nulled = T1WNul( T1_CSF, T1_WM, T1_MS )
    T1z = Tset(1,TisSel);                     % What T1 are we concerned with nulling? TODO: Figure this out?! T1z1?
    TI2_max = log(IE2/(1 + IE*exp(-sTR_ef/T1z1)))*T1z1; % Max TI2 nulling T1z1, usually CSF (at max TI1 = TReff)
%     TI2_max = T1z1*log(2/(1 + 2*exp(-sTR/T1z1) - exp(-sTR_ef/T1z1))); % Old TI2_max calculation without InvEff
    TI2 = 1:1:TI2_max;                              % Vector of TI2 values 1 ms apart (could've used linspace() here)
    LEN = length(TI2);                              % Length = Array size in the biggest dimension
    TI1 = zeros(LEN,1);                             % Vector of TI1 values for all TI2s; replace by TI1s(TI2s)?
    S0_DIR = zeros(LEN,2);                          % Array of two DIR signals@TE=0 for TI2s
    
%     TI1s(t) = CalcTI1(t,T1z1);                      % TI1(TI2) calc.
    for i = 1:LEN
        TI1(i)      = CalcTI1D(TI2(i),T1z1);            % TI1(TI2) calc. for T1_CSF
        S0_DIR(i,1) = CalcS0([TI1(i),TI2(i)],T1_T1);    % S0(TI2) calc. (was CalcS0D(TI1(i),TI2(i),T1_T1) )
        S0_DIR(i,2) = CalcS0([TI1(i),TI2(i)],T1_T2);    % --"--
    end
    
%     S0_DIFF = abs(S0_DIR(:,1)-S0_DIR(:,2));         % Calculate signal difference between two T1s
%     [~,minind] = min(S0_DIFF);                      % Find the T1Wnull crossing graphically (as [minval,minind])
%     TI2_n = TI2(minind);
    TI2_n = uint16(FindT1WnTI2(T1_T1, T1_T2));      % Find [TI1,TI2] to null the DIR T1 weighting between two T1s
    TI1_n = uint16(CalcTI1D(TI2_n,T1z1));           % TI1(TI2) calc. by function
%     TI2_n = uint16(TI2_n);                        % Make figure text nicer w/o e+ notation
%     TI1_n = uint16(TI1s( TI2_n ));                % TI1(TI2) calc. by symbolic expr. (slower?!)
    
    T1_n2 = uint16(FindT1nulled(TI1_n,TI2_n));       % Second tissue nulled by DIR. UInt for display.
    T1z2 = T1_n2; % = T1z2
    newS = S0_DIR(TI2_n,1);                         % WM signal at readout; = WML signal here
    
    hold on
    plot( TI2, S0_DIR(:,1), 'r-',               ...
        'LineWidth', 1.0                        )   % Plot the signals calculated above
    plot( TI2, S0_DIR(:,2), 'b-' ,              ...
        'LineWidth', 1.0                        )
%     plot(TI2,TI1);                                  % Debug: Plot the TI1 calculated for each TI2
%     fplot( t, TI1s(t) );                            % Debug: Plot the symbolic TI1(TI2) - WIP
%     plot(TI2,S0_DIFF);                              % Debug: Plot the WM/MS signal difference
    hold off
    
    whitebg( 'white' );
%     whitebg( [ 0.70 1.00 0.70 ] );                  % Debug: See white plot/GUI elements
    txt =                                       ...
        ['T1-DIR (T1W nulling) at ' B0tag ', for TR_{eff} = ' num2str(sTR_ef) ' ms'];
    title(  txt,                                ...
            'FontSize', 16                      );
    ylabel('Rel. MR signal (TE = 0 ms)', 'FontSize', 12 )
    xlabel('TI_2 inversion time (ms)',   'FontSize', 12 )
    xlm = [ 0 200*ceil( TI2_max/200 ) ]; xlim(xlm); % Round up to 200 on the x axis
    ylm = [ -1 1 ];                      ylim(ylm); % Rel. S0 is in [-1,1]
    
    Tis = [ 1 4 ];                                  % Make a legend for tissue 1(WM) & 4(MS)...
    NT1s = length(Tis); LegS = strings(1,NT1s);
    for i = 1:NT1s, j = Tis(i);                     % ...showing the tissue short names and T1 times
        LegS(i) = sprintf('%s (T1 = %i ms)', Ttag(1,j), Tset(1,j) );
    end
    legend( LegS(1), LegS(2),                   ... % Was: legend( 'WM', 'MS', ...
            'Location', 'NorthWest'             );  % 'Best' may conflict with UI and TextBox
    
    lxlm = [ TI2_n TI2_n ];                         % S0 crossing x (=TI2) value
%     lylm = [ ylm(1) S0_DIR(lxlm(1),1) ];            % Lower/upper y of line (=S0 curve)
    lylm = [ ylm(1) 0 ];                            % Lower/upper y of line (= 0)
    line( lxlm, lylm,                           ... % Vertical line @ crossing
        'Color', 'black',                       ...
        'LineWidth', 1.25,                      ...
        'DisplayName', 'TI_{2} nullT1W',        ... % (Tip: If declared after legend, the line appears in it)
        'LineStyle', ':'                        );
    
%     fStr = ['\\itTissue T1 times:\n\\rm'        ...
%             '  T1_{WM} = %4.i ms\n'             ...
%             '  T1_{MS} = %4.i ms\n\n'           ...
    fStr = ['\\itT1W-0 DIR times:\n\\rm\\bf'    ... % Tex \bf\it\rm = bold/italic/normal; escape \ for sprintf!
            '  TI_{1}    = %4.i ms\n'           ...
            '  TI_{2}    = %4.i ms\n'           ...
            '  T1_{nt}  = %4.i ms'              ];
    fVar = [ TI1_n, TI2_n, T1z2 ];                  % Was: T1_WM, T1_MS, etc
    figTextBox( fStr, fVar );                       % Display an info text box on the figure

    Ps.pStr = ['\n'                             ...
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
    Ps.pVar = [ T1z1, T1z2, TI1_n, TI2_n, (TI1_n-TI2_n) ];
    
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
        S0_IR = S0_IR - (-IE)^i*IE2*exp(-Ti(NInv-i)./T1); % Note: TIv() are defined as time to readout!
    end
    S0_IR = S0_IR - (-IE)^NInv*exp(-sTR_ef./T1);    % One signal value per given T1
    % TODO: Include sT2p in the above! exp(-sT2p./T2) – need to pass T2 times too (or all times in one array)!
end % fn

% function S0_IR = CalcS0S( Ti, T1 )                  % Calculate single IR relative S0 (def: S0 = S/(E2*M0) at TE=0)
%     Ti = FoS(Ti); T1 = FoS(T1);                     % Cast input as float unless it's a sym
%     S0_IR = 1 - IE2*exp(-Ti/T1) + IE*exp(-sTR_ef/T1);
% end % fn

% function S0_DIR = CalcS0D( Ti1, Ti2, T1 )           % Calculate DIR relative S0 (def: S0 = S/(E2*M0) at TE=0)
%     Ti1 = FoS(Ti1); Ti2 = FoS(Ti2); T1 = FoS(T1);   % Cast input as float unless it's a sym
%     S0_DIR = 1 - IE2*exp(-Ti2/T1) + IE*IE2*exp(-Ti1/T1) - IE^2*exp(-sTR_ef/T1); % (New def.: TI1_new=TI1_old+TI2_old)
%     S0_DIR = CalcS0( [ Ti1, Ti2 ], T1 );            % New way: Generic IR signal function
% end % fn

function TI1_DIR = CalcTI1D( Ti2, T1 )              % Calc. TI1 (new TI1, the old one was - TI2) nulling T1 given TI2
    Ti2 = FoS(Ti2); T1 = FoS(T1);                   % Cast input as float unless it's a sym
    TI1_DIR = log(IE*IE2/(IE2*exp(-Ti2/T1) + IE^2*exp(-sTR_ef/T1) - 1))*T1; % This may be a sym, so don't cast to int
end % fn

function T1WnTI2 = FindT1WnTI2( T1w1, T1w2 )        % Find [TI1,TI2] to null the DIR T1 weighting between two T1s
    T1w1 = FoS(T1w1); T1w2 = FoS(T1w2);             % Cast input as float unless it's a sym
%     T1WnTI2 = vpasolve((CalcS0D(CalcTI1D(t,T1_n1),t,T1w1) ...   % Solve S0DIR == 0 for T1
%                       - CalcS0D(CalcTI1D(t,T1_n1),t,T1w2)), t); % Return TI2 (then calculate TI1(TI2) later)
    T1WnTI2 = vpasolve((CalcS0([CalcTI1D(t,T1_n1),t],T1w1) ...   % Solve S0DIR == 0 for T1
                      - CalcS0([CalcTI1D(t,T1_n1),t],T1w2)), t); % Return TI2 (then calculate TI1(TI2) later)
end % fn

function T1n_DIR = FindT1nulled( Ti1o, Ti2o )       % Find T1 time nulled by DIR
    Ti1o = FoS(Ti1o); Ti2o = FoS(Ti2o);             % Cast input as float unless it's a sym
%     T1n_DIR = vpasolve(CalcS0D(Ti1o,Ti2o,t), t)     % Solve S0DIR == 0 for T1
    T1n_DIR = vpasolve(CalcS0([Ti1o,Ti2o],t), t);   % Solve S0DIR == 0 for T1
end % fn

function [ Ti1n, Ti2n ] = Find2TIn2T1( T1n1, T1n2 ) % Find TI1/TI2 that null T1n1/2
    T1n1 = FoS(T1n1); T1n2 = FoS(T1n2);             % Cast input as float unless it's a sym
%     Ti2n = vpasolve(CalcS0D(CalcTI1D(t,T1n1),t,T1n2), t); % (Used VPAsolve for numeric output but still got syms?)
    Ti2n = vpasolve(CalcS0([CalcTI1D(t,T1n1),t],T1n2), t); % (Used VPAsolve for numeric output but still got syms?)
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

function SetMode( mode )                            % Sets the run mode and resets the Nulling selection (TisSel)
%     if ( mode == 0 ); mode = 1; end     % WIP: Make SetMode able to handle startup, to-from etc; merge with SetT1n()?!
    
    T1_n1 = Tset(1,3);                              % Default for all modes: Null out CSF as the longest T1
    switch mode
        case 1                                      % 1IR default: Null out CSF (=FLAIR)
            SetT1n( 1, 3 );
        case 2                                      % DIR default: Null out WM as tissue 2
            if ( Ps.Mode == 3 )
                SetT1n( 2, -1 );                    % A trick to notify SetT1n() that we've done T1 nulling
            else
                SetT1n( 2, 1 );
            end % if
        case 3                                      % T1W-DIR default: (Null out CSF as the longest T1)
    end % switch
    Ps.Mode = mode;                                 % Make the mode setting global
end % fn

function SetT1n( mode, T1s )                        % Sets the Nulling selection and T1_n#
    if ( T1s > 0 )
        switch mode
            case 1                                  % 1IR (only one tissue nulled)
                T1_n1 = Tset(1,T1s);
            otherwise                               % DIR (set second T1 to null; the first is CSF by default)
                if Tset(1,T1s) >= T1_n1, T1s=1; end % Don't null T1_n1/CSF twice; use WM as tissue 2 instead then.
                T1_n2 = Tset(1,T1s);
        end % switch
    else
        T1s = -T1s;
    end % if
    TisSel = T1s;
    
    switch mode                                     % Reset which T1 times to plot, based on mode. TODO: Change this?!?
        case 1                                      % 1IR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4;
        case 2                                      % DIR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4;
        case 3                                      % T1W-DIR: Plot the signal for WM and WML
            i = [ 1 4 ];                            % NOTE: This is not implemented yet; Mode 3 uses its own plot
    end % switch
%     T1Plt = Tset(1,i); TsLeg = Ttag(1,i);           % Update the T1 times to plot and their label texts
    TsPlt = Tset(:,i); TsLeg = Ttag(1,i);           % Update the T1 times to plot and their label texts
end % fn

function setRelaxTimes( B0opt )
    B0sel = B0opt;                                  % Ensures right selection in the UI
    B0str = [   "1.5 T (Alsop)", "3.0 T (Alsop)",                   ... % Array of B0/sources for relax. times
                "3.0 T (GE)", "1.5 T (MRIQue)", "3.0 T (Lalande)",  ...
                "1.5 T (Visser)", "7.0 T (Visser)"  ];  %   (used in the UI)
    B0tag = char(B0str(B0sel)); B0tag = B0tag(1:5);     % Char array of only the field strength, for figure titles
    TsStr = [   " WM     ", " GM     ",   ...           % Array of tissue names
                " CSF    ", " WML/MS ", " Fat    " ];   %   (for formatted output)
%     T1 = split(TsStr); T1 = char(T1(:,:,3));        % Array of only the tissue type (3rd word) for legends,...
%     TsTag = string(T1(:,1:3,:));                    % ...trimmed to length 3 (e.g., "CSF" or "WML")
    TsTag = [   "WM ",  "GM ",                  ... % Explicit declaration for sanity
                "CSF",  "WML",  "Fat"           ];
    
    switch B0opt                % Set T1 times based on B0 and sources
        case 1                  % 1.5 T (Alsop); Madhuranthakam et al, Mag Res Med 2012 67(1):81-88
            sB0    = 1.5;
            T1_WM  =  650.0;            % 1.5 T Alsop
            T1_GM  = 1300.0;            % 1.5 T --"--
            T1_CSF = 4200.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ? T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
            TsStr(5) = "(Fat)   ";      % (This value was cited from elsewhere)
        case 2                  % 3.0 T (Alsop)
            sB0    = 3.0;
            T1_WM  =  750.0;            % 3.0 T Alsop
            T1_GM  = 1400.0;            % 3.0 T litteratur
            T1_CSF = 4200.0;            % 3.0 T litteratur
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean)
            TsStr(2) = "(GM)    ";      % (This value was cited from elsewhere)
            TsStr(5) = "(Fat)   ";      % (This value was cited from elsewhere)
        case 3                  % 3.0 T (GE)
            sB0    = 3.0;
            T1_WM  =  825.0;            % 3.0 T litteratur (GE har 825 ms, andre 850 ms?)
            T1_GM  = 1400.0;            % 3.0 T --"--
            T1_CSF = 4200.0;            % 3.0 T --"--
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean)
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
            TsStr(5) = "(Fat)   ";      % (This value was cited from elsewhere)
        case 4                  % 1.5 T (MRIQ); MRIquestions.com and other online resources
            sB0    = 1.5;
            T1_WM  =  580.0;            % 1.5 T MRIQ
            T1_GM  =  940.0;            % 1.5 T --"--
            T1_CSF = 3600.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ? T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQ
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
        case 5                  % 3.0 T (Lalande); Bojorquez et al, MRI 35(2017) 69–80
            sB0    = 3.0;
            T1_WM  =  940.0; %842.0;    % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR); Stikov (2015) has 940(860-992) ms.
            T1_GM  = 1425.0;            % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR)
            T1_CSF = 4300.0;            % 3.0 T Lalande (trimmed mean of 4(6) reports)
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean; 405.0 in Barral 2010 w/ SE-IR)
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
        case 6                  % 1.5 T (Visser); Visser 2010 MagResMed
            sB0    = 1.5;
            T1_WM  =  656.0;            % 1.5 T Rooney
            T1_GM  = 1188.0;            % 1.5 T --"--
            T1_CSF = 4329.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ? T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
            TsStr(5) = "(Fat)   ";      % (This value was cited from elsewhere)
        case 7                  % 7.0 T (Visser); Visser 2010 MagResMed
            sB0    = 7.0;
            T1_WM  = 1220.0;            % 7.0 T Rooney
            T1_GM  = 2132.0;            % 7.0 T --"--
            T1_CSF = 4400.0;            % Estimated from the 1.5 T value...
            T1_MS  =  T1_GM;            % Estimate T1_WML ? T1_GM
            T1_Fat =  600.0;            % 7.0 T estimated(?!?)
            TsStr(3) = "(CSF)   ";      % (This value was cited from elsewhere)
            TsStr(4) = "(MS/WML)";      % (This value was cited from elsewhere)
            TsStr(5) = "(Fat)   ";      % (This value was cited from elsewhere)
        otherwise
            error('ERROR: Undefined B0/source!');
    end % switch B0opt

    %   case 6          % 1.5 T (Visser); Visser 2010 MagResMed – used here as fallback values (1.5 & 3.0 T)
            T2_WM  =   74.0;            % 1.5 T Wehrli, Yacoub
            T2_GM  =   87.0;            % 1.5 T Wehrli, Yacoub
            T2_CSF = 2280.0;            % 1.5 T Helms
            T2_MS  =  100.0;            % 1.5 T Alsop
            T2_Fat =   40.0;            % NB: PURE GUESSWORK! BEWARE!
    switch B0opt                % Set T2 times based on B0 and sources, if we found any
        case 1                  % 1.5 T (Alsop); Madhuranthakam et al, Mag Res Med 2012 67(1):81-88
            T2_WM  =   70.0;            % 1.5 T Alsop
            T2_CSF = 2000.0;            % 1.5 T Alsop
            T2_MS  =  100.0;            % 1.5 T --"--
        case 5                  % 3.0 T (Lalande); Bojorquez et al, MRI 35(2017) 69–80
            T2_WM  =   75.0;            % 3.0 T (Lu 2005; jMRI 22(1):13–22)
            T2_GM  =   83.0;            % 3.0 T Lalande
        case 7                  % 7.0 T (Visser); Visser 2010 MagResMed
            T2_WM  =   46.0;            % 7.0 T Wehrli, Yacoub (= a factor 0.63 from 1.5 T)
            T2_GM  =   55.0;            % 7.0 T Wehrli, Yacoub (= a factor 0.63 from 1.5 T)
            T2_CSF = 1800.0;            % Estimated from the 1.5 T value (> factor 0.63...)
            T2_MS  =   70.0;            % Estimated from the 1.5 T value (> factor 0.63...)
    end % switch B0opt
    
    Tset  = [ [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ];  ...
              [ T2_WM, T2_GM, T2_CSF, T2_MS, T2_Fat ] ];
    Ttag  = [ TsTag; TsStr ];
    T1_n1 = T1_CSF;                     % Default for all modes: Null out CSF as the longest T1
    switch Ps.Mode                      % When a scheme is selected, adjust the nulling selection accordingly.
        % TODO: Remove this switch; just call SetT1n( TisSel ) and let the function handle Ps.Mode!
        case 1                          % 1IR T1 nulling
            SetT1n( 1, TisSel );
        case 2                          % DIR T1 nulling
            SetT1n( 2, TisSel );
        case 3                          % T1W nulling
    end % switch
end % setRelaxTimes

%% PLOT

% Plot the magnetization Mz of several T1s under multiple inversion time points (NB: Absolute TI times are used here!)
function MagZ = PlotMagZ( TIabs )                   % ( TIabs, T1ts, T1cs ); TIabs = [ dTI1 dTI1+dTI2 Tacq DUR ]
    NIT = length(TIabs)-2;                          % # of inversion time points, without readout time and duration
    TIn = TIabs(1:NIT);                             % The actual inversion time points
    DUR = TIabs(NIT+2);                             % Total duration for the plot
    TIR = TIabs(NIT+1); xRo = double([ TIR TIR ]);  % Readout time (=TI)
%     sZRo = 0.02*sETD;                               % NOTE: sZRo is irrelevant as spin lock dictates readout from TIR?
%     TRo = xRo + [ -sZRo, (sETD-sZRo) ];             % Readout time range. Linear(sZRo=0.5) or center-out (sZRo small)?
    TRo = xRo + [ 0, sETD ];                        % Readout time range. Due to spin lock, this should start at TIR.
    dT = 1.0; Tp = (0:dT:DUR-dT);                   % Sampling interval dT (ms, Int), Tp is the time point vector
                                                    % TODO: Allow non-integer sampling intervals?
    
    T1Plt = TsPlt(1,:); T2Plt = TsPlt(2,:);
    Ntis = length( TsPlt(1,:) );                    % # of relaxation times for which to calculate/plot
%     if ( sT2p == 0 )                                % WIP: Plotting extra tissues is hard if T2 in T2p is unknown!
%         if ~ismember( T1z, T1Plt )                     % Add T1z if not in the preselected T1 ensemble
%             Ntis  = Ntis + 1;
%             T1Plt = double( [ T1Plt, T1z ] );
%             T2Plt = double( [ T2Plt, T1z ] );       % WIP: Estimate T2 from known T1-T2 relationships (or look it up)?!?
%             TsLeg = [ TsLeg, sprintf('T1#%s', num2str(Ntis)) ];
%         end % if ~ismember T1z
%     end % if sT2p
    Mz = ones( Ntis, DUR );                         % MagZ for all tissues and time points (was zeros(NT1s+1...)
    
    % WIP: Starting Z magnetization calculation.
%     Mz0 = ones(1,length(T1Plt));                   % Test/debug: Start with relaxed signal (TR >> T1)
%     Mz0 = -IE*Mz0;                                  % --"--: Start with inversion
%     Mz0 = 0*Mz0;                                    % --"--: Start with saturation
%     Mz0 = -CalcS0( [ sTR_ef, 0 ], [ T1Plt ] )      % --"--: CalcS0 uses TI differences.
%     Mz0 = CalcMagZ( [ 0 ], sTR_ef, T1Plt )         % --"--: [1 - exp(-sTR_ef./T1Plt)]
    % TODO: Use derivation instead of iteration to arrive at SS? Maybe not so critical. Also, nice to simulate it.
    RTI = zeros( 1, sETL);
    for i = 1:sETL
        RTI(i) = TRo(1) + round(i*sESP);            % Array of readout refocus times. TODO: Avoid rounding?
    end % for
    
    its=3; MzIt = ones(Ntis,its);                   % Iterate to get to a steady state
    for i = 1:its                                   % Iteration loop (Note: Needs only two TR if Mz(TRo)=0)
        Mz(:,1) = Mz(:,DUR);                            % The starting magnetization (before first inv.) is MagZ(TR)
        for j = 2:DUR                                   % Step through all time points
            Ts = (j-2)*dT;                              % Was (j-1)*dT, but I want the time point 0 to be defined
            Mz1 = Mz(:,j-1);                            % MagZ at the previous time point
            if ismember( Ts, TIn )                      % Inversion time, so invert the Mz (accounting for Inv.Eff.)
                Mz(:,j) = -IE*Mz1;                      % (was any(Ts == TIn(:)) before)
            elseif ( Ts == TRo(1) )                     % Readout time
                MzTI    = Mz1;                          % MagZ at readout
                SRo     = sin(sFA*pi/180)*Mz1;          % Signal strength at readout
                Mz(:,j) = cos(sFA*pi/180)*Mz1;
                MzRo    = MzTI - Mz(:,j);
%             elseif ( Ts > TRo(1) ) && ( Ts < TRo(2) )   % Readout with CPMG train, so signal doesn't evolve (so much)
            elseif ismember( Ts, RTI )                  % CPMG refocus time (rounded to the nearest ms...)
                Mz(:,j) = -sIET*Mz1;
%                 for k = 1:NT1s                        % TODO: Evolve slowly during readout with adjusted T1s? Or?!?
%                     Mz(k,j) = 1 - (1 - Mz(k,j-1))*exp(-dT*T1rhF/T1Plt(k));
%                 end % for k
            elseif ( Ts == TRo(2) )                     % Readout end
                Mz(:,j) = Mz1 + sRew*MzRo;              % A rewinder reclaims some excitation signal as Mz? Or can it?!?
%             elseif ( Ts > sTR - sT2p )                  % T2 preparation (effectively decays Mz with T2)
            elseif ( Ts < sT2p )                        % T2 preparation (effectively decays Mz with T2)
                for k = 1:Ntis                          % Can't assign all tissues at once with exp, so use a loop
                    Mz(k,j) = Mz(k,j-1)*exp(-dT/T2Plt(k)); % 
                end % for k
            else                                        % Signal evolvement ( if(t~=TI(:)) )
                for k = 1:Ntis                          % Can't assign all tissues at once with exp, so use a loop
    %                 Mz(k,j) = Mz(k,j-1)*exp(-dT/T1v(k)) + (1-exp(-dT/T1v(k)));
                    Mz(k,j) = 1 - (1 - Mz(k,j-1))*exp(-dT/T1Plt(k));
                end % for k
            end % if Ts
        end % for j
        MzIt(:,i) = Mz(:,DUR);                      % The end Z magnetization for all T1s, for this iteration
    end % for i
    MagZ = Mz(:,DUR)';                              % Return the end Z magnetization for all T1s as a vector
    
%     Mz = [ Mz; zeros(1,DUR) ];                    % Add a zero line to the plot
    plot( Tp, Mz, '-',                          ... % Plot the (steady state) Z magnetization
        'LineWidth', 1.5                        )
    
    xlabel( 't (ms)'      , 'FontSize', 12 )
    ylabel( 'Rel. Z magn.', 'FontSize', 12 )
    xlm = [  0  DUR ]; xlim(xlm);
    ylm = [ -1   1  ]; ylim(ylm);
%     set( gca, 'XTick', [] );                        % Remove x axis ticks
%     set( gca, 'YTick', [] );                        % Remove y axis ticks

    line( xlm, [ 0 0 ],                         ... % Zero line
        'Color', 'black',                       ...
        'LineWidth', 0.7,                       ...
        'LineStyle', ':'                        );
    line( TRo, [ -0.975 -0.975 ],               ... % Horizontal line representing readout
        'DisplayName', 'Echo train',            ...
        'Color', [ 0.3 0.6 0.3 ],               ... % 'black'
        'LineWidth', 8,                         ...
        'LineStyle', '-'                        );
    
    line( [ 0 sT2p ], [ -0.975 -0.975 ],        ... % Horizontal line representing T2prep
        'DisplayName', 'T2 prep.',              ...
        'Color', [ 0.9 0.6 0.3 ],               ... %
        'LineWidth', 8,                         ...
        'LineStyle', '-'                        );
    
    for i = 1:NIT
    line( [ TIn(i) TIn(i) ], [ -1.00 -0.92 ],   ... % Vertical line representing inversion
        'DisplayName', 'Inversion',             ...
        'Color', [ 0.5 0.3 0.7 ],               ...
        'LineWidth', 5,                         ...
        'LineJoin', 'chamfer',                  ... % Note: This doesn't round corners unless they're joined.
        'LineStyle', '-'                        );
    end % for i
    
%     line( xRo, ylm,                             ... % Vertical line at TEeff (=TI)
%         'Color', 'black',                       ... % (Tip: If declared after legend, the line appears in it)
%         'DisplayName', 'TI_{1}',                ... % 'TE_{eff}'
%         'LineStyle', ':'                        );
    
    LegS = strings(1,Ntis);                             % Figure legends
    for i = 1:Ntis
        LegS(i) = sprintf('%4s (T1: %4i ms)', TsLeg(i), TsPlt(1,i) );   %legend('WM', 'GM', 'CSF'...
    end % for i
    legend( LegS, 'Location', 'NorthWest');
end % PlotMagZ

function figTextBox( fStr, fVar )
    txt = sprintf( fStr, fVar );                    % Formatted string "print"
    lin = count( fStr, '\n' ) + 1;                  % Count lines in the string
    xlm = xlim; xt =  0.800*double(xlm(2));         % Max value of axis; (x,y) is a plot point!
    ylm = ylim; yt = -0.900*double(ylm(2));         % --"--
    yt = yt + 0.035*lin;                            % Compensate for number of lines
    text( xt, yt, txt,                          ... % Note: annotation('textbox' etc doesn't scale with text
        'EdgeColor', 'black',                   ...
        'Margin', 4,                            ...
        'LineWidth', 1,                         ...
        'LineStyle', '-'                        );
end % fn

% function saveFigAsImg()
%     imgName = ['TReff_' int2str(sTR_ef) '_TI1_' int2str(TI1_n) '_TI2_'  int2str(TI2_n) ];
%     if exist('Figures','dir')~=7, mkdir('Figures'), end
%     cd Figures
%     % TODO: What's the difference between print to image and the imwrite command?
%     print('-djpeg', '-r300', imgName )
%     cd '..'
% end % fn

%% USER INTERFACE

function createUI()                                 % Display a button to show/hide the UI controls
%     uiX = 74.0*fSz; uiY = 50.0*fSz;                 % Default start position (in px) for masterUI controls
    eW = 0.065 ; eH = 0.065; eS = eH + 0.025;       % Button width/height/spc for masterUI controls (was [50 40 60] px)
    myAx = gca;                                     % A handle to the current graphics axes (plot)
    uiX = myAx.OuterPosition(3) - 0.015 - eW;       % Use the axis' outer x position (def. [ 0 0 1 1 ]) for UI pos.
    uiY = myAx.Position(4)      + 0.010     ;       % Use the axis' inner y position to determine UI pos.
    
    uicontrol( 'Style', 'pushbutton',           ... % For all modes:
        'String', 'Panel',                      ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[s] Show/hide UI',    ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-0*eS eW eH ],     ...
        'Callback', @masterUI_callback          );  % UI to show/hide the UI control panel
    
    switch Ps.Mode
        case 1                                      % 1IR T1 nulling
            b2txt = 'DIR';
            b4txt = '';
        case 2                                      % DIR T1 nulling
            b2txt = '(FLA)IR';
            b4txt = 'T1-null';
        case 3                                      % T1W nulling
            b2txt = '';
            b4txt = 'Back';
    end % switch Ps.Mode
    if ismember( Ps.Mode, [ 1 2 ] )                   % If 1IR/DIR:
    uicontrol( 'Style', 'pushbutton',           ... % uiMode
        'String', b2txt,                        ... % 'togglebutton'
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[m] Switch mode',     ...
        'Units', 'normalized',                  ... % Make the UI resize automatically with the figure
        'Position', [ uiX uiY-1*eS eW eH ],     ...
        'Callback', {@cMode_callback,[ 2 1 0 ]} );  % UI to set run mode (for mode 1–3, switch to...)
    
    uicontrol( 'Style', 'pushbutton',           ... % If 1IR/DIR:
        'String', 'Plot T2',                    ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[p] T2 decay plot',   ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-2*eS eW eH ],     ...
        'Callback', @plotT2_callback            );  % UI to create T2 vs TE plot
    end % if
    
    if ismember( Ps.Mode, [ 2 3 ] )                   % If DIR/T1W_n:
    uicontrol( 'Style', 'pushbutton',           ...
        'String', b4txt,                        ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[t] Switch mode',     ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-3*eS eW eH ],     ...
        'Callback', {@cMode_callback,[ 0 3 2 ]} );  % UI to change run mode to T1W nulling plot
    end % if
    
    createUIPanel();
    switch UIs.Show
        case true ; UIsBox.Visible = 'on';
        case false; UIsBox.Visible = 'off';
    end % switch UIs.Show
end % createUI

function createUIPanel()
    % Make UI control buttons for saving images instead of the variable, like with printVars?
    %      Probably not necessary, since you can save the image from the File menu.
    
    boxClr = figClr;                                % UI panel background color (same as fig background)
%     uiX = 58.3*fSz; uiY = 54.5*fSz;                 % Default start position (in px) for UI controls
    mX = 10; mY =  6;                               % x/y margins in UI panel
    eN = 11; eW = 120; eH = 30;                     % Number of UI elements; default element width/height
    UIs.Hpx = eN*eH + 2*mY; uiT = UIs.Hpx - mY +  3;    % UI panel height in pixels; top position for elements
    UIs.Wpx = eW + 2*mX;                              % UI panel width in pixels
    
    UIsBox = uipanel( Ps.Fig1,                  ... % UIsBox
        'BackgroundColor', boxClr,              ... %  [ 0.8 0.8 0.4 ]
        'BorderType', 'etchedin',               ... % default 'etchedin'
        'BorderWidth', 1.0,                     ... % default 1
        'Clipping', 'off',                      ... % default 'on'; 'off' is good for debugging
        'Units', 'normalized',                  ... % Normalized units facilitate resizing to axes
        'ResizeFcn', @placeUiBox                );  % Resize function for the UI panel
%         'Units', 'pixels',                      ... % With normalized units (def.), got Pos. [ .726 .420 .172 .500 ]
%         'Position', [ uiX uiY-uiH uiW uiH ],    ... % UI panel container for the controls (pixel units)
%         'Title', 'Controls',                    ... % If specifying a title, it'll show up on the panel
%         'FontSize', 10,                         ...
%     uistack( UIsBox, 'top' );                       % Bring the panel to the top of the graphic stack (not necessary)
    placeUiBox();                                   % Adjust the UI panel position to the figure window
    
    for i = 1:eN; eY = uiT - i*eH;                    % Generate each UI element on a new row in the panel
        switch i
        case  1
        uicontrol( 'Parent',    UIsBox,             ... % ui pT1s
            'Position', [ mX eY+0 eW 24 ],          ...
            'Style',            'pushbutton',       ...
            'String',           'Info',             ...
            'ToolTipString',    '[f] Print info',   ... % 'Print T1 times'
            'Callback', @printInfo_callback         );  % UI to print settings to the command window

        case 2
        uicontrol( 'Parent',    UIsBox,             ... % ui Ps.pVar
            'Position', [ mX eY+0 eW 24 ],          ...
            'Style',            'pushbutton',       ...
            'String',           'Results',          ...
            'ToolTipString',    '[w] Print results',...
            'Callback', @printVars_callback         );  % UI to print results to the command window
        
        case  3
        uicontrol( 'Parent', UIsBox,                ... % ui TR[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Repetition time',     ...
            'String', 'TR:                             ms');
        UIs.TRS = uicontrol( 'Parent', UIsBox,      ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'String',   sTR,                        ...
            'Value',    sTR,                        ...
            'Callback', @tr_set_callback            );  % UI to set TR

        case  4
        uicontrol( 'Parent', UIsBox,                ... % ui TRo[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'TSE Readout time',    ...
            'String', 'T_Ro:                         ms');
        uicontrol( 'Parent', UIsBox,                ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'ToolTipString', '= ETL * ES',          ...
            'String',   sETD,                       ...
            'Value',    sETD,                       ...
            'Callback', @tro_set_callback           );  % UI to set T_Read

        case  5
        uicontrol( 'Parent', UIsBox,                ... % ui IE[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Inversion efficiency',...
            'String', 'Inv.eff.:                      %');
        uicontrol( 'Parent', UIsBox,                ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'String',   sIEf,                       ...
            'Value',    sIEf,                       ...
            'ToolTipString', '-cos(FA_inv)',        ...
            'Callback', @ie_set_callback            );  % UI to set inversion efficiency

        case  6
        uicontrol( 'Parent', UIsBox,                ... % ui FAex[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Excitation pulse FA', ...
            'String', 'FAex:                         °'); % Note: UIControl text can't support Tex/LaTeX
        UIs.FAexS = uicontrol( 'Parent', UIsBox,    ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 10000, 'Max', 10000,             ...
            'String',   sFA,                        ...
            'Value',    sFA,                        ...
            'Callback', @fa_set_callback            );  % UI to set FA (in degrees)

        case  7
        uicontrol( 'Parent', UIsBox,                ... % ui T2p[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'T2 preparation dur.', ...
            'String', 'T2prep:                      ms');
        UIs.T2pS = uicontrol( 'Parent', UIsBox,     ...
            'Position', [ mX+40 eY+0 60 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 10000, 'Max', 10000,             ...
            'ToolTipString', 'WIP: Unfinished!',    ...
            'String',   sT2p,                       ...
            'Value',    sT2p,                       ...
            'Callback', @t2p_set_callback           );  % UI to set T2prep duration (in ms)
        UIs.T2pS.ForegroundColor  = myGray;             % Gray out the T2prep text for now since it's WIP

        case  8
        UIs.B0T = uicontrol( 'Parent', UIsBox,      ... % ui B0[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Ref. values for T1s', ...
            'String', 'B0:'                         );
        UIs.B0S = uicontrol( 'Parent', UIsBox,      ...
            'Position', [ mX+22 eY+0 eW-22 23 ],    ... % Note: Height setting won't change the actual height here
            'Style', 'popup',                       ...
            'String',   B0str,                      ...
            'Value',    B0sel,                      ...
            'Callback', @b0_sel_callback            );  % UI to set B0 (relaxation times)

        case  9
        UIs.Tn1T = uicontrol( 'Parent', UIsBox,     ... % ui Tn1[T|S]
            'Position',     [ mX eY+1 eW 16 ],      ...
            'Style',        'text',                 ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Longer T1 to null',   ...
            'String', 'T1_n1:                       ms'); % Note: UIControl text can't support Tex/LaTeX
        uicontrol(              'Parent', UIsBox,       ...
            'Position',         [ mX+40 eY+0 60 20 ],   ...
            'Style',            'edit',                 ...
            'Min', 10000,       'Max', 10000,           ...
            'String',           T1_n1,                  ...
            'Value',            T1_n1,                  ...
            'Callback',         @t1n1_set_callback      );  % UI to set T1_n1 (longer T1, typically CSF)

        case  10
            if ismember( Ps.Mode, [ 2 3 ] )                   % DIR T1 nulling and T1-DIR null two T1
            UIs.Tn2T = uicontrol( 'Parent', UIsBox,     ... % ui Tn2[T|S]
                'Position',     [ mX eY+1 eW 16 ],      ...
                'Style',        'text',                 ...
                'BackgroundColor', boxClr,              ...
                'HorizontalAlignment', 'Left',          ...
                'ToolTipString', 'Shorter T1 to null',  ...
                'String', 'T1_n2:                       ms');
            UIs.Tn2S = uicontrol( 'Parent', UIsBox,     ...
                'Position',     [ mX+40 eY+0 60 20 ],   ...
                'Style',        'edit',                 ...
                'Min', 10000,   'Max', 10000,           ...
                'String',       T1_n2,                  ...
                'Value',        T1_n2,                  ...
                'Callback',     @t1n2_set_callback      );  % UI to set T1_n2 (shorter T1; GM/WM/etc)
            end % if

        case 11
            UIs.STnT = uicontrol( 'Parent', UIsBox,     ... % ui STn[T|S]
                'Position',     [ mX eY+1 eW 16 ],      ...
                'Style',        'text',                 ...
                'BackgroundColor',      boxClr,         ...
                'HorizontalAlignment',  'Left',         ...
                'ToolTipString', 'Ref. tissue to null', ...
                'String',       'Null T1:'              );
            if ismember( Ps.Mode, [ 1 2 ] )                   % 1IR and DIR T1 nulling use tissue selector
            UIs.STnS = uicontrol( 'Parent', UIsBox,     ...
                'Position',     [ mX+40 eY+0 eW-40 22 ],...
                'Style',        'popup',                ...
                'String',       Ttag(2,:),              ...
                'Value',        TisSel,                 ...
                'Callback',     @t1_sel_callback        );  % UI to set T1_n(2) by tissue
            if T1z ~= Tset(1,TisSel)
                UIs.STnS.ForegroundColor  = myGray;         % Gray out the tissue selector if entering values manually
            end % if T1z
            else
                UIs.STnT.String = 'Rel. S0:                     %';
                if isempty(oldS); oldS = newS; end % if         % For use in the relSNR uicontrol box
                if isempty(oldB); oldB = sB0 ; end % if         % 
                if ( debugInfo )
                    fprintf("\t"); fprintf("%-10s", ["new S0","old S0","(rel. units)"] ); fprintf("\n");     % Legends
                    disp( [ newS, oldS ] )
                end % if debug
                relS0 = uint16(100*(newS/oldS)*(sB0/oldB));     % Relative S0 is proportional to CNR, and rel. SNR here
                oldS = newS;
                uicontrol(      'Parent', UIsBox,       ... % Text box to show rel. S0 (instead of tissue selector)
                'Position',     [ mX+40 eY+0 60 20 ],   ...
                'Style',        'edit',                 ...
                'Min', 10000,   'Max', 10000,           ...
                'String',       relS0,                  ...
                'Enable',       'Inactive'              );
            end % if

% TODO: Make a button to show a list box for selecting what to plot, in a separate GUI panel?
%         case 11
%             UIs.ShTT = uicontrol( 'Parent', UIsBox,     ... % ui ShT[T|S]
%                 'Position', [ mX eY+1 eW 16 ],          ...
%                 'Style', 'text',                        ...
%                 'BackgroundColor', boxClr,              ...
%                 'HorizontalAlignment', 'Left',          ...
%                 'ToolTipString', 'Tissues to show', ...
%                 'String', 'Showing:'                    );
%             UIs.ShTS = uicontrol( 'Parent', UIsBox,     ...
%                 'Position', [ mX+40 eY-60 eW-40 82 ],   ...
%                 'Style', 'listbox',                     ...
%                 'Min',0,'Max',1,                        ... % trick to select multiple listbox items
%                 'String',   Ttag(1,:),                  ...
%                 'Value',    TsLeg,                      ...
%                 'Callback', @t_disp_callback            );  % UI to set which tissues to display/calculate

         end % switch i
    end % for i
    
    switch Ps.Mode
        case 1                                      % 1IR T1 nulling
        UIs.Tn1T.String  = 'T1_n:                         ms';  % There is only one T1 to null in 1IR
        UIs.STnT.Position = UIs.STnT.Position + [ 0 eH 0 0 ];   % Move the tissue selector one line up
        UIs.STnS.Position = UIs.STnS.Position + [ 0 eH 0 0 ];   % --"--
%         UIs.Tn2S.Visible = 'off';                   % Hide the tissue 2 edit text/box (ui###T/S, if created)
        case 2                                      % DIR T1 nulling
%         UIs.B0S.ForegroundColor = myGray;           % Gray out the B0 select text/box (ui###T/S, if created)
        case 3                                      % T1W nulling
        UIs.Tn2S.Enable           = 'inactive';
        UIs.Tn2S.ForegroundColor  = myGray;         % Gray out the tissue 2 edit text (it's set automatically)
        UIs.FAexS.Enable          = 'off';          % Gray out the FA edit box (FA is irrelevant in this mode)
    end % switch Ps.Mode

end % createUIPanel

function val = str2val( str )                       % Function to process numeric input for UI controls
    val = str2double( str );
    if ( val ~= val )                               % test for NaN from malformed input
        val = 0;
    end % if
end % fn

function masterUI_callback(~,~)                     % UI that shows/hides the UI control panel
    UIs.Show = ~UIs.Show;
    switch UIs.Show
        case true ; UIsBox.Visible = 'on';
        case false; UIsBox.Visible = 'off';
    end
end % fn

function cMode_callback(~,~,toMode)                 % Switch between calculation modes
    switch Ps.Mode
        case 1                                      % 1IR T1 nulling    ->|
            SetMode( toMode(1) );
        case 2                                      % DIR T1 nulling    |<-
            SetMode( toMode(2) );
        case 3                                      % T1W nulling       <<-
            SetMode( toMode(3) );
    end                                             % Note: The figure title shows the active mode
    main();
end % fn

function plotT2_callback(~,~)                       % Run T2 plot (in separate script)
    IR_T2_Subplot( SRo(1:4), Tset(2,:), Ttag(1,1:4), sETD, [ 4 1 2 ] )   % S0/T2/tags for a set of tissues to plot
    % TODO: Make this more flexible and sensible!
end % fn

function b0_sel_callback(src,~)                     % UI to set B0 (rel. times)
    setRelaxTimes( src.Value );                     % Was get(src,'Value') pre-R2014b
    main();
end % fn

function tr_set_callback(src,~)                     % UI to set TR
    sTR = str2val( src.String );                    % str2double( src.String );
    Lim = sETD + 5; if sTR < Lim, sTR = Lim; end    % TR > T_Read
    main();
end % fn

function tro_set_callback(src,~)                    % UI to set T_Read
    sETD = str2val( src.String );
    Lim = sTR - 5; if sETD > Lim, sETD = Lim; end % TR > T_Read
    main();
end % fn

function t1n1_set_callback(src,~)                   % UI to set T1_n1 (longer T1)
    T1_n1 = str2val( src.String );
    if Ps.Mode > 1                                  % Using multiple inversions...
        Lim = T1_n2 + 1;
        if T1_n1 < Lim, T1_n1 = Lim; end            % T1_n1 is the longest T1 time to null
    end
    main();
end % fn

function t1n2_set_callback(src,~)                   % UI to set T1_n2 (shorter T1)
    T1_n2 = str2val( src.String );
    Lim = T1_n1 - 1;
    if T1_n2 > Lim, T1_n2 = Lim; end                % T1_n1 is the longest T1 time to null
    main();
end % fn

function t1_sel_callback(src,~)                     % UI to select a tissue T1 to null
    SetT1n( Ps.Mode, src.Value );
    main();
end % fn

function fa_set_callback(src,~)                     % UI to select the flip angle (in degrees)
    sFA = str2val( src.String );
    main();
end % fn

function t2p_set_callback(src,~)                    % UI to select the T2 preparation time (in ms)
    sT2p = str2val( src.String );
    main();
end % fn

function ie_set_callback(src,~)                     % UI to select the inversion efficiency (in percent)
    sIEf = str2val( src.String );
    main();
end % fn

function printVars_callback(~,~)                    % UI that prints results to the command window
    fprintf( Ps.pStr, Ps.pVar );
    if debugInfo && ismember( Ps.Mode, [ 1 2 ] )    % Debug info for plot
%         len = length(MzIt);
%         fprintf('\nMagZ(TR) over the first %i repetitions (starting at 1):\n', len);
%         disp( MzIt );                               % Show the iterative Z magnetization over the first TRs
%         disp( MzIt(:,2:its+1) );
        fStr = [ '\nMagZ and Signal(0) by T1 at the'...
            ' last readout with flip angle %3.1f°'  ...
            ' [S0=%1.3f*MagZ(TI)]:\n'               ];
        fprintf(fStr, round(sFA,1,'significant'), sin(sFA*pi/180)); % Heading
        fprintf("\t "); fprintf("%-10s", TsLeg); fprintf("\n");     % Tissue legends
        if ( sFA == 90 )
            fStr =          SRo'  ;
        else
            fStr = [ MzTI'; SRo' ];
        end % if
        disp ( fStr );                              % Mz/S (after FA° pulse) at readout
    end % if
end % fn

function printInfo_callback(~,~)                    % UI that prints settings to the command window
    fStr = ['\n'                                ...
            '*******************************\n' ... % fprintf( fStr, B0str(B0sel), T1set );
            '***    Relaxation times     ***\n' ...
            '***    for %-15s',       '  ***\n' ...
            '*******************************\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ... % [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ]
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*******************************\n' ...
            '*         Settings:           *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s %%   *\n'... % Used %4.i when printing T1s directly
            '*******************************\n' ]; 
    T1info = [ Ttag(2,1:5); num2cell(Tset(1,1:5)) ]; % Use this trick to weave strings and numbers into a string array
    T1pars = [ ["TR    ","T_Readout","Inv.Eff."]; num2cell([ sTR,sETD,sIEf]) ];
    fprintf( fStr, B0str(B0sel), T1info{:}, T1pars{:} );
end % fn

function keyPress_callback(~,evt)
    switch evt.Key
        case 's'                                    % Show/hide UI panel
            masterUI_callback([],[])                % Calling with dummy values
        case 'p'                                    % Plot T2 decay
            plotT2_callback([],[])
        case 't'                                    % T1-nulling, plot T2 decay, focus back on this fig.
            if ismember( Ps.Mode, [ 2 3 ] )
                cMode_callback([],[],[ 0 3 2 ])
                cMode_callback([],[],[ 0 3 2 ])
                plotT2_callback([],[])
                figure( Ps.Fig1 );
            end % if
        case 'f'                                    % Print info
            printInfo_callback([],[])
        case 'w'                                    % Print results
            printVars_callback([],[])
        case 'm'                                    % Change mode
            if ismember( Ps.Mode, [ 1 2 ] )
                cMode_callback([],[],[ 2 1 0 ])
            end % if
        case 'a'                                    % Change active figure
            if ishandle( Ps.Fig2 )
                figure( Ps.Fig2 );
            end % if
        case 'q'                                    % Focus on TR setting uicontrol
            uicontrol( UIs.TRS );
%         case 'c'                                    % Change focus to this figure (WIP)
%             figure( Ps.Fig1 );
    end % switch
end % fn

function placeUiBox(~,~)                            % Function to resize the UI panel (automatically called)
    if ishandle(UIsBox)                             % (The first time this fn is called, the panel doesn't exist yet)
        pxFig = getpixelposition(Ps.Fig1);          % Figure size in px, for manual normalizing
        myAx = gca;                                 % A handle to the current graphics axes/plot
        uiW = UIs.Wpx / pxFig(3);                   % Normalized UI panel width
        uiH = UIs.Hpx / pxFig(4);                   % Normalized UI panel height
        uiX = myAx.Position(3) - uiW + 0.115;       % The axis' inner right x pos. (but why the "fudge factor"?)
        uiY = myAx.Position(4) - uiH + 0.090;       % The axis' inner upper y pos. (--"--)
        UIsBox.Position = [ uiX uiY uiW uiH ];
    end % if
end % fn

end % script fn