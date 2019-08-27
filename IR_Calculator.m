function IR_Calculator()
%% (D)IR MRI TI/T1++ calculations
%   - In mode 1 & 2, find TI time(s) that null given T1(s) in DIR & 1IR (FLAIR/STIR/etc).
%   - In mode 3, calculate/plot relative signal strength for WM & WML as a function of TI2 in a DIR sequence,
%       - find TI2 and TI1(TI2) such that there is no T1 weighting between WM & WML while nulling CSF,
%       - and find a T1 that's nulled at the resulting TIs, corresponding to a fictive tissue.
%   - Also plot the T2 decay at readout (in a separate file), to determine the optimal TE.
% 
%  NOTE:
%   - TI times are defined vendor-style here! That is, TI1 is the full time and TI2 only the time to readout.
%       - So TI1 = time_TI1, and TI2 = time_TI1 - time_TI2. This makes the signal formulae simpler.
%   - Pulse durations are calculated into Ti times, but not T2prep except when reporting Siemens times. Vendor specific.
%   - Relaxation times vary widely in literature. We have used values from GE's scanner implementation, and other sources.
%       - Lalande et al, MRI 2016: https://www.sciencedirect.com/science/article/pii/S0730725X16301266
%   - There's a good explanation of simple inversion mathematics at http://xrayphysics.com/contrast.html
%   - About (single) incomplete inversion, see Miho Kita et al, MRI 31(9) 2013 p 1631–1639.
%       - Kita et al, MRI 2013: https://www.sciencedirect.com/science/article/pii/S0730725X13002312
%       - https://www.seichokai.or.jp/fuchu/dept1603.php (partially in Japanese; contains formulae)
%       - OliG: 95-99% is realistic, depending on the pulses. Nonselective IR is better. (HypSec ~98%?)
%   - Simulated TI_FLAIR/TI1 were lower than the ones used by Siemens/GE, typically by 150 ms @ TR 5 s w/ T2prep
%       - Intentional, accepting a little CSF signal to gain WML signal and contrast? Seems not. CSF is really nulled.
%       - Simulating 5-6% residual CSF signal does give +20% signal/contrast for WML! Warrants further study.
%       - Maybe a somewhat lower effective (or "apparent") TE could be a further tweak to improve SNR/CNR.
%       - For Philips, the simulation seems in agreement with their values.
%       - Upon checking GE 3D-FSE source, turns out the reason was compensations for CSF Mz @ end of readout.
%       - GE uses a Mz_end "Mz at the end of XETA/CUBE refocusing train, assuming 130-deg max flip", only for FLAIR
%       - GE seemingly uses artificially high T1csf for DIR in lieu of explicit Mz_end calculations?
% 
%  TOFIX/WIP:
%   - Set TR/T_Inv etc on vendor change; back to former when leaving. In/Out.
%   - When setting TE in WM-DIR, RelC got the wrong sign on the Results printout.
%   - Phi adds T2P time to TI?! Add Phi vendor with this.
%   - Add relax. times based on Neema, Hanson etc reporting increased T1/T2 in NAWM/NAGM in serious MS patients!?
%       - Check out a range of these to see the effect on parameters and signal/contrast.
%   - For 7T w/o T2Prep, no T1-nulled solution is found below TReff ~> 12 s. But there is one!? Find it?
%   - Correct SNR/CNR for field strength (the nominal increase is max. prop. to B0... or?)
%   - Philips center pulses around TIs, not starting there (so their t=0 is not ours!). Correct their times accordingly.
%   - SetRelaxTimes(), SetMode() and setT1n() are an interdependent mess! Merge at least the latter two?!?
% 
%  TODO:
%   - Move all setup fcns to a separate file (so one fcn calls the others in that file?)?
%   - Same with UI creation? Look to GUIDE for guidance?! The UI control values should all be in an accessible struct.
%   - A tool to calculate optimal TR w/ respect to CNR/t given the settings?
%   - If the GUI panel were wider, could have two fields per line?
%   - Better Bloch simulation of curves! Simply let Mz -> [Mz,Mt], and move Mz->Mt for excitation etc.
%       - Plot Mt with dotted lines?
%       - Do we have spoiling between refocus pulses? Do Phi use hard pulses with +/- phase (semi-balanced)?
%   - THS: T2 effect during adiabatic inversion!
%       - Compensate by adding a little T2prep during inv.? E.g., half the pulse duration? But that's just guesstimate.
%       - If we had [Mz,Mt] sim. and the adiabatic trajectory, we could get this more accurately.
%   - Calculate times/settings for all three major vendors, switchable with iS.Vnd
%   - Fix the T1null curves for T2prep. The TI values are right but the red/blue curves don't follow.
%   - Rewrite with handles structure ( use guihandles(), guidata() etc ) like GUIDE; check GUIDE first.
%   - Use less globals! Maybe also a struct/cell array for a set of properties/settings? Doesn't work for some vars?
%   - Include AtleB's ASL BS TI optimalization? Use a press button for it, as it takes some time to run.
%   - Implement incomplete spin lock during readout? TSE refocus pulses may be 120-150° or even 40-60° instead of 180°.
%       - So readout CPMG inversion efficiency is then -cos(FA)=50-87% using 120°. It's not (only) inversion!
%       - Merely simulating incomplete inversions doesn't lead to a change in behaviour, as spin lock is still strong.
%       - We'll need a more thorough Bloch simulation to account for all the spin components.
% 
%  DONE:
%   - Vendor setting, so sequence details and time calculations can be adjusted accordingly.
%   - Pulse durations. HypSec inv. pulses are typically ~15 ms, which is significant (GE uses pw_rf0 = 16ms)
%       - No T1 during pulses?! So subtract pulse dur. from TReff
%       - Simulate the first pulse as starting at time = 0 (unlike Phi; correct for that in their times)
%       - In the mag plot, just let the z mag drop (sigmoidally?) to its new value during the pulse duration. -tan(t)?
%   - A setting for residual CSF percentage. Let's say we accept 5% (DIR CSF signal = -0.05) to get better contrast...?
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
%       - Without reclaiming (drive)? What is the case in our sequences?
%       - Add a GUI for effective FA (90° for full excitation, 0° if refocused). Maybe also add drive?
%       - Retired the unused excitation FA UI again. It's commented out but still there should anyone need it.
%   - Implement incomplete inversion using the Inversion Efficiency = -cos(FA_inv)
%   - Added a generic S0 formula CalcS0() instead of CalcS0S() and CalcS0D().
%   - Made a T2-TE figure that can be updated from the GUI panel for the current signal/plot settings
%   - Correct CNR for sB0
%   - Simulate "T2 prep", in which shorter T2 times are dephased so their signal becomes nearly saturation recovery.
%       - Allow T2 entry? Crowds the UI. Alternatively, only emulate T2prep if times aren't set manually? Or ignore.
%   - Found and implemented T1/T2 values for WM/GM/CSF from GE 3dfse source code.
%
%  ONHOLD:
%   - Don't make signal equations generally symbolic rather than using functions, as symexpr are slower.
%   - Anonymous function handle for CalcTI1, for easier passing to other functions such as vpasolve()? Likely not.


%% INIT
% global debug; if isempty(debug); debug = 0; end % if
% if( debug ~= 0 )
%     clear all
% %     saveImg = true;                                 % Save figure as image in Figures dir? (Can use the menu instead)
% end % if
debugInfo  = 1;                                     % Show extra info for debugging in the command window?

global iPr iUI iS Ts iT iR                          % Program/UI/system/tissue/times/results data structs
global ircInit iUIBox IE                            % Globals that didn't thrive in structs or need short names

[pathstr,~,~]=fileparts(mfilename('fullpath'));     % [pathstring,name,extension]
addpath(genpath(pathstr));                          % Add all subfolders to search path
rmpath(genpath([pathstr filesep '.git']));          % Remove unnecessary search path
clear pathstr

% sSz = get(0,'ScreenSize');                          % Get the screen size (actually, in R2015b+ pixels are 1/96")
fSz = 10.80; % sSz(4)/100.00;                       % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])    % Default fig. position & size
wh = [ 1 1 1 ]; iPr.FigClr = 0.92*wh; iPr.Gray = 0.50*wh; % Colors of fig. background and grayed-out UI elements
% fig1 = get(groot, 'CurrentFigure');                 % Get handle of the current fig. ('gcf' creates one)
figName = 'MR-IR inversion time calculator';
iPr.Fig1 = findobj( 'Type', 'Figure',           ... % Find this script's figure (possibly among others)
        'Name', figName                         );
if isempty(iPr.Fig1)                                % Iff my figure does not yet exist... (was ~ishghandle)
    iPr.Fig1 = figure(                          ... % ...make a new figure!
        'Name', figName,                        ...
        'Color', iPr.FigClr,                    ...
        'KeyPressFcn', @keyPress_callback,      ...
        'NumberTitle', 'off',                   ... % Removes 'Figure #: ' from the window title
        'Position', [ 90 25 80 60 ]*fSz         );  % At pos. (90,25)% of ScreenSize, etc.
else
    figure(iPr.Fig1);                               % If my figure exists, make it the CurrentFigure
end % if

t = sym('t');                                       % Set the symbolic explicitly as it's used in nested functions
% TI2s = sym('TI2s'); TI1s(TI2s) = TI2s;              % Temporary assignment globalizes the sym expression
% DIR_S0(t) = ( 1 - 2*exp(-TI2/t) + 2*exp(-TI1/t) - exp(-TReff/t) ); DIR S0(T1) as a sym expression
% TI1s(TI2s) = T1*log(2/(2*exp(-TI2s/T1) + exp(-iS.TRef/T1) - 1));    % TI1(TI2) nulling T1_n1, using syms
% % Note: Do not replace CalcTI1(TI2,T1) with TI1s(t) above: It's much slower.

% Parameter presets for the first run (most can be changed in the UI)
if isempty( ircInit )
    ircInit = true;                                 % The program has been initialized

iS.Vnd  = ""        ;                               % Vendor name (for sequence adjustments)
iS.TR   = 6000.0    ; % ms                          % Repetition time, as float
iS.ETD  =  800.0    ; % ms                          % Readout time T_Ro = ETL * ESP (GE: Check CVs after Download)
iS.IEf  =  100      ; % %                           % Inversion efficiency in percent (are GE's SECH pulses ~98%?)
iS.IPD  =   16.0    ; % ms                          % Inversion pulse duration (GE/Phi used 16/17 ms HypSec in tests)
iS.T2pOn = true     ;                               % WIP: Turn this off when T1n is set manually (as T2 is then unknown)?
iS.T2p  =    0.0    ; % ms                          % T2 preparation time
iS.TOb  =   25.0    ; % ms  % WIP: Relate to iS.IPD % TI1 outboard (Siemens: Add at least 23 ms; +10? ms for T2p?)
iS.FAx  =   90      ; % °                           % (MEX = cos(FAx*pi/180) is MagZ after excitation with angle FAx)
iS.ESP  = 3.65; iS.IET = 0.50; % iS.ETL = 122;      % Echo Spacing, Echo Train Length and Inv.Eff. for TSE readout
                                                    % WIP: Just assuming incomplete inversion doesn't lead to anything!
                                                    % A proper Bloch simulation, at least w/ [ Mz, Mt ] is required.
                                                    % Actual flip angles are swept from about 10-170° to preserve Mz!
iR.RCS  =    0.0    ; % fraction                    % Residual CSF signal at readout. Partial nulling may increase CNR.
iS.Rew  =    0.0    ; %                             % TODO: Rewinder after readout, regaining some of the pre-readout Mz?
iUI.Show    = true  ;                               % Show the UI control panel
iPr.Mode    =    1  ;                               % Start in mode 1, then keep the run mode over reruns
iPr.TsSel   =    3  ;                               % Start with a value (CSF), then keep the selection over reruns
iPr.B0sel   =    4  ;                               % Start with GE 3T setup, w/ relaxation times from literature/source
SetRelaxTimes( iPr.B0sel )                          % Need to initialize relaxation times before mode.
SetMode( iPr.Mode )                                 % Initialize mode settings
end % if ircInit

%% MAIN
main();
function main()
    iS.TRef     = iS.TR   - iS.ETD - iS.IPD;        % Recalculate based on current parameter settings
    if ( iS.T2pOn )
        iS.TRef = iS.TRef - iS.T2p;
    end % if
    IE = iS.IEf/100; % IE2 = 1 + IE;                % MagZ to recover after incomplete inv. (was log(2/...) etc)
    iS.ETL = ceil(iS.ETD/iS.ESP);                   % Echo train length (# pulses) for TSE readout (ceil or floor best?)
    
    clf                                             % Clear the figure so we can change it
    switch iPr.Mode
        case 1                                      % Mode 1 (1IR tissue nulling)
            mainSIRNul( iT.Tn1 );
        case 2                                      % Mode 2 (DIR tissue nulling)
            iS.TRef = iS.TRef - iS.IPD;             % Subtract one more inversion pulse duration for DIR        WIP
            mainDIRNul( iT.Tn1, iT.Tn2 );
        case 3                                      % Mode 3 (T1W nulling)
            iS.TRef = iS.TRef - iS.IPD;             % Subtract one more inversion pulse duration for DIR        WIP
            Tn2 = mainT1WNul();                     % ( T1_CSF, T1_WM, T1_MS )  % TODO: Set iT.Tn2 here?
    end
    IR_CreateUI();
%     if saveImg, saveFigAsImg(),end                  % Automatically save the fig. as an image (or use the fig. menu)
end % main

% 1) Calculate TI in a 1IR/STIR/FLAIR sequence, nulling a specified T1 time
function mainSIRNul( T1z1 )
    iT.T1z = T1z1;                                  % What T1 are we concerned with nulling?
    Tis = Ts.T(:,iPr.TsSel);                        % What tissue is selected? (Needed for T2prep)
    E2  = CalcPE2( Tis );                           % Calculate E2 if using T2prep (=1 otherwise)
    if ( iS.Vnd == "GE" )
        Ts.Mze  = [ 0.0, 0.0, 0.045, 0.0, 0.0 ];    % GE "Mz @ end of CUBE train w/ 130° max FA"; only used for FLAIR
    else
        Ts.Mze  = [ 0.0, 0.0, 0.0  , 0.0, 0.0 ];    % [ WM, GM, CSF, MS, Fat ] Mz @ readout end (only used by GE or...?)
    end % if
    Mze = Ts.Mze( iPr.TsSel );                      % End Mz for the tissue (GE value for XETA/Cube refocusing train)
    iT.Ti1n = CalcTInT1( iT.T1z, E2, 1-Mze );       % Calc. TI to null T1z
    iT.Ti1n = uint16(round( iT.Ti1n + iS.IPD ));    % Add the inversion pulse duration; make value int for plot display
    
    TInts = [ 0, iT.Ti1n, iS.TR ];                  % SIR: Inversion is at t = 0, then readout at TI
    IR_T1_Magplot( TInts );                         % Plot Z magnetization (Mz) of several tissue T1s over time
    
    figTitle( "Single IR tissue nulling" )          % Add an informative title to the plot
    fStr = ['\\itTissue T1 time:\n\\rm'         ...
            '  T1_{n}   = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling IR:\n\\rm\\bf'  ...
            '  TI_{n}    = %4.i ms'             ];
    fVar = [ uint16(iT.T1z), uint16(iT.Ti1n)    ];
    if ( iS.Vnd == "Siemens" )
        fStr = [ fStr,                          ...
            '\n\\rm  TI_{Sie} = %4.i ms'        ];
        fVar = [ fVar, SieT(iT.Ti1n)            ];
    end % if
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    iPr.pStr = [                           '\n' ...
            '*******************************\n' ...
            '***    IR tissue nulling:   ***\n' ...
            '*******************************\n' ...
            '*  Tissue T1       = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  Nulling TI      = %4.i ms  *\n' ...
            '*******************************\n' ];
    iPr.pVar = [ iT.T1z, iT.Ti1n ];
    if ( iS.Vnd == "Siemens" )
        iPr.pStr = strcat( iPr.pStr, [          ...
            '***      For Siemens:       ***\n' ...
            '*******************************\n' ...
            '*  TI              = %4.i ms  *\n' ...
            '*******************************\n' ] );
        iPr.pVar = horzcat( iPr.pVar, SieT(iT.Ti1n) );  % Take into account T2p and some outboard time for Siemens
    end % if
%     if ( iS.Vnd == "GE" )
%         iPr.pStr = strcat( iPr.pStr, [          ...
%             '***      For GE (DV25):     ***\n' ...
%             '*******************************\n' ...
%             '*  TI              = %4.i ms  *\n' ...
%             '*******************************\n' ] );
%         iPr.pVar = horzcat( iPr.pVar, iT.Ti1n + 82 );   % Estimated from experiments by WN/ØBG 2018-09-26
%     end % if                                            % NOTE: We do better now, with the new GE corrections!

end % mainSIRNul

% 2) Calculate TI1+TI2 in a DIR sequence, nulling two specified T1 times (T1_n1 = T1_CSF by default)
function mainDIRNul( T1z1, T1z2 )
    iT.T1z = T1z2;                                  % What T1 are we concerned with nulling?
    Tis = [ Ts.T(:,3), Ts.T(:,iPr.TsSel) ];         % CSF and selected tissue (needed for T2prep)
    E2 = CalcPE2( Tis );                            % Calculate E2 if using T2prep (=1 otherwise)
    [ iT.Ti1n, iT.Ti2n ] = Find2TIn2T1([T1z1,T1z2],E2); % Inversion times that null 2 T1
    iT.Ti1n = uint16(round( iT.Ti1n + iS.IPD*2 ));  % Add the inversion pulse durations
    iT.Ti2n = uint16(round( iT.Ti2n + iS.IPD   ));  % Cast as int for plot display; ceil()/round() don't change syms!
    
    TInts = [ 0, (iT.Ti1n-iT.Ti2n), iT.Ti1n, iS.TR ];   % DIR: Vector of all time points (inversion, readout and end)
    IR_T1_Magplot( TInts );
%     PlotMagZ( TInts );                                  % Plot Z magnetization over time [returned end MagZ vector w/ Mzt = ]

    figTitle( "DIR dual tissue nulling" )           % Add an informative title to the plot
    
    fStr = ['\\itTissue T1 times:\n\\rm'        ...
            '  T1_{n1}  = %4.i ms\n'            ...
            '  T1_{n2}  = %4.i ms\n'            ...
            '  \n'                              ...
            '\\itTissue nulling DIR:\n\\rm\\bf' ...
            '  TI_{1}    = %4.i ms\n'           ...
            '  TI_{2}    = %4.i ms'             ];
    fVar = [ T1z1, T1z2, iT.Ti1n, iT.Ti2n ];
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    
    iPr.pStr = [                           '\n' ...
            '*******************************\n' ...
            '***   DIR tissue nulling:   ***\n' ...
            '*******************************\n' ...
            '*  T1 tissue 1     = %4.i ms  *\n' ...
            '*  T1 tissue 2     = %4.i ms  *\n' ...
            '*******************************\n' ...
            '*  TI1             = %4.i ms  *\n' ...
            '*  TI2             = %4.i ms  *\n' ...
            '*******************************\n' ];
%             '*  TI1 - TI2       = %4.i ms  *\n' ...
    iPr.pVar = [ T1z1, T1z2, iT.Ti1n, iT.Ti2n ];    % , (iT.Ti1n-iT.Ti2n)
    if ( iS.Vnd == "GE" )
        iPr.pStr = strcat( iPr.pStr, [          ...
            '***   Times for GE (DV25):  ***\n' ...
            '*******************************\n' ...
            '*  TI1             = %4.i ms  *\n' ... % Estimated from experiments by WN/ØBG 2018-09-26
            '*  TI2             = %4.i ms  *\n' ... % --"--
            '*******************************\n' ] );
        iPr.pVar = horzcat( iPr.pVar, [ uint16( iT.Ti1n + iS.TRef*0.038 - 80 ), iT.Ti2n + 12 ] );
    end % if

end % mainDIRNul

% 3) Calculate & plot MRI S0 of WM vs MS/WML in a DIR sequence, and TI1+TI2 nulling their T1 weighting
function T1z2 = mainT1WNul()                        % T1_nulled = T1WNul( T1_CSF, T1_WM, T1_MS )
    T1C = Ts.T(1,3); T1W = Ts.T(1,1); T1L = Ts.T(1,4);      % (CSF, WM, WML)
%     T2C = Ts.T(2,3); T2W = Ts.T(2,1); T2L = Ts.T(2,4);    % (--"--)
    iT.T1z = iT.Tn1; % Ts.T(1,iPr.TsSel);           % What T1 are we concerned with nulling? TODO: T1C? Allow T1_n1.
    T1C = iT.T1z;  % Override T1C with the set value (although we keep T2C from selection)
    Tis = Ts.T(:, [ 3 1 4 ] );                      % CSF, WM, WML (needed for T2prep)
    E2 = CalcPE2( Tis );                            % Calculate E2 if using T2prep (=1 otherwise)
    E2C = E2(1)    ; E2W = E2(2)    ; E2L = E2(3);  % (CSF, WM, WML)
    iT.Ti2max = T1C*log( (1+IE)/(1+IE*exp(-iS.TRef/T1C)) ); % Max TI2 nulling T1C, usually CSF (@ max TI1 = TReff)
%     iT.Ti2max = T1C*log(2/(1 + 2*exp(-iS.TR/T1C) - exp(-iS.TRef/T1C))); % Old Ti2max calculation without InvEff/T2P
    Ti2 = 1:1:iT.Ti2max; DUR = length(Ti2);         % Vector of TI2 values 1 ms apart (could've used linspace() here)
    Ti1 = zeros(DUR,1);                             % Vector of TI1 values for all TI2s; replace by TI1s(TI2s)?
    S0_DIR = zeros(DUR,2);                          % Array of two DIR signals@TE=0 for TI2s
    
%     TI1s(t) = CalcTI1(t,T1C);                       % TI1(TI2) calc.
    for i = 1:DUR
        Ti1(i)      = CalcTi1D(Ti2(i),T1C,E2C);         % TI1(TI2) calc. for T1_CSF
        S0_DIR(i,1) = CalcS0([Ti1(i),Ti2(i)],T1W,E2W);  % S0(TI2) calc. (was CalcS0D(TI1(i),TI2(i),T1W) )
        S0_DIR(i,2) = CalcS0([Ti1(i),Ti2(i)],T1L,E2L);  % --"--
    end
    
    iT.Ti2n = uint16( FindTi2nT1W(Tis(1,:),E2) );   % Find TI2 to null the DIR T1 weighting between two T1s (T1C, T1W, T1L)
    iT.Ti1n = uint16( CalcTi1D(iT.Ti2n,T1C,E2C)  ); % TI1(TI2) calc. by function (uint to avoid e+ notation in fig. text)
%     iT.Ti1n = uint16(TI1s( iT.Ti2n ));              % TI1(TI2) calc. by symbolic expr. (slower?!)
    
    % NOTE: GE operates with a guessed T2 when T1 is set manually. DV25: T2wm; RX27: T2gm. But is it used for T2prep?
    % WIP: For now, let's ignore T2prep w/ E=1. Judging from observation, that's what GE does too!?
    iT.Tn2 = uint16( FindT1null([iT.Ti1n,iT.Ti2n],1) ); % Second tissue nulled by DIR. UInt for display.
    T1z2 = iT.Tn2;                                  % This function outputs the fictive T1 nulled
    iS.newS = S0_DIR(iT.Ti2n,1);                    % WM signal at readout; = WML signal here
    
    hold on
    plot( Ti2, S0_DIR(:,1), 'r-',               ...
        'LineWidth', 1.0                        )   % Plot the signals calculated above
    plot( Ti2, S0_DIR(:,2), 'b-' ,              ...
        'LineWidth', 1.0                        )
%     S0W_t   = CalcS0([iT.Ti1n,iT.Ti2n],T1W,E2W);    % DEBUG: S0(Ti1,Ti2) calc.
%     S0L_t   = CalcS0([iT.Ti1n,iT.Ti2n],T1L,E2L);    % --"--
%     disp( "    S0W,      S0L      (Debug T1null)" )
%     disp( [ S0W_t, S0L_t ] )
%     plot(Ti2,Ti1);                                  % DEBUG: Plot the TI1 calculated for each TI2
%     S0_DIFF = abs(S0_DIR(:,1)-S0_DIR(:,2));         % DEBUG: Calculate signal difference between two T1s
%     [~,minind] = min(S0_DIFF);                      % Find the T1Wnull crossing graphically (as [minval,minind])
%     iT.Ti2n = Ti2(minind);
%     plot(Ti2,S0_DIFF);                              % DEBUG: Plot the WM/MS signal difference
    hold off
    
%     whitebg( [ 0.70 1.00 0.70 ] );                  % DEBUG: See white plot/GUI elements
%     whitebg( 'white' );                             % NOTE: Setting this also resets the figure color.
    figTitle( "T1-DIR (T1W nulling)" )                  % Add an informative title to the plot
    ylabel('Rel. MR signal (TE = 0 ms)', 'FontSize', 12 )
    xlabel('TI_2 inversion time (ms)',   'FontSize', 12 )
    xlm = [ 0 200*ceil( iT.Ti2max/200 ) ]; xlim(xlm);   % Round up to 200 on the x axis
    ylm = [ -1 1 ];                      ylim(ylm);     % Rel. S0 is in [-1,1]
    
    LegT = [ 1 4 ];                                     % Make a legend for tissue 1(WM) & 4(MS)...
    NT1s = length(LegT); LegS = strings(1,NT1s);
    for i = 1:NT1s, j = LegT(i);                        % ...showing the tissue short names and T1 times
        LegS(i) = sprintf('%s (T1 = %i ms)', Ts.Tag(1,j), Ts.T(1,j) );
    end
    legend( LegS(1), LegS(2),                   ... % Was: legend( 'WM', 'MS', ...
            'Location', 'NorthWest'             );  % 'Best' may conflict with UI and TextBox
    
    lxlm = [ iT.Ti2n iT.Ti2n ];                         % S0 crossing x (=TI2) value
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
    fVar = [ iT.Ti1n, iT.Ti2n, T1z2 ];              % Was: T1_WM, T1_MS, etc
    figTextBox( fStr, fVar );                       % Display an info text box on the figure

    iPr.pStr = [                           '\n' ...
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
    iPr.pVar = [ T1C, T1z2, iT.Ti1n, iT.Ti2n, (iT.Ti1n-iT.Ti2n) ];
    
end % mainT1WNul

%% FUNCTIONS
function FltOrSym = FoS( InNr )                     % Cast as double for, e.g., exp() unless it's used as a sym
    if      isa( InNr, 'sym' ), FltOrSym = InNr;
    else,                       FltOrSym = double(InNr);
    end % if
end % fcn

function SieTI = SieT( TI )                         % Add outboard time estimated for the Siemens sequence
    SieTI = TI + iS.TOb;
    if ( iS.T2pOn ) && ( iS.T2p > 0 )               % Add a little extra outboard (10 ms?) for T2 prep?
        SieTI = SieTI + iS.T2p;                     % Was + 10 for T2prep, but that's iffy...
    end % if
    SieTI = uint16(SieTI);
end % fcn

function E2 = CalcPE2( T1T2 )                       % Calculate E2 for T2prep for a (vector of) T2 time(s)
    if ( iS.T2pOn ) && ( size(T1T2,1) > 1 )         % Only use the T2prep formula if T2 is known
        E2 = exp( -iS.T2p./T1T2(2,:) );
    else
        E2 = ones(1,size(T1T2,2));
    end % if
    E2 = double(E2);
end % fcn

function S0_IR      = CalcS0( Ti, T1, E2 )          % Calculate generic IR relative S0 (def: S/M0 = E2*S0 @ TE=0 ms)
    Ti = FoS(Ti); T1 = FoS(T1);                     % Cast input as float unless it's a sym
    NI = length( Ti ); S0_IR = ones( 1,length(T1) );            % Ti() are defined as time to readout
    for i = 0:NI-2    % (was 0:NInv-1)
        S0_IR = S0_IR - (-IE)^i*(1+IE).*exp(-Ti(NI-i)./T1);     % Step backwards through the inversions
    end
    S0_IR = S0_IR - (-IE)^(NI-1)*(1+IE*E2).*exp(-Ti(1)./T1);    % The last step includes E2 for T2prep
    S0_IR = S0_IR - (-IE)^(NI-0)*E2.*exp(-iS.TRef./T1);         % Gives one signal value per given T1
end % fcn                                           % NOTE: Inv. pulse dur. isn't incorporated in these calculations

% function S0_IR      = CalcS0S( Ti, T1, E2 )         % Calculate single IR relative S0 (def: S0 = S/(E2*M0) @ TE=0)
%     Ti = FoS(Ti); T1 = FoS(T1);                     % Cast input as float unless it's a sym
%     S0_IR = 1 - (1+IE*E2).*exp(-Ti./T1) + IE*E2.*exp(-iS.TRef./T1);
% end % fcn

% function S0_DIR     = CalcS0D( Ti1, Ti2, T1, E2 )   % Calculate DIR relative S0 (def: S0 = S/(E2*M0) @ TE=0)
%     Ti1 = FoS(Ti1); Ti2 = FoS(Ti2); T1 = FoS(T1);   % Cast input as float unless it's a sym
%     S0_DIR = 1 - (1+IE).*exp(-Ti2./T1) + IE*(1+IE*E2).*exp(-Ti1./T1) - IE^2*E2.*exp(-iS.TRef./T1);
% end % fcn

function TI_SIR     = CalcTInT1( T1, E2, MC )       % Calc. TI nulling T1 in a single IR sequence w/ T2prep and Mz corr.
    T1 = FoS(T1);                                   % Cast input as float unless it's a sym
%     TI_SIR = log(2)*T1;                             % The simple case, TR>>T1, gives TIn = ln(2)*T1 ~ 0.69*T1
    TI_SIR = T1*log( (1+IE.*E2) ./ ( 1-iR.RCS + IE*MC.*E2.*exp(-iS.TRef./T1) ) );  % 1IR TI nulling T1
end % fcn

function TI1_DIR    = CalcTi1D( Ti2, T1, E2 )       % Calc. TI1 (new TI1, the old one was - TI2) nulling T1 given TI2
    Ti2 = FoS(Ti2); T1 = FoS(T1);                   % Cast input as float unless it's a sym
    TI1_DIR = T1.*log( IE*(1+IE*E2)./( -(1+iR.RCS) + (1+IE).*exp(-Ti2./T1) + IE^2*E2.*exp(-iS.TRef./T1) ) );  % (May be a sym)
end % fcn

function S0_nD      = CalcS0for0( Ti2, T1d, E2d )   % Calc. S0 from a DIR T1 found, combining CalcS0D+TI1D for sanity
    Ti2 = FoS(Ti2); T1d = FoS(T1d);                 % Cast input as float unless it's a sym (also, using 1x1 vars here)
    T1 = T1d(1); E2 = E2d(1);
    Ti1 = T1*log( IE*(1+IE*E2)/( -(1+iR.RCS) + (1+IE)*exp(-Ti2/T1) + IE^2*E2*exp(-iS.TRef/T1) ) ); % No .* ./ here: Only one T1
    T1 = T1d(2); E2 = E2d(2);
%     S0_nD  = CalcS0( [ Ti1, Ti2 ], T1(2), E2(2) );
    S0_nD = 1 - (1+IE)*exp(-Ti2/T1) + IE*(1+IE*E2)*exp(-Ti1/T1) - IE^2*E2*exp(-iS.TRef/T1);
end % fcn

function [ Ti1n, Ti2n ] = Find2TIn2T1( T1d, E2 )    % Find TI1/TI2 that null two T1 (e.g., CSF and WM) in DIR
    T1 = FoS(T1d);                                  % Cast input as float unless it's a sym
%     Ti2n = vpasolve( CalcS0( [ CalcTi1D(t,T1(1),E2(1)) , t ],T1(2),E2(2) ), t );
    Ti2n = vpasolve( CalcS0for0( t, T1, E2 ), t );
    Ti1n = CalcTi1D( Ti2n,T1(1),E2(1) );
end % fcn

function T1WnTi2    = FindTi2nT1W( T1_CWL, E2 )     % Find TI2 to null CSF and the DIR T1 weighting between WM and WML
    T1 = FoS(T1_CWL);                               % Cast input as float unless it's a sym
    T1WnTi2 = vpasolve( ( CalcS0for0( t, [ T1(1) T1(2) ], [ E2(1) E2(2) ] ) ...     % Solve S0DIR == 0 for this T1
                        - CalcS0for0( t, [ T1(1) T1(3) ], [ E2(1) E2(3) ] ) ), t ); % Return TI2 (calc. TI1(TI2) later)
end % fcn

function T1n_DIR    = FindT1null( Ti12, E2 )        % Find the second T1 time nulled by DIR
    Ti12 = FoS(Ti12);                               % Cast input as float unless it's a sym
%     T1n_DIR = vpasolve(CalcS0D(Ti1o,Ti2o,t), t)     % Solve S0DIR == 0 for T1
    T1n_DIR = vpasolve( CalcS0(Ti12,t,E2), t);      % Solve S0DIR == 0 for T1
end % fcn

function SetMode( mode )                            % Sets the run mode and resets the Nulling selection (iPr.TsSel)
    % WIP: Merge SetMode with SetT1n()?!
    oldMode = iPr.Mode;                             % WIP: If oldMode = mode, just skip to SetT1n!? But not on first run?
%     iT.Tn1 = Ts.T(1,3);                             % Default for all modes: Null out CSF as the longest T1
    switch mode
        case 1                                      % 1IR default: Null out CSF (=FLAIR)
            iT.Tn1 = Ts.T(1,3);
            SetT1n( 1, 3 );
        case 2                                      % DIR default: Null out WM as tissue 2
            if ( oldMode == 3 )
                SetT1n( 2, -1 );                    % A trick to notify SetT1n() that we've done T1 nulling
            else
                SetT1n( 2, 1 );
            end % if
        case 3                                      % T1W-DIR default: (Null out CSF as the longest T1)
    end % switch
    iPr.Mode = mode;                                % Make the mode setting global
    if ( mode ~= 3 ) && ( oldMode ~= 3 )
        SetRelaxTimes( iPr.B0sel )                  % Vendor implementation and mode may affect relaxation time settings
    end % if
end % fcn

function SetT1n( mode, T1s )                        % Sets the Nulling selection and T1_n#
    if ( T1s > 0 )
        switch mode
            case 1                                  % 1IR (only one tissue nulled)
                iT.Tn1 = Ts.T(1,T1s);
            otherwise                               % DIR (set second T1 to null; the first is CSF by default)
                if Ts.T(1,T1s) >= iT.Tn1, T1s=1; end % Don't null Tn1/CSF twice; use WM as tissue 2 instead then.
                iT.Tn2 = Ts.T(1,T1s);
        end % switch
    else
        T1s = -T1s;
    end % if
    iPr.TsSel = T1s;
    
    switch mode                                     % Reset which T1 times to plot, based on mode. TODO: Change this?!?
        case 1                                      % 1IR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4;
        case 2                                      % DIR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4;
        case 3                                      % T1W-DIR: Plot the signal for WM and WML
            i = [ 1 4 ];                            % NOTE: This is not implemented yet; Mode 3 uses its own plot
    end % switch
    Ts.Plt = Ts.T(:,i); Ts.Leg = Ts.Tag(1,i);           % Update the T1 times to plot and their label texts
end % fcn

function SetRelaxTimes( B0opt )
    T1_CSF = IR_SetRelaxTimes( B0opt ); % This is in a separate file to make this one shorter and things clearer
    
    iT.Tn1 = T1_CSF;                    % Default for all modes: Null out CSF as the longest T1
    switch iPr.Mode                     % When a scheme is selected, adjust the nulling selection accordingly.
        % TODO: Remove this switch; just call SetT1n( iPr.TsSel ) and let the function handle iPr.Mode!
        case 1                          % 1IR T1 nulling
            SetT1n( 1, iPr.TsSel );
        case 2                          % DIR T1 nulling
            SetT1n( 2, iPr.TsSel );
        case 3                          % T1W nulling
    end % switch
%     toggleT2prep( true );                           % T2prep is compatible with tissue selection
end % setRelaxTimes

%% PLOT

% Note: The PlotMagZ function was moved to another file since this file is so long

function figTitle( startStr )                       % Add an informative title to the plot
    txt = [ char(startStr) ' at ' iPr.B0tag ', TR_{eff} = ' num2str(iS.TRef) ' ms' ];
    if ( iS.Vnd ~= "" )                             % If a vendor is specified, show it
        txt = [ txt ' (for ' char(iS.Vnd) ')' ];
    end % if
    title(  txt,                                ...
            'FontSize', 16                      );
end % fcn

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
end % fcn

% function saveFigAsImg()
%     imgName = ['TReff_' int2str(iS.TRef) '_TI1_' int2str(iT.Ti1n) '_TI2_'  int2str(iT.Ti2n) ];
%     if exist('Figures','dir')~=7, mkdir('Figures'), end
%     cd Figures
%     % TODO: What's the difference between print to image and the imwrite command?
%     print('-djpeg', '-r300', imgName )
%     cd '..'
% end % fcn

%% USER INTERFACE

function IR_CreateUI()                                 % Display a button to show/hide the UI controls
%     uiX = 74.0*fSz; uiY = 50.0*fSz;                 % Default start position (in px) for masterUI controls
    eW = 0.065 ; eH = 0.065; eS = eH + 0.025;       % Button width/height/spc for masterUI controls (was [50 40 60] px)
    myAx = gca;                                     % A handle to the current graphics axes (plot)
    uiX = myAx.OuterPosition(3) - 0.015 - eW;       % Use the axis' outer x position (def. [ 0 0 1 1 ]) for UI pos.
    uiY = myAx.Position(4)      + 0.010     ;       % Use the axis' inner y position to determine UI pos.
    
    uicontrol( 'Style', 'pushbutton',           ... % For all modes:
        'String', 'Panel',                      ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[p] Show/hide UI',    ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-4*eS eW eH ],     ...
        'Callback', @masterUI_callback          );  % UI to show/hide the UI control panel
    
    switch iPr.Mode
        case 1                                      % 1IR T1 nulling
            b2txt = 'DIR';
            b4txt = '';
        case 2                                      % DIR T1 nulling
            b2txt = '(FLA)IR';
            b4txt = 'T1-null';
        case 3                                      % T1W nulling
            b2txt = '';
            b4txt = 'Back';
    end % switch iPr.Mode
    if ismember( iPr.Mode, [ 1 2 ] )                % If 1IR/DIR:
    uicontrol( 'Style', 'pushbutton',           ... % uiMode
        'String', b2txt,                        ... % 'togglebutton'
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[w] Switch mode',     ...
        'Units', 'normalized',                  ... % Make the UI resize automatically with the figure
        'Position', [ uiX uiY-0*eS eW eH ],     ...
        'Callback', {@cMode_callback,[ 2 1 0 ]} );  % UI to set run mode (for mode 1–3, switch to...)
    
    uicontrol( 'Style', 'pushbutton',           ... % If 1IR/DIR:
        'String', 'Plot T2',                    ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[s] T2 decay plot',   ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-1*eS eW eH ],     ...
        'Callback', @plotT2_callback            );  % UI to create T2 vs TE plot
    end % if
    
    if ismember( iPr.Mode, [ 2 3 ] )                % If DIR/T1W_n:
    uicontrol( 'Style', 'pushbutton',           ...
        'String', b4txt,                        ...
        'FontWeight', 'bold',                   ...
        'ToolTipString', '[t] Switch mode',     ...
        'Units', 'normalized',                  ...
        'Position', [ uiX uiY-2*eS eW eH ],     ...
        'Callback', {@cMode_callback,[ 0 3 2 ]} );  % UI to change run mode to T1W nulling plot
    end % if
    
    iUI.Vnd = uicontrol( 'Style', 'pushbutton', ... % ui Vnd
        'String'    , 'Vendor'              ,   ...
        'FontWeight', 'bold'                ,   ...
        'ToolTipString', '[v] Select vendor',   ...
        'Units'     , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-5*eS eW eH ],   ...
        'Callback'  , @vendor_callback          );  % UI to select vendor implementation
%     if ( iS.Vnd == "" )
        iUI.Vnd.String = "Vendor";
%     else
%         iUI.Vnd.String = iS.Vnd;
%     end % if
    
    iUI.Stp = uicontrol( 'Style', 'pushbutton', ... % ui Stp
        'String'    , 'Step TR'             ,   ...
        'FontWeight', 'bold'                ,   ...
        'ToolTipString', '[d] [DEBUG]'      ,   ...
        'Units'     , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-7*eS eW eH ],   ...
        'Callback'  , @debug_callback           );  % UI to debug by increasing TR in steps

    createUIPanel();
    switch iUI.Show
        case true ; iUIBox.Visible = 'on';
        case false; iUIBox.Visible = 'off';
    end % switch iUI.Show
end % IR_CreateUI

function createUIPanel()
    % Make UI control buttons for saving images instead of the variable, like with printVars?
    %      Probably not necessary, since you can save the image from the File menu.
    
    boxClr = iPr.FigClr;                            % UI panel background color (same as fig background)
%     uiX = 58.3*fSz; uiY = 54.5*fSz;                 % Default start position (in px) for UI controls
    eN = 11; eW = 130; eH = 30;                     % Num. of UI elements; default element width/height (were 120/30)
    mX = 10; mY =  6;                               % x/y margins in UI panel
    eW2 = eW/2 - 3; mX2 = mX + eW/2 + 3;            % Width/position of half-width elements
    iUI.Hpx = eN*eH + 2*mY; uiT = iUI.Hpx - mY + 3; % UI panel height in pixels; top position for elements
    iUI.Wpx = eW + 2*mX;                            % UI panel width in pixels
    
    iUIBox = uipanel( iPr.Fig1,                 ... % iUIBox
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
%     uistack( iUIBox, 'top' );                       % Bring the panel to the top of the graphic stack (not necessary)
    placeUiBox();                                   % Adjust the UI panel position to the figure window
    
    for i = 1:eN; eY = uiT - i*eH;                    % Generate each UI element on a new row in the panel
        switch i
        case  1
        uicontrol( 'Parent',    iUIBox,             ... % ui pT1s
            'Position', [ mX eY+0 eW2 24 ],         ...
            'Style',            'pushbutton',       ...
            'String',           'Info',             ...
            'ToolTipString',    '[f] Print info',   ... % 'Print T1 times'
            'Callback', @printInfo_callback         );  % UI to print settings to the command window
        uicontrol( 'Parent',    iUIBox,             ...
            'Position', [ mX2 eY+0 eW2 24 ],        ...
            'Style',            'pushbutton',       ...
            'String',           'Results',          ...
            'ToolTipString',    '[r] Print results',...
            'Callback', @printVars_callback         );  % UI to print results to the command window

        case 2
        iUI.B0T = uicontrol( 'Parent', iUIBox,      ... % ui B0[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Ref. values for T1s', ...
            'String', 'T1s:'                        );
        iUI.B0S = uicontrol( 'Parent', iUIBox,      ...
            'Position', [ mX+24 eY+0 eW-34 23 ],    ... % Note: Height setting won't change the actual height here
            'Style', 'popup',                       ...
            'String',   iPr.B0str,                  ...
            'Value',    iPr.B0sel,                  ...
            'Callback', @b0_sel_callback            );  % UI to set relaxation times based on B0 and ref.
        
        case  3
        uicontrol( 'Parent', iUIBox,                ... % ui TR[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', '[q] Repetition time', ...
            'String', 'TR:                          ms');
        iUI.TRS = uicontrol( 'Parent', iUIBox,      ...
            'Position', [ mX+40 eY+0 50 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'ToolTipString', 'TR = TR_ef + T_Ro + T_T2p', ...
            'String',   iS.TR,                      ...
            'Value',    iS.TR,                      ...
            'Callback', @tr_set_callback            );  % UI to set TR

        case  4
        uicontrol( 'Parent', iUIBox,                ... % ui TRo[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'TSE Readout time',    ...
            'String', 'T_Ro:                      ms');
        uicontrol( 'Parent', iUIBox,                ...
            'Position', [ mX+40 eY+0 50 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'ToolTipString', 'ETD = ETL * ES',      ...
            'String',   iS.ETD,                     ...
            'Value',    iS.ETD,                     ...
            'Callback', @tro_set_callback           );  % UI to set T_Read

        case  5
        uicontrol( 'Parent', iUIBox,                ... % ui IE[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Inversion efficiency',...
            'String', 'Inv.eff.:                   %');
        uicontrol( 'Parent', iUIBox,                ...
            'Position', [ mX+40 eY+0 50 20 ],       ...
            'Style', 'edit',                        ...
            'String',   iS.IEf,                     ...
            'Value',    iS.IEf,                     ...
            'ToolTipString', 'IE = -cos(FA_inv)',   ...
            'Callback', @ie_set_callback            );  % UI to set inversion efficiency

        case  6
        uicontrol( 'Parent', iUIBox,                ... % ui IPD[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Inv. pulse dur.',     ...
            'String', 'T_Inv.:                     ms');
        uicontrol( 'Parent', iUIBox,                ...
            'Position', [ mX+40 eY+0 50 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 50000, 'Max', 50000,             ...
            'ToolTipString', 'Inv. pulse dur.',     ...
            'String',   iS.IPD,                     ...
            'Value',    iS.IPD,                     ...
            'Callback', @ipd_set_callback           );  % UI to set inversion pulse duration

        case  7
        iUI.T2pT = uicontrol( 'Parent', iUIBox,     ... % ui T2p[T|S]
            'Position', [ mX eY+1 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'T2 preparation',      ...
            'String', 'T2p:                        ms');
%             'Enable',        'Inactive',            ... % TODO: Make the text clickable, reflecting state?
%             'ButtonDownFcn', 'toggleT2prep( -1 )'   );  % Problem: Fcn will be run in main data space so no local fns.
        iUI.T2pS = uicontrol( 'Parent', iUIBox,     ...
            'Position', [ mX+40 eY+0 50 20 ],       ...
            'Style', 'edit',                        ...
            'Min', 10000, 'Max', 10000,             ...
            'ToolTipString', 'E2p = exp(-T2p/T2)',  ...
            'String',   iS.T2p,                     ...
            'Value',    iS.T2p,                     ...
            'Callback', @t2p_set_callback           );  % UI to set T2prep duration (in ms)
%         iUI.T2pS.ForegroundColor  = iPr.Gray;           % Gray out the T2prep text for now since it's WIP
        iUI.T2pB = uicontrol( 'Parent', iUIBox,     ... % ui T2p[T|S]
            'Position', [ mX+23 eY+3 14 14 ],       ...
            'Style', 'radio',                       ...
            'ToolTipString', '[2] T2 prep. On/Off', ...
            'Callback', @t2p_toggle_callback        );  % UI to set T2prep duration (in ms)
    if ( iS.T2pOn )
        iUI.T2pS.Enable         = 'on';
        iUI.T2pB.Value          = 1;
    else
        iUI.T2pS.Enable         = 'off';
        iUI.T2pB.Value          = 0;
    end % if

        case  8
        iUI.Tn1T = uicontrol( 'Parent', iUIBox,     ... % ui Tn1[T|S]
            'Position',     [ mX eY+1 eW 16 ],      ...
            'Style',        'text',                 ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Longer T1 to null',   ...
            'String', 'T1_n1:                    ms'); % Note: UIControl text can't support Tex/LaTeX
        uicontrol(              'Parent', iUIBox,       ...
            'Position',         [ mX+40 eY+0 50 20 ],   ...
            'Style',            'edit',                 ...
            'Min', 10000,       'Max', 10000,           ...
            'String',           iT.Tn1,                  ...
            'Value',            iT.Tn1,                  ...
            'Callback',         @t1n1_set_callback      );  % UI to set Tn1 (longer T1, typically CSF)

        case  9
            if ismember( iPr.Mode, [ 2 3 ] )                % DIR T1 nulling and T1-DIR null two T1
            iUI.Tn2T = uicontrol( 'Parent', iUIBox,     ... % ui Tn2[T|S]
                'Position',     [ mX eY+1 eW 16 ],      ...
                'Style',        'text',                 ...
                'BackgroundColor', boxClr,              ...
                'HorizontalAlignment', 'Left',          ...
                'ToolTipString', 'Shorter T1 to null',  ...
                'String', 'T1_n2:                    ms');
            iUI.Tn2S = uicontrol( 'Parent', iUIBox,     ...
                'Position',     [ mX+40 eY+0 50 20 ],   ...
                'Style',        'edit',                 ...
                'Min', 10000,   'Max', 10000,           ...
                'String',       iT.Tn2,                  ...
                'Value',        iT.Tn2,                  ...
                'Callback',     @t1n2_set_callback      );  % UI to set Tn2 (shorter T1; GM/WM/etc)
            end % if

        case 10
            iUI.STnT = uicontrol( 'Parent', iUIBox,     ... % ui STn[T|S]
                'Position',     [ mX eY+1 eW 16 ],      ...
                'Style',        'text',                 ...
                'BackgroundColor',      boxClr,         ...
                'HorizontalAlignment',  'Left',         ...
                'ToolTipString', 'Ref. tissue to null', ...
                'String',       'Null T1:'              );
            if ismember( iPr.Mode, [ 1 2 ] )                % 1IR and DIR T1 nulling use tissue selector
            iUI.STnS = uicontrol( 'Parent', iUIBox,     ...
                'Position',     [ mX+40 eY+0 eW-50 22 ],...
                'Style',        'popup',                ...
                'String',       Ts.Tag(2,:),            ...
                'Value',        iPr.TsSel,              ...
                'Callback',     @t1_sel_callback        );  % UI to set T1_n(2) by tissue
            if iT.T1z ~= Ts.T(1,iPr.TsSel)
                iUI.STnS.ForegroundColor  = iPr.Gray;       % Gray out the tissue selector if entering values manually
            end % if iT.T1z
            else
                iUI.STnT.String = 'Rel. S0:                     %';
                if ~isfield(iS,'oldS'); iS.oldS = iS.newS; end % if     % For use in the relSNR uicontrol box
                if ~isfield(iS,'oldB'); iS.oldB = iS.B0  ; end % if     % (w/ normal global var., was if isempty)
                relS0 = uint16(100*(iS.newS/iS.oldS)*(iS.B0/iS.oldB));  % Rel.S0 is prop. to CNR, and rel.SNR here
                if ( debugInfo )
                    fprintf("   "); fprintf("%-10s", ["new S0L","old S0L","ratio","(rel. units)"] ); fprintf("\n");
                    disp( [ iS.newS, iS.oldS, iS.newS/iS.oldS ] )
                end % if debug
                iS.oldS = iS.newS;
                uicontrol(      'Parent', iUIBox,       ... % Text box to show rel. S0 (instead of tissue selector)
                'Position',     [ mX+40 eY+0 60 20 ],   ...
                'Style',        'edit',                 ...
                'Min', 10000,   'Max', 10000,           ...
                'String',       relS0,                  ...
                'Enable',       'Inactive'              );
            end % if

        case  11
        uicontrol( 'Parent', iUIBox,                ... % ui RS[T|S]
            'Position', [ mX eY+2 eW 16 ],          ...
            'Style', 'text',                        ...
            'BackgroundColor', boxClr,              ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString', 'Residual signal',     ...
            'String', 'Res.S:                 %'    );
        uicontrol( 'Parent', iUIBox,                ...
            'Position', [ mX+40 eY+0 40 20 ],       ...
            'Style', 'edit',                        ...
            'String',   100*iR.RCS,                 ...
            'Value',    100*iR.RCS,                 ...
            'ToolTipString', 'May improve CNR?',    ...
            'Callback', @rs_set_callback            );  % UI to set desired residual CSF signal

%         case  11
%         uicontrol( 'Parent', iUIBox,                ... % ui FAex[T|S]
%             'Position', [ mX eY+1 eW 16 ],          ...
%             'Style', 'text',                        ...
%             'BackgroundColor', boxClr,              ...
%             'HorizontalAlignment', 'Left',          ...
%             'ToolTipString', 'Excitation pulse FA', ...
%             'Visible',          'off',              ... % Disable this for now (not needed)
%             'String', 'FAex:                         °'); % Note: UIControl text can't support Tex/LaTeX
%         iUI.FAexS = uicontrol( 'Parent', iUIBox,    ...
%             'Position', [ mX+40 eY+0 60 20 ],       ...
%             'Style', 'edit',                        ...
%             'Min', 10000, 'Max', 10000,             ...
%             'String',   iS.FAx,                     ...
%             'Value',    iS.FAx,                     ...
%             'Visible',       'off',                 ... % Disable this for now (not needed)
%             'Callback', @fa_set_callback            );  % UI to set FA (in degrees)
        

% TODO: Make a button to show a list box for selecting what to plot, in a separate GUI panel?
%         case 11
%             iUI.ShTT = uicontrol( 'Parent', iUIBox,     ... % ui ShT[T|S]
%                 'Position', [ mX eY+1 eW 16 ],          ...
%                 'Style', 'text',                        ...
%                 'BackgroundColor', boxClr,              ...
%                 'HorizontalAlignment', 'Left',          ...
%                 'ToolTipString', 'Tissues to show', ...
%                 'String', 'Showing:'                    );
%             iUI.ShTS = uicontrol( 'Parent', iUIBox,     ...
%                 'Position', [ mX+40 eY-60 eW-40 82 ],   ...
%                 'Style', 'listbox',                     ...
%                 'Min',0,'Max',1,                        ... % trick to select multiple listbox items
%                 'String',   Ts.Tag(1,:),                  ...
%                 'Value',    Ts.Leg,                      ...
%                 'Callback', @t_disp_callback            );  % UI to set which tissues to display/calculate

         end % switch i
    end % for i
    
    switch iPr.Mode
        case 1                                      % 1IR T1 nulling
        iUI.Tn1T.String  = 'T1_n:                      ms';  % There is only one T1 to null in 1IR
        iUI.STnT.Position = iUI.STnT.Position + [ 0 eH 0 0 ];   % Move the tissue selector one line up
        iUI.STnS.Position = iUI.STnS.Position + [ 0 eH 0 0 ];   % --"--
%         iUI.Tn2S.Visible = 'off';                   % Hide the tissue 2 edit text/box (ui###T/S, if created)
        case 2                                      % DIR T1 nulling
%         iUI.B0S.ForegroundColor = iPr.Gray;         % Gray out the B0 select text/box (ui###T/S, if created)
        case 3                                      % T1W nulling
        iUI.Tn2S.Enable           = 'inactive';     % The tissue 2 selector is not editable
        iUI.Tn2S.ForegroundColor  = iPr.Gray;       % Gray out the tissue 2 edit text (it's set automatically)
%         iUI.FAexS.Enable          = 'off';          % Gray out the FA edit box (FA is irrelevant in this mode)
    end % switch iPr.Mode

end % createUIPanel

function val = str2val( str )                       % Function to process numeric input for UI control callback fcns
    val = str2double( str );
    if ( val ~= val )                               % test for NaN from malformed input
        val = 0;
    end % if
end % fcn

function masterUI_callback(~,~)                     % UI that shows/hides the UI control panel
    iUI.Show = ~iUI.Show;
    switch iUI.Show
        case true ; iUIBox.Visible = 'on';
        case false; iUIBox.Visible = 'off';
    end
end % fcn

function cMode_callback(~,~,toMode)                 % Switch between calculation modes
    switch iPr.Mode
        case 1                                      % 1IR T1 nulling    ->|
            SetMode( toMode(1) );
        case 2                                      % DIR T1 nulling    |<-
            SetMode( toMode(2) );
        case 3                                      % T1W nulling       <<-
            SetMode( toMode(3) );
    end % switch                                    % Note that the figure title shows the active mode
    main();
end % fcn

function vendor_callback(~,~)                       % Switch between vendor implementations
    switch iS.Vnd
        case ""                                     % ->
            iS.Vnd = "GE";
            iS.T2pOld = iS.T2p;                     % Store the current T2 prep time
            iS.T2p = 200.0;                         % GE uses this value consistently
        case "GE"                                   % ->
            iS.Vnd = "Siemens";
            iS.T2p = 170.0;                         % Sie uses this value consistently
%         case "GE"                                   % ->
%             iS.Vnd = "Philips";                     % NOTE: No special corrections exist for Philips at this point
%             iS.T2p = 125.0;                         % Phi uses this value as default
        case "Siemens"                              % <<-
            iS.Vnd = "";
            iS.T2p = iS.T2pOld;                     % Retrieve the previous T2 prep time
    end % switch                                    % Note that the figure title shows the active vendor
    SetRelaxTimes( iPr.B0sel )                      % Vendor implementation and mode may affect relaxation time settings
    main();
end % fcn

function plotT2_callback(~,~)                       % Run T2 plot (in separate script)
    iT2s        =   struct(                     ...
        'Fig1'  ,   iPr.Fig1                ,   ...
        'S0'    ,   iR.S0(1:4)              ,   ...
        'T2'    ,   Ts.T(2,:)               ,   ...
        'TsTag' ,   Ts.Tag(1,1:4)           ,   ...
        'T2tis' ,   [ 4 1 2 ]               );      % Tissues to T2 plot: WML, WM, GM
    iPr.Fig2 = IR_T2_Subplot( iT2s, iS );
    % iR.S0(1:4), Ts.T(2,:), Ts.Tag(1,1:4), [ 4 1 2 ], iS )    % S0/T2/tags for a set of tissues to plot
end % fcn

function b0_sel_callback(src,~)                     % UI to set rel. times based on B0/ref.
    SetRelaxTimes( src.Value );                     % Was get(src,'Value') pre-R2014b
    main();
end % fcn

function tr_set_callback(src,~)                     % UI to set TR
    iS.TR = str2val( src.String );                  % str2double( src.String );
    Lim = iS.ETD + 5; if iS.TR < Lim, iS.TR = Lim; end  % TR > T_Read
    main();
end % fcn

function debug_callback(src,~)                     % UI to step TR (DEBUG)
    iS.TR = iS.TR + 1000;
    main();
end % fcn

function tro_set_callback(src,~)                    % UI to set T_Read
    iS.ETD = str2val( src.String );
    Lim = iS.TR - 5; if iS.ETD > Lim, iS.ETD = Lim; end % TR > T_Read
    main();
end % fcn

function ipd_set_callback(src,~)                    % UI to set Inv. pulse dur.
    Lim = (iS.TRef + iS.IPD) - 1000;
    if ( iPr.Mode ~= 1 )
        Lim = double( iT.Ti2n ) - 50;               % TODO: The second TI needs recalculation!
    end
    iS.IPD = str2val( src.String );
    if iS.IPD >= Lim, iS.IPD = Lim; end
    if iS.IPD <  0  , iS.IPD = 0  ; end
    main();
end % fcn

function t1n1_set_callback(src,~)                   % UI to set Tn1 (longer T1)
    iT.Tn1 = str2val( src.String );
    if iPr.Mode > 1                                 % Using multiple inversions...
        Lim = iT.Tn2 + 1;
        if iT.Tn1 < Lim, iT.Tn1 = Lim; end          % Tn1 is the longest T1 time to null
    end
%     toggleT2prep( false );                          % T2prep is incompatible with setting T1 manually (no known T2)
    main();
end % fcn

function t1n2_set_callback(src,~)                   % UI to set Tn2 (shorter T1)
    iT.Tn2 = str2val( src.String );
    Lim = iT.Tn1 - 1;
    if iT.Tn2 > Lim, iT.Tn2 = Lim; end              % Tn1 is the longest T1 time to null
%     toggleT2prep( false );                          % T2prep is incompatible with setting T1 manually (no known T2)
    main();
end % fcn

function t1_sel_callback(src,~)                     % UI to select a tissue T1 to null
    SetT1n( iPr.Mode, src.Value );
%     toggleT2prep( true );                           % T2prep is compatible with tissue selection
    main();
end % fcn

% function fa_set_callback(src,~)                     % UI to select the flip angle (in degrees)
%     iS.FAx = str2val( src.String );
%     main();
% end % fcn

function t2p_set_callback(src,~)                    % UI to select the T2 preparation time (in ms)
    iS.T2p = str2val( src.String );
    main();
end % fcn

function ie_set_callback(src,~)                     % UI to select the inversion efficiency (in percent)
    iS.IEf = str2val( src.String );
    if iS.IEf <= 0, iS.IEf = 1; end                 % Zero/negative Inv.Eff. crashes the program
    main();
end % fcn

function rs_set_callback(src,~)                     % UI to select the desired residual CSF signal (in percent)
    iR.RCS = str2val( src.String )/100;
    Lim = 0.66;                                     % Too high residual signal crashes the program
    if iR.RCS >  Lim, iR.RCS =  Lim; end
    if iR.RCS < -Lim, iR.RCS = -Lim; end
    main();
end % fcn

function printVars_callback(~,~)                    % UI callback that prints results to the command window
    fprintf( iPr.pStr, iPr.pVar );
    if debugInfo && ismember( iPr.Mode, [ 1 2 ] )   % Debug info for plot
%         len = length(iR.MzIt);
%         fprintf('\nMagZ(TR) over the first %i repetitions (starting at 1):\n', len);
%         disp( iR.MzIt );                            % Show the iterative Z magnetization over the first TRs
%         disp( iR.MzIt(:,2:its+1) );
        fStr = [ '\nMagZ and Signal(0) by T1 at'    ...
            ' readout with flip angle %3.1f°'       ...
            ' [S0=%1.2f*MagZ(TI)]:\n'               ];
        fprintf(fStr, round(iS.FAx,1,'significant'), sin(iS.FAx*pi/180));   % Heading
        fprintf("\t "); fprintf("%-10s", Ts.Leg); fprintf("\n");            % Tissue legends
        if ( iS.FAx == 90 )                         % Report resulting S0 iff angle less than 90°
            fStr =             iR.S0'   ;
        else
            fStr = [ iR.MzTI'; iR.S0' ] ;
        end % if
        disp ( fStr );                              % Mz/S (after FA° pulse) at readout
    end % if
end % fcn

function printInfo_callback(~,~)                    % UI callback that prints settings to the command window
    fStr = [                               '\n' ...
            '*******************************\n' ... % fprintf( fStr, iPr.B0str(iPr.B0sel), T1set );
            '***    Relaxation times     ***\n' ...
            '***    for %-15s',       '  ***\n' ...
            '*******************************\n' ... % [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ]
            '*  T1 %-8s', '     = %+4s ms  *\n' ... % Used %4.i when printing T1s directly
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*  T1 %-8s', '     = %+4s ms  *\n' ...
            '*******************************\n' ...
            '*  T2 %-8s', '     = %+4s ms  *\n' ...
            '*  T2 %-8s', '     = %+4s ms  *\n' ...
            '*  T2 %-8s', '     = %+4s ms  *\n' ...
            '*  T2 %-8s', '     = %+4s ms  *\n' ...
            '*  T2 %-8s', '     = %+4s ms  *\n' ...
            '*******************************\n' ...
            '*         Settings:           *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*  %-11s',   '     = %+4s %%   *\n'...
            '*  %-11s',   '     = %+4s ms  *\n' ...
            '*******************************\n' ]; 
    T1info = [ Ts.Tag(2,1:5); num2cell(Ts.T(1,1:5)) ]; % Use this trick to weave strings and numbers into a string array
    T2info = [ Ts.Tag(2,1:5); num2cell(Ts.T(2,1:5)) ];
    if ( iS.T2pOn )
        T2prep = iS.T2p;
    else
        T2prep = "--";
    end % if
    T1pars = [  [  "TR    ", "T_Readout", "Inv.Eff.", "T_T2prep"];    ...
        num2cell([  iS.TR  ,  iS.ETD    ,  iS.IEf   ,  T2prep   ])    ];
    fprintf( fStr, iPr.B0str(iPr.B0sel), T1info{:}, T2info{:}, T1pars{:} );
end % fcn

function keyPress_callback(~,evt)
    switch evt.Key
        case 'p'                                    % Show/hide UI panel
            masterUI_callback([],[])                % Calling with dummy values
        case 's'                                    % Plot T2 decay
            plotT2_callback([],[])
        case 'v'                                    % Switch vendor implementation
            vendor_callback([],[])
        case 'd'                                    % Step TR (DEBUG)
            debug_callback([],[])
        case 't'                                    % Switch mode
            switch iPr.Mode
                case 1
                    cMode_callback([],[],[ 2 1 0 ])
                case 2                              % T1-nulling, plot T2 decay, focus back on this fig.
                    cMode_callback([],[],[ 0 3 2 ])
                    cMode_callback([],[],[ 0 3 2 ])
                    plotT2_callback([],[])
                    figure( iPr.Fig1 );
                case 3
                    cMode_callback([],[],[ 0 3 2 ])
            end % switch Mode
        case 'f'                                    % Print info
            printInfo_callback([],[])
        case 'r'                                    % Print results
            printVars_callback([],[])
        case 'w'                                    % Switch mode
            if ismember( iPr.Mode, [ 1 2 ] )
                cMode_callback([],[],[ 2 1 0 ])
            end % if
        case 'a'                                    % Change active figure
            if isfield( iPr, 'Fig2' )               % ishandle( iPr.Fig2 )
                figure( iPr.Fig2    );
            end % if
        case 'q'                                    % Focus on TR setting uicontrol
            uicontrol( iUI.TRS );
        case '2'                                    % Toggle T2prep on/off
            t2p_toggle_callback;
%         case 'c'                                    % Change focus to this figure (WIP)
%             figure( iPr.Fig1 );
    end % switch
end % fcn

function placeUiBox(~,~)                            % Function to resize the UI panel (automatically called)
    if ishandle(iUIBox)                             % (The first time this fn is called, the panel doesn't exist yet)
        pxFig = getpixelposition(iPr.Fig1);         % Figure size in px, for manual normalizing
        myAx = gca;                                 % A handle to the current graphics axes/plot
        uiW = iUI.Wpx / pxFig(3);                   % Normalized UI panel width
        uiH = iUI.Hpx / pxFig(4);                   % Normalized UI panel height
        uiX = myAx.Position(3) - uiW + 0.115;       % The axis' inner right x pos. (but why the "fudge factor"?)
        uiY = myAx.Position(4) - uiH + 0.090;       % The axis' inner upper y pos. (--"--)
        iUIBox.Position = [ uiX uiY uiW uiH ];
    end % if
end % fcn

function t2p_toggle_callback(~,~)                   % UI callback that toggles T2 prep on/off
    toggleT2prep(-1);
    main();
end % fcn

function toggleT2prep( set )
    if ( set == -1 )
        iS.T2pOn = ~iS.T2pOn;
    else
        iS.T2pOn = set;
    end % if set
%     if ( iS.T2pOn )                                 % Note: Can't change UI here, as the GUI is redrawn later.
%     else
%     end % if
end % fcn

end % script fn