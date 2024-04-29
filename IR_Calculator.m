function IR_Calculator()
%% (D)IR MRI TI/T1++ calculation and optimizing
%   - In mode 1 & 2, find TI time(s) that null given T1(s) in DIR & 1IR (FLAIR/STIR/etc).
%   - In mode 3, calculate/plot relative signal strength for WM & WML as a function of TI2 in a DIR sequence,
%       - find TI2 and TI1(TI2) such that there is no T1 weighting between WM & WML while nulling CSF,
%       - and find a T1 that's nulled at the resulting TIs, corresponding to a fictive tissue.
%   - Also plot the T2 decay at readout (in a separate file), to determine the optimal TE.
% 
%  NOTE:
%   - If changes are made to the files, make sure to `clear all` before running anew. To this end, see `DEBUG` below.
%   - TI times are defined vendor-style here: That is, TI1 is the full time and TI2 only the time to readout.
%       - So TI1 = time_TI1, and TI2 = time_TI1 - time_TI2. This makes the signal formulae simpler.
%   - Pulse durations are calculated into Ti times, but not T2prep except when reporting Siemens times. Vendor specific.
%   - Typical readout for our GE CUBE-FLAIR is ETL 161, ES 5.29 ms. 
%       - We've used 160*5.29 = 845 ? 850 ms as our standard T_Ro.
%   - Relaxation times vary widely in literature. We have used values from GE's scanner implementation, and other sources.
%       - Lalande et al, MRI 2016: https://www.sciencedirect.com/science/article/pii/S0730725X16301266
%   - There's a good explanation of simple inversion mathematics at https://xrayphysics.com/contrast.html
%   - About (single) incomplete inversion, see Miho Kita et al, MRI 31(9) 2013 p 1631–1639.
%       - Kita et al, MRI 2013: https://www.sciencedirect.com/science/article/pii/S0730725X13002312
%       - https://www.seichokai.or.jp/fuchu/dept1603.php (partially in Japanese; contains formulae; not working anymore?)
%       - https://seichokai.jp/fuchu/null_point_english/ (Miho Kita's online tool for calculating TIs for various sequences)
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
%   - For all fns in this and other files, pass the now global struct byref instead of making it global:
%       - iC = someFn( fnVars, iC )
%       - In newer MatLab, this is suggested as preferable to globals!
%   - Make the custom times a separate, non-version-controlled file for user interaction.
%       - IR_UserRelaxTimes.m should be at program root level.
%       - The program can generate it if it isn't there. And SetRelaxTimes should read from it.
%       - Make it one struct? Could have User1.iS.B0, User1.T1_WM etc.
%   - Set TR/T_Inv etc on vendor change; back to former when leaving. In/Out.
%   - When setting TE in WM-DIR, RelC got the wrong sign on the Results printout.
%   - Add relax. times based on Neema, Hanson etc reporting increased T1/T2 in NAWM/NAGM in serious MS patients!?
%       - Check out a range of these to see the effect on parameters and signal/contrast.
%   - For 7T w/o T2Prep, no T1-nulled solution is found below TReff ~> 12 s. But there is one!? Find it?
%   - Correct SNR/CNR for field strength (the nominal increase is max. prop. to B0... or?)
%   - Philips center pulses around TIs, not starting there (so their t=0 is not ours!). Correct their times accordingly.
%   - SetRelaxTimes(), SetMode() and setT1n() are an interdependent mess! Merge at least the latter two?!?
% 
%  TODO:
%   - For regression, could we generate a set of TrueT2 TI1/TI2/T1n based on a matrix of T1_WM/T1_WML/TR values?
%   - OliG, 2021-03-09: Can we simulate a range of B1 inhomogeneity? Affects inversion eff., or adiabatic enough? Readout!
%   - Instead of the set of globals, have one 'iC' global and then iC.UI_###, iC.Pr_### etc.
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
%   - Calculate times/settings for all three major vendors, switchable with iC.S.Vnd
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
%   - Refresh T1/T2 shortcut (press 'r') for when custom values were added in IR_SetRelaxTimes().
%   - Vendor setting, so sequence details and time calculations can be adjusted accordingly.
%       - Phi and Sie add T2P time to TI. GE doesn't.
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


%% DEBUG
% global debug; if isempty(debug); debug = 1; end % if    % If changes are made to the files, activate this part to ensure a clear data space.
% if( debug ~= 0 )
%     clear all                                       ;   % A saveImg var to save figure as image in Figures dir? Can use the fig. menu instead.
% end % if debug
debugInfo  = false                                  ;   % Show extra info debugging info in the command window?

%% INIT
global iC                                               % IR-Calc S_ystem, P_rogram, M_atter(Tissue), T_imes, R_esults data structs
global iC_Ini1 F1_Pan UI1                               % Globals that didn't thrive in structs; UI for Fig 1   % TODO: Move UI1 to MagPlot.m?

myPath = mfilename('fullpath')                      ;   % The path to this script, including its file name
[myPath,~,~] = fileparts( myPath )                  ;   % [pathString,name,extension]
cd(myPath)                                          ;
addpath(genpath(myPath))                            ;   % Add all subfolders to search path
rmpath(genpath([myPath filesep '.git']))            ;   % Remove unnecessary search path (the .git folder)
clear myPath

% sSz = get(0,'ScreenSize')                         ;   % Get the screen size (actually, in R2015b+ pixels are 1/96")
fSz = 10.80; % sSz(4)/100.00;                           % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])	  % Default fig. position & size
wh = [ 1 1 1 ]; iC.P.FigClr = 0.92*wh; iC.P.Gray = 0.50*wh;	% Colors of fig. background and grayed-out UI elements
% fig1 = get(groot, 'CurrentFigure')                ;   % Get handle of the current fig. ('gcf' creates one)
figName = 'MR-IR inversion time calculator'     ;
iC.P.Fig1 = findobj( 'Type' , 'Figure'          ,   ... % Find this script's figure (possibly among others)
                    'Name'  , figName           )   ;
if isempty(iC.P.Fig1)                                   % Iff my figure does not yet exist... (was ~ishghandle)
    iC.P.Fig1 = figure(                             ... % ...make a new figure!
        'Name'          , figName               ,   ...
        'Color'         , iC.P.FigClr           ,   ...
        'NumberTitle'   , 'off'                 ,   ... % Removes 'Figure #: ' from the window title
        'Position'      , [ 90 25 80 60 ]*fSz   ,   ... % At pos. (90,25)% of ScreenSize, etc.
        'KeyPressFcn'   , @keyPress_callback    )   ;
else
    figure(iC.P.Fig1)                               ;   % If my figure already exists, make it the CurrentFigure
end % if

t = sym('t');                                           % Set the symbolic explicitly as it's used in nested functions
% TI2s = sym('TI2s'); TI1s(TI2s) = TI2s;                % Temporary assignment globalizes the sym expression
% DIR_S0(t) = ( 1 - 2*exp(-TI2/t) + 2*exp(-TI1/t) - exp(-TReff/t) ) ;   % DIR S0(T1) as a sym expression
% TI1s(TI2s) = T1*log(2/(2*exp(-TI2s/T1) + exp(-iC.S.TRef/T1) - 1)) ;   % TI1(TI2) nulling T1_n1, using syms
% % Note: Do not replace CalcTI1(TI2,T1) with TI1s(t) above: It's much slower.


if isempty( iC_Ini1 )                                   % Parameter presets for the first run (most can be changed in the UI)
                                                        % NB: Use `clear all` before running to ensure this part is run
    iC.S.Vnd    = ""        ;                           % Vendor name (for sequence adjustments)
    iC.S.TR     = 6000.0    ; % ms                      % Repetition time, as float
    iC.S.ETD    =  850.0    ; % ms                      % Readout time T_Ro = ETL * ESP (GE: Check CVs after Download)
    iC.S.IEf    =  100      ; % %                       % Inversion efficiency in percent (are GE's SECH pulses ~98%?)
    iC.S.IPD    =   16.0    ; % ms                      % Inversion pulse duration (GE/Phi used 16/17 ms HypSec in tests)
    iC.S.T2pOn  = true      ;                           % WIP: Turn this off when T1n is set manually (as T2 is then unknown)?
    iC.S.T2p    =    0.0    ; % ms                      % T2 preparation time
    iC.S.TOb    =   25.0    ; % ms  % WIP: Relate to iC.S.IPD   % TI1 outboard (Siemens: Add around 23 ms; +10? ms for T2p?)
    iC.S.FAx    =   90      ; % °                       % (MEX = cos(FAx*pi/180) is MagZ after excitation with angle FAx)
    iC.S.ESP    = 3.65; iC.S.IET = 0.50; % iC.S.ETL = 122;  % Echo Spacing, Echo Train Length and Inv.Eff. for TSE readout
                                                        % WIP: Just assuming incomplete inversion doesn't lead to anything!
                                                        % A proper Bloch simulation, at least w/ [ Mz, Mt ] is required.
                                                        % Actual flip angles are swept from about 10-170° to preserve Mz!
    iC.P.RCS    =   0.0     ; % fraction                % Residual CSF signal at readout. Partial nulling may increase CNR.
    iC.S.Rew    =   0.0     ; %                         % TODO: Rewinder after readout, regaining some of the pre-readout Mz?
    UI1.Show    =  true     ;                           % Show the UI control panel
    iC.P.Mode   =     1     ;                           % Start in mode 1, then keep the run mode over reruns
    iC.P.TsSel  =     3     ;                           % Start with a value (CSF), then keep the selection over reruns
    iC.P.B0sel  =     2     ;                           % Start with GE 3T setup, w/ relaxation times from literature/source
    SetRelaxTimes( iC.P.B0sel )                         % Need to initialize relaxation times before mode
    SetMode( iC.P.Mode )                                % Initialize mode settings
    iC_Ini1 = true;                                     % The program has been initialized (this var will stay upon rerunning)
end % if initialized

%% MAIN
main();
function main()
    iC.S.TRef = iC.S.TR - iC.S.ETD - iC.S.IPD       ;   % Recalculate based on current parameter settings
    if ( iC.S.T2pOn )
        iC.S.TRef = iC.S.TRef - iC.S.T2p            ;
    end % if
    iC.S.IE   = iC.S.IEf/100                        ;   % MagZ to recover after incomplete inv. (was log(2/...) etc)
    iC.S.ETL  = ceil(iC.S.ETD/iC.S.ESP)             ;   % Echo train length (# pulses) for TSE readout (ceil or floor best?)
    
    clf                                                 % Clear the figure so we can change it
    switch iC.P.Mode
        case 1                                          % Mode 1 (1IR tissue nulling)
            mainSIRNul( iC.T.Tn1 )                  ;
        case 2                                          % Mode 2 (DIR tissue nulling)
            iC.S.TRef = iC.S.TRef - iC.S.IPD        ;   % Subtract one more inversion pulse duration for DIR        WIP
            mainDIRNul( iC.T.Tn1, iC.T.Tn2 )        ;
        case 3                                          % Mode 3 (T1W nulling)
            iC.S.TRef = iC.S.TRef - iC.S.IPD        ;   % Subtract one more inversion pulse duration for DIR        WIP
            Tn2 = mainT1WNul()                      ;   % ( T1_CSF, T1_WM, T1_MS )  % TODO: Set iC.T.Tn2 = mainT1... here?
    end
    IR_CreateUI();
%     if saveImg, saveFigAsImg(),end                    % Automatically save the fig. as an image (or use the fig. menu)
end % main

% 1) Calculate TI in a 1IR/STIR/FLAIR sequence, nulling a specified T1 time
function mainSIRNul( T1z1 )
    iC.T.T1z = T1z1                                 ;   % What T1 are we concerned with nulling?
    Tis = iC.M.T(:,iC.P.TsSel)                      ;   % What tissue is selected? (Needed for T2prep)
    E2  = CalcPE2( Tis )                            ;   % Calculate E2 if using T2prep (=1 otherwise)
    if ( iC.S.Vnd == "GE" )
        iC.M.Mze  = [ 0.0, 0.0, 0.045, 0.0, 0.0 ]   ;   % GE "Mz @ end of CUBE train w/ 130° max FA"; only used for FLAIR
    else
        iC.M.Mze  = [ 0.0, 0.0, 0.0  , 0.0, 0.0 ]   ;   % [ WM, GM, CSF, MS, Fat ] Mz @ readout end (only used by GE or...?)
    end % if
    Mze = iC.M.Mze( iC.P.TsSel )                    ;   % End Mz for the tissue (GE value for XETA/Cube refocusing train)
    iC.T.Ti1n = CalcTInT1( iC.T.T1z, E2, 1-Mze )    ;   % Calc. TI to null T1z
    iC.T.Ti1n = uint16(round( iC.T.Ti1n + iC.S.IPD ))   ;   % Add the inversion pulse duration; make value int for plot display
    
    TInts = [ 0, iC.T.Ti1n, iC.S.TR ]               ;   % SIR: Inversion is at t = 0, then readout at TI
    IR_T1_MagPlot( TInts )                          ;   % Plot Z magnetization (Mz) of several tissue T1s over time
    
    figTitle( "Single IR tissue nulling" )              % Add an informative title to the plot
    fStr = ['\\itTissue T1 time:\n\\rm'             ...
            '  T1_{n}   = %4.i ms\n'                ...
            '  \n'                                  ...
            '\\itTissue nulling IR:\n\\rm\\bf'      ...
            '  TI_{n}    = %4.i ms'             ]   ;
    fVar = [ uint16(iC.T.T1z), uint16(iC.T.Ti1n) ]  ;
    if ( iC.S.Vnd == "Siemens" )
        fStr = [ fStr,                              ...
            '\n\\rm  TI_{Sie} = %4.i ms'        ]   ;
        fVar = [ fVar, SieT(iC.T.Ti1n)  ]           ;
    end % if
    figTextBox( fStr, fVar )                        ;   % Display an info text box on the figure
    
    iC.P.pStr = [                           '\n'    ...
            '*******************************\n'     ...
            '***    IR tissue nulling:   ***\n'     ...
            '*******************************\n'     ...
            '*  Tissue T1       = %4.i ms  *\n'     ...
            '*******************************\n'     ...
            '*  Nulling TI      = %4.i ms  *\n'     ...
            '*******************************\n' ]   ;
    iC.P.pVar = [ iC.T.T1z, iC.T.Ti1n ]             ;
    if ( iC.S.Vnd == "Siemens" )
        iC.P.pStr = strcat( iC.P.pStr, [            ...
            '***      For Siemens:       ***\n'     ...
            '*******************************\n'     ...
            '*  TI              = %4.i ms  *\n'     ...
            '*******************************\n' ] ) ;
        iC.P.pVar = horzcat( iC.P.pVar,SieT(iC.T.Ti1n) ) ;  % Take into account T2p and some outboard time for Siemens
    end % if
%     if ( iC.S.Vnd == "GE" )
%         iC.P.pStr = strcat( iC.P.pStr, [            ...
%             '***      For GE (DV25):     ***\n'     ...
%             '*******************************\n'     ...
%             '*  TI              = %4.i ms  *\n'     ...
%             '*******************************\n' ] ) ;
%         iC.P.pVar = horzcat( iC.P.pVar, iC.T.Ti1n + 82 );   % Estimated from experiments by WN/ØBG 2018-09-26
%     end % if                                                % NOTE: We do better now, with the new GE corrections!

end % mainSIRNul

% 2) Calculate TI1+TI2 in a DIR sequence, nulling two specified T1 times (T1_n1 = T1_CSF by default)
function mainDIRNul( T1z1, T1z2 )
    iC.T.T1z = T1z2                                 ;   % What T1 are we concerned with nulling?
    Tis = [ iC.M.T(:,3), iC.M.T(:,iC.P.TsSel) ]     ;   % CSF and selected tissue (needed for T2prep)
    E2 = CalcPE2( Tis )                             ;   % Calculate E2 if using T2prep (=1 otherwise)
    [ iC.T.Ti1n, iC.T.Ti2n ] = Find2TIn2T1([T1z1,T1z2],E2); % Inversion times that null 2 T1
    iC.T.Ti1n = uint16(round( iC.T.Ti1n + iC.S.IPD*2 )) ;   % Add the inversion pulse durations
    iC.T.Ti2n = uint16(round( iC.T.Ti2n + iC.S.IPD   )) ;   % Cast as int for plot display; ceil()/round() don't change syms!
    
    TInts = [ 0, (iC.T.Ti1n-iC.T.Ti2n), iC.T.Ti1n, iC.S.TR ] ;  % DIR: Vector of all time points (inversion, readout and end)
    IR_T1_MagPlot( TInts )                          ;   % Plot Z magnetization over time [PlotMagZ used to return end MagZ vector]

    figTitle( "DIR dual tissue nulling" )               % Add an informative title to the plot
    
    fStr = ['\\itTissue T1 times:\n\\rm'            ...
            '  T1_{n1}  = %4.i ms\n'                ...
            '  T1_{n2}  = %4.i ms\n'                ...
            '  \n'                                  ...
            '\\itTissue nulling DIR:\n\\rm\\bf'     ...
            '  TI_{1}    = %4.i ms\n'               ...
            '  TI_{2}    = %4.i ms'             ]   ;
    fVar = [ T1z1, T1z2, iC.T.Ti1n, iC.T.Ti2n ]     ;
    figTextBox( fStr, fVar )                        ;   % Display an info text box on the figure
    
    iC.P.pStr = [                           '\n'    ...
            '*******************************\n'     ...
            '***   DIR tissue nulling:   ***\n'     ...
            '*******************************\n'     ...
            '*  T1 tissue 1     = %4.i ms  *\n'     ...
            '*  T1 tissue 2     = %4.i ms  *\n'     ...
            '*******************************\n'     ...
            '*  TI1             = %4.i ms  *\n'     ...
            '*  TI2             = %4.i ms  *\n'     ...
            '*******************************\n' ]   ;
%             '*  TI1 - TI2       = %4.i ms  *\n'     ...
    iC.P.pVar = [ T1z1, T1z2, iC.T.Ti1n, iC.T.Ti2n ]    ;   % , (iC.T.Ti1n-iC.T.Ti2n)
    if ( iC.S.Vnd == "GE" )
        iC.P.pStr = strcat( iC.P.pStr, [            ...
            '***   Times for GE (DV25):  ***\n'     ...
            '*******************************\n'     ...
            '*  TI1             = %4.i ms  *\n'     ... % Estimated from experiments by WN/ØBG 2018-09-26
            '*  TI2             = %4.i ms  *\n'     ... % --"--
            '*******************************\n' ] ) ;
        iC.P.pVar = horzcat( iC.P.pVar, [ uint16( iC.T.Ti1n + iC.S.TRef*0.038 - 80 ), iC.T.Ti2n + 12 ] ) ;
    end % if

end % mainDIRNul

% 3) Calculate & plot MRI S0 of WM vs MS/WML in a DIR sequence, and TI1+TI2 nulling their T1 weighting
function T1z2 = mainT1WNul()                            % T1_nulled = T1WNul( T1_CSF, T1_WM, T1_MS )
    T1W = iC.M.T(1,1) ; T1L = iC.M.T(1,4) ; % T1C = iC.M.T(1,3) ;   % (CSF, WM, WML)
%     T2C = iC.M.T(2,3) ; T2W = iC.M.T(2,1) ; T2L = iC.M.T(2,4) ;   % (--"--)
    iC.T.T1z = iC.T.Tn1 ;   % iC.M.T(1,iC.P.TsSel)  ;               % What T1 are we concerned with nulling? TODO: T1C? Allow T1_n1.
    T1C = iC.T.T1z                                  ;   % Override T1C with the set value (we keep T2C from selection)
    Tis = iC.M.T(:, [ 3 1 4 ] )                     ;   % CSF, WM, WML (needed for T2prep)
    E2 = CalcPE2( Tis )                             ;   % Calculate E2 if using T2prep (=1 otherwise)
    E2C = E2(1)    ; E2W = E2(2)    ; E2L = E2(3)   ;   % (CSF, WM, WML)
    iC.T.Ti2max = T1C*log( (1+iC.S.IE)/(1+iC.S.IE*exp(-iC.S.TRef/T1C)) ); % Max TI2 nulling T1C, usually CSF (@ max TI1 = TReff)
%     iC.T.Ti2max = T1C*log(2/(1 + 2*exp(-iC.S.TR/T1C) - exp(-iC.S.TRef/T1C))); % Old Ti2max calculation without InvEff/T2P
    Ti2 = 1:1:iC.T.Ti2max; DUR = length(Ti2)        ;   % Vector of TI2 values 1 ms apart (could've used linspace() here)
    Ti1 = zeros(DUR,1)                              ;   % Vector of TI1 values for all TI2s; replace by TI1s(TI2s)?
    S0_DIR = zeros(DUR,2)                           ;   % Array of two DIR signals@TE=0 for TI2s
    
%     TI1s(t) = CalcTI1(t,T1C);                         % TI1(TI2) calc.
    for i = 1:DUR
        Ti1(i)      = CalcTi1D(Ti2(i),T1C,E2C)      ;   % TI1(TI2) calc. for T1_CSF
        S0_DIR(i,1) = CalcS0([Ti1(i),Ti2(i)],T1W,E2W) ; % S0(TI2) calc. (was CalcS0D(TI1(i),TI2(i),T1W) )
        S0_DIR(i,2) = CalcS0([Ti1(i),Ti2(i)],T1L,E2L) ; % --"--
    end
    
    iC.T.Ti2n = uint16( FindTi2nT1W(Tis(1,:),E2) )  ;   % Find TI2 to null the DIR T1 weighting between two T1s (T1C, T1W, T1L)
    iC.T.Ti1n = uint16( CalcTi1D(iC.T.Ti2n,T1C,E2C) ) ; % TI1(TI2) calc. by function (uint to avoid e+ notation in fig. text)
%     iC.T.Ti1n = uint16(TI1s( iC.T.Ti2n ));            % TI1(TI2) calc. by symbolic expr. (slower?!)
    
    % NOTE: GE operates with a guessed T2 when T1 is set manually. DV25: T2wm; RX27: T2gm. But is it used for T2prep?
    % WIP: For now, let's ignore T2prep w/ E=1. Judging from observation, that's what GE does too!?
    iC.T.Tn2 = uint16( FindT1null([iC.T.Ti1n,iC.T.Ti2n],1) ); % Second tissue nulled by DIR. UInt for display.
    T1z2 = iC.T.Tn2                                 ;   % This function outputs the fictive T1 nulled
    iC.S.newS = S0_DIR(iC.T.Ti2n,1)                 ;   % WM signal at readout; = WML signal here
    
    hold on
    plot( Ti2, S0_DIR(:,1), 'r-', 'LineWidth', 1.0 )    % Plot the signals calculated above
    plot( Ti2, S0_DIR(:,2), 'b-', 'LineWidth', 1.0 )
%     S0W_t   = CalcS0([iC.T.Ti1n,iC.T.Ti2n],T1W,E2W)     ;   % DEBUG: S0(Ti1,Ti2) calc.
%     S0L_t   = CalcS0([iC.T.Ti1n,iC.T.Ti2n],T1L,E2L)     ;   % --"--
%     disp( "    S0W,      S0L      (Debug T1null)" )
%     disp( [ S0W_t, S0L_t ] )
%     plot(Ti2,Ti1)                                   ;   % DEBUG: Plot the TI1 calculated for each TI2
%     S0_DIFF = abs(S0_DIR(:,1)-S0_DIR(:,2))          ;   % DEBUG: Calculate signal difference between two T1s
%     [~,minind] = min(S0_DIFF)                       ;   % Find the T1Wnull crossing graphically (as [minval,minind])
%     iC.T.Ti2n = Ti2(minind)                         ;
%     plot(Ti2,S0_DIFF)                               ;   % DEBUG: Plot the WM/MS signal difference
    hold off
    
%     whitebg( [ 0.70 1.00 0.70 ] )                   ;   % DEBUG: See white plot/GUI elements
%     whitebg( 'white' )                              ;   % NOTE: Setting this also resets the figure color.
    figTitle( "T1-DIR (T1W nulling)" )                  % Add an informative title to the plot
    ylabel('Rel. MR signal (TE = 0 ms)' , 'FontSize' , 12 )
    xlabel('TI_2 inversion time (ms)'   , 'FontSize' , 12 )
    xlm = [ 0 200*ceil( iC.T.Ti2max/200 ) ]; xlim(xlm); % Round up to 200 on the x axis
    ylm = [ -1 1 ];                          ylim(ylm); % Rel. S0 is in [-1,1]
    
    LegT = [ 1 4 ]                                  ;   % Make a legend for tissue 1(WM) & 4(MS)...
    NT1s = length(LegT); LegS = strings(1,NT1s)     ;
    for i = 1:NT1s, j = LegT(i)                     ;   % ...showing the tissue short names and T1 times
        LegS(i) = sprintf('%s (T1 = %i ms)', iC.M.Tag(1,j), iC.M.T(1,j) );
    end
    legend( LegS(1) , LegS(2) ,                     ... % Was: legend( 'WM', 'MS', ...
            'Location'  , 'NorthWest'           )   ;   % 'Best' may conflict with UI and TextBox
    
    lxlm = [ iC.T.Ti2n iC.T.Ti2n ]                  ;   % S0 crossing x (=TI2) value
%     lylm = [ ylm(1) S0_DIR(lxlm(1),1) ];              % Lower/upper y of line (=S0 curve)
    lylm = [ ylm(1) 0 ]                             ;   % Lower/upper y of line (= 0)
    line( lxlm, lylm,                               ... % Vertical line @ crossing
        'Color'         , 'black'               ,   ...
        'LineWidth'     , 1.25                  ,   ...
        'DisplayName'   , 'TI_{2} nullT1W'      ,   ... % (Tip: If declared after legend, the line appears in it)
        'LineStyle'     , ':'                   )   ;
    
%     fStr = ['\\itTissue T1 times:\n\\rm'            ...
%             '  T1_{WM} = %4.i ms\n'                 ...
%             '  T1_{MS} = %4.i ms\n\n'               ...
    fStr = ['\\itT1W-0 DIR times:\n\\rm\\bf'        ... % Tex \bf\it\rm = bold/italic/normal; escape \ for sprintf!
            '  TI_{1}    = %4.i ms\n'               ...
            '  TI_{2}    = %4.i ms\n'               ...
            '  T1_{nt}  = %4.i ms'              ]   ;
    fVar = [ iC.T.Ti1n, iC.T.Ti2n, T1z2 ]               ;   % Was: T1_WM, T1_MS, etc
    figTextBox( fStr, fVar )                        ;   % Display an info text box on the figure

    iC.P.pStr = [                           '\n'    ...
            '*******************************\n'     ...
            '***    DIR T1-W nulling:    ***\n'     ...
            '*******************************\n'     ...
            '*  T1 nulled (CSF) = %4.i ms  *\n'     ...
            '*  T1 also nulled  = %4.i ms  *\n'     ...
            '*******************************\n'     ...
            '*  TI1             = %4.i ms  *\n'     ...
            '*  TI2             = %4.i ms  *\n'     ...
            '*  TI1 - TI2       = %4.i ms  *\n'     ...
            '*******************************\n' ]   ;
    iC.P.pVar = [ T1C, T1z2, iC.T.Ti1n, iC.T.Ti2n, (iC.T.Ti1n-iC.T.Ti2n) ];
    
end % mainT1WNul

%% FUNCTIONS
function FltSym = FoS( inNr )                           % Cast as double for, e.g., exp() unless it's used as a sym
    if  isa( inNr , 'sym' ) , FltSym = inNr         ;   % TODO: This exists both in the main script and T2_SigPlot now.
    else,                     FltSym = double(inNr) ;
    end % if
end % fcn

function SieTI = SieT( TI )                             % Add outboard time estimated for the Siemens sequence
    SieTI = TI + iC.S.TOb                           ;   % TODO: Rework this for Sie+Phi? Can use iC.S.TOb for Sie but not Phi?
    if ( iC.S.T2pOn ) && ( iC.S.T2p > 0 )               % Add a little extra outboard (10 ms?) for T2 prep?
        SieTI = SieTI + iC.S.T2p                    ;   % Was + 10 for T2prep, but that's iffy...
    end % if
    SieTI = uint16(SieTI)                           ;
end % fcn

function E2 = CalcPE2( T1T2 )                           % Calculate E2 for T2prep for a (vector of) T2 time(s)
    if ( iC.S.T2pOn ) && ( size(T1T2,1) > 1 )           % Only use the T2prep formula if T2 is known
        E2 = exp( -iC.S.T2p./T1T2(2,:) )            ;
    else
        E2 = ones(1,size(T1T2,2))                   ;
    end % if
    E2 = double(E2)                                 ;
end % fcn

function S0_IR      = CalcS0( Ti, T1, E2 )              % Calculate generic IR relative S0 (def: S/M0 = E2*S0 @ TE=0 ms)
    Ti = FoS(Ti); T1 = FoS(T1)                      ;   % Cast input as float unless it's a sym
    NI = length( Ti ); S0_IR = ones( 1,length(T1) ) ;   % Ti() are defined as time to readout
    IE  = iC.S.IE                                   ;   % Inv.Eff.
    for i = 0:NI-2    % (was 0:NInv-1)
        S0_IR = S0_IR - (-IE)^i*(1+IE).*exp(-Ti(NI-i)./T1)  ;   % Step backwards through the inversions
    end
    S0_IR = S0_IR - (-IE)^(NI-1)*(1+IE*E2).*exp(-Ti(1)./T1) ;   % The last step includes E2 for T2prep
    S0_IR = S0_IR - (-IE)^(NI-0)*E2.*exp(-iC.S.TRef./T1)    ;   % Gives one signal value per given T1
end % fcn                                               % NOTE: Inv. pulse dur. isn't incorporated in these calculations

% function S0_IR      = CalcS0S( Ti, T1, E2 )             % Calculate single IR relative S0 (def: S0 = S/(E2*M0) @ TE=0)
%     Ti = FoS(Ti); T1 = FoS(T1)                      ;   % Cast input as float unless it's a sym
%     S0_IR = 1 - (1+iC.S.IE*E2).*exp(-Ti./T1) + iC.S.IE*E2.*exp(-iC.S.TRef./T1) ;
% end % fcn

% function S0_DIR     = CalcS0D( Ti1, Ti2, T1, E2 )       % Calculate DIR relative S0 (def: S0 = S/(E2*M0) @ TE=0)
%     Ti1 = FoS(Ti1); Ti2 = FoS(Ti2); T1 = FoS(T1)    ;   % Cast input as float unless it's a sym
%     IE  = iC.S.IE                                   ;   % Inv.Eff.
%     S0_DIR = 1 - (1+IE).*exp(-Ti2./T1) + IE*(1+IE*E2).*exp(-Ti1./T1) - IE^2*E2.*exp(-iC.S.TRef./T1) ;
% end % fcn

function TI_SIR     = CalcTInT1( T1, E2, MC )           % Calc. TI nulling T1 in a single IR sequence w/ T2prep and Mz corr.
    T1 = FoS(T1)                                    ;   % Cast input as float unless it's a sym
%     TI_SIR = log(2)*T1                            ;   % The simple case, TR>>T1, gives TIn = ln(2)*T1 ~ 0.69*T1
    TI_SIR =    ...
        T1*log( (1+iC.S.IE.*E2) ./ ( 1-iC.P.RCS + iC.S.IE*MC.*E2.*exp(-iC.S.TRef./T1) ) );  % 1IR TI nulling T1
end % fcn

function TI1_DIR    = CalcTi1D( Ti2, T1, E2 )           % Calc. TI1 (new TI1, the old one was - TI2) nulling T1 given TI2
    Ti2 = FoS(Ti2); T1 = FoS(T1)                    ;   % Cast input as float unless it's a sym
    IE  = iC.S.IE                                   ;   % Inv.Eff.
    TI1_DIR =	...
        T1.*log( IE*(1+IE*E2)./( -(1+iC.P.RCS) + (1+IE).*exp(-Ti2./T1) + IE^2*E2.*exp(-iC.S.TRef./T1) ) );  % (May be a sym)
end % fcn

function S0_nD      = CalcS0for0( Ti2, T1d, E2d )       % Calc. S0 from a DIR T1 found, combining CalcS0D+TI1D for sanity
    Ti2 = FoS(Ti2); T1d = FoS(T1d)                  ;   % Cast input as float unless it's a sym (also, using 1x1 vars here)
    T1 = T1d(1); E2 = E2d(1)                        ;
    IE  = iC.S.IE                                   ;   % Inv.Eff.
    Ti1 =       ...
        T1*log( IE*(1+IE*E2)/( -(1+iC.P.RCS) + (1+IE)*exp(-Ti2/T1) + IE^2*E2*exp(-iC.S.TRef/T1) ) ); % No .* ./ here: Only one T1
    T1 = T1d(2); E2 = E2d(2)                        ;
%     S0_nD  = CalcS0( [ Ti1, Ti2 ], T1(2), E2(2) ) ;
    S0_nD =     ...
        1 - (1+IE)*exp(-Ti2/T1) + IE*(1+IE*E2)*exp(-Ti1/T1) - IE^2*E2*exp(-iC.S.TRef/T1);
end % fcn

function [ Ti1n, Ti2n ] = Find2TIn2T1( T1d, E2 )        % Find TI1/TI2 that null two T1 (e.g., CSF and WM) in DIR
    T1 = FoS(T1d)                                   ;   % Cast input as float unless it's a sym
%     Ti2n = vpasolve( CalcS0( [ CalcTi1D(t,T1(1),E2(1)) , t ],T1(2),E2(2) ), t );
    Ti2n = vpasolve( CalcS0for0( t, T1, E2 ), t )   ;
    Ti1n = CalcTi1D( Ti2n,T1(1),E2(1) )             ;
end % fcn

function T1WnTi2    = FindTi2nT1W( T1_CWL, E2 )         % Find TI2 to null CSF and the DIR T1 weighting between WM and WML
    T1 = FoS(T1_CWL);                                   % Cast input as float unless it's a sym
    T1WnTi2 = vpasolve(                             ... % Solve S0DIR == 0 for this T1
        ( CalcS0for0( t, [ T1(1) T1(2) ], [ E2(1) E2(2) ] )         ...
        - CalcS0for0( t, [ T1(1) T1(3) ], [ E2(1) E2(3) ] ) ), t )  ;   % Return TI2 (calc. TI1(TI2) later)
end % fcn

function T1n_DIR    = FindT1null( Ti12, E2 )            % Find the second T1 time nulled by DIR
    Ti12 = FoS(Ti12)                                ;   % Cast input as float unless it's a sym
%     T1n_DIR = vpasolve(CalcS0D(Ti1o,Ti2o,t), t)         % Solve S0DIR == 0 for T1
    T1n_DIR = vpasolve( CalcS0(Ti12,t,E2), t)       ;   % Solve S0DIR == 0 for T1
end % fcn

function SetMode( mode )                                % Sets the run mode and resets the Nulling selection (iC.P.TsSel)
    % WIP: Merge SetMode with SetT1n()?!
    oldMode = iC.P.Mode                             ;   % WIP: If oldMode = mode, just skip to SetT1n!? But not on first run?
%     iC.T.Tn1 = iC.M.T(1,3)                        ;   % Default for all modes: Null out CSF as the longest T1
    switch mode
        case 1                                          % 1IR default: Null out CSF (=FLAIR)
            iC.T.Tn1 = iC.M.T(1,3)                  ;
            SetT1n( 1, 3 )                          ;
        case 2                                          % DIR default: Null out WM as tissue 2
            if ( oldMode == 3 )
                SetT1n( 2, -1 )                     ;   % A trick to notify SetT1n() that we've done T1 nulling
            else
                SetT1n( 2, 1 )                      ;
            end % if
        case 3                                          % T1W-DIR default: (Null out CSF as the longest T1)
    end % switch
    iC.P.Mode = mode                                ;   % Make the mode setting global
    if ( mode ~= 3 ) && ( oldMode ~= 3 )
        SetRelaxTimes( iC.P.B0sel )                     % Vendor implementation and mode may affect relaxation time settings
    end % if
end % fcn

function SetT1n( mode, T1s )                            % Sets the Nulling selection and T1_n#
    if ( T1s > 0 )
        switch mode
            case 1                                      % 1IR (only one tissue nulled)
                iC.T.Tn1 = iC.M.T(1,T1s)            ;
            otherwise                                   % DIR (set second T1 to null; the first is CSF by default)
                if iC.M.T(1,T1s) >= iC.T.Tn1, T1s=1; end    % Don't null Tn1/CSF twice; use WM as tissue 2 instead then.
                iC.T.Tn2 = iC.M.T(1,T1s)            ;
        end % switch
    else
        T1s = -T1s                                  ;
    end % if
    iC.P.TsSel = T1s                                ;
    
    switch mode                                         % Reset which T1 times to plot, based on mode. TODO: Change this?!?
        case 1                                          % 1IR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4                                 ;
        case 2                                          % DIR: Plot Mz for brain tissues (WM, GM, CSF, WML?)
            i = 1:4                                 ;
        case 3                                          % T1W-DIR: Plot the signal for WM and WML
            i = [ 1 4 ]                             ;   % NOTE: This is not implemented yet; Mode 3 uses its own plot
    end % switch
    iC.M.Plt = iC.M.T(:,i); iC.M.Leg = iC.M.Tag(1,i)        ;   % Update the T1 times to plot and their label texts
end % fcn

function SetRelaxTimes( B0opt )
    T1_CSF = IR_SetRelaxTimes( B0opt )              ;   % This is in a separate file to make this one shorter and things clearer
    
    iC.T.Tn1 = T1_CSF                               ;   % Default for all modes: Null out CSF as the longest T1
    switch iC.P.Mode                                    % When a scheme is selected, adjust the nulling selection accordingly.
        % TODO: Remove this switch; just call SetT1n( iC.P.TsSel ) and let the function handle iC.P.Mode!
        case 1                                          % 1IR T1 nulling
            SetT1n( 1, iC.P.TsSel )                 ;
        case 2                                          % DIR T1 nulling
            SetT1n( 2, iC.P.TsSel )                 ;
        case 3                                          % T1W nulling
    end % switch
%     toggleT2prep( true )                          ;   % T2prep is compatible with tissue selection
end % setRelaxTimes

%% PLOT

function figTitle( startStr )                           % Add an informative title to the MagZ plot
    txt = [ char(startStr) ' at ' iC.P.B0tag ', TR_{eff} = ' num2str(iC.S.TRef) ' ms' ];
    if ( iC.S.Vnd ~= "" )                               % If a vendor is specified, show it
        txt = [ txt ' (for ' char(iC.S.Vnd) ')' ]   ;
    end % if
    title(  txt              , 'FontSize' , 16 )    ;
end % fcn

function figTextBox( fStr, fVar )
    txt = sprintf( fStr, fVar );                        % Formatted string "print"
    lin = count( fStr, '\n' ) + 1;                      % Count lines in the string
    xlm = xlim; xt =  0.800*double(xlm(2));             % Max value of axis; (x,y) is a plot point!
    ylm = ylim; yt = -0.900*double(ylm(2));             % --"--
    yt = yt + 0.035*lin;                                % Compensate for number of lines
    text( xt, yt, txt                           ,   ... % Note: annotation('textbox' etc doesn't scale with text
        'EdgeColor' , 'black'                   ,   ...
        'Margin'    , 4                         ,   ...
        'LineWidth' , 1                         ,   ...
        'LineStyle' , '-'                       )   ;
end % fcn

% function saveFigAsImg()
%     imgName = ['TReff_' int2str(iC.S.TRef) '_TI1_' int2str(iC.T.Ti1n) '_TI2_'  int2str(iC.T.Ti2n) ];
%     if exist('Figures','dir')~=7, mkdir('Figures'), end
%     cd Figures
%     % TODO: What's the difference between print to image and the imwrite command?
%     print('-djpeg', '-r300', imgName )
%     cd '..'
% end % fcn

%% USER INTERFACE

function IR_CreateUI()                                  % Display a button to show/hide the UI controls
%     uiX = 74.0*fSz; uiY = 50.0*fSz            ;       % Default start position (in px) for masterUI controls
    eW = 0.065 ; eH = 0.065; eS = eH + 0.025    ;       % Button width/height/spc for masterUI controls (was [50 40 60] px)
    myAx = gca                                  ;       % A handle to the current graphics axes (plot)
    uiX = myAx.OuterPosition(3) - 0.015 - eW    ;       % Use the axis' outer x position (def. [ 0 0 1 1 ]) for UI pos.
    uiY = myAx.Position(4)      + 0.010         ;       % Use the axis' inner y position to determine UI pos.
    
    uicontrol( 'Style'  , 'pushbutton'          ,   ... % For all modes:
        'ToolTipString' , '[h] Show/hide UI'    ,   ...
        'String'        , 'Panel'               ,   ...
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-4*eS eW eH ]    ,   ...
        'Callback'  , @masterUI_callback        )   ;   % UI to show/hide the UI control panel
    
    switch iC.P.Mode
        case 1                                          % 1IR T1 nulling
            b2txt = 'DIR'                           ;
            b4txt = ''                              ;
        case 2                                          % DIR T1 nulling
            b2txt = '(FLA)IR'                       ;
            b4txt = 'T1-null'                       ;
        case 3                                          % T1W nulling
            b2txt = ''                              ;
            b4txt = 'Back'                          ;
    end % switch iC.P.Mode
    if ismember( iC.P.Mode, [ 1 2 ] )                    % If 1IR/DIR:
        cPar = [ 2 1 0 ]                            ;   % Switch to this mode from mode 1-3
        uicontrol( 'Style'  , 'pushbutton'      ,   ... % uiMode
        'ToolTipString' , '[w] Switch mode'     ,   ...
        'String'        , b2txt                 ,   ... % 'togglebutton'
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ... % Make the UI resize automatically with the figure
        'Position'  , [ uiX uiY-0*eS eW eH ]    ,   ...
        'Callback'  , {@cMode_callback,cPar}    )   ;   % UI to set run mode
        
        uicontrol( 'Style'  , 'pushbutton'      ,   ... % If 1IR/DIR:
        'ToolTipString' , '[s] T2 decay plot'   ,   ...
        'String'        , 'Plot T2'             ,   ...
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-1*eS eW eH ]    ,   ...
        'Callback'  , @plotT2_callback          )   ;   % UI to create T2 vs TE plot
    end % if
    
    if ismember( iC.P.Mode, [ 2 3 ] )                    % If DIR/T1W_n:
        cPar = [ 0 3 2 ]                            ;   % Switch to this mode from mode 1-3
        uicontrol( 'Style'  , 'pushbutton'      ,   ...
        'ToolTipString' , '[t] Switch mode'     ,   ...
        'String'        , b4txt                 ,   ...
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-2*eS eW eH ]    ,   ...
        'Callback'  , {@cMode_callback,cPar}    )   ;   % UI to change run mode to T1W nulling plot
    end % if
    
    UI1.Vnd = uicontrol( 'Style', 'pushbutton'  ,   ... % ui Vnd
        'ToolTipString' , '[v] Select vendor'   ,   ...
        'String'        , 'Vendor'              ,   ...
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-5*eS eW eH ]    ,   ...
        'Callback'  , @vendor_callback          )   ;   % UI to select vendor implementation
%     if ( iC.S.Vnd == "" )
        UI1.Vnd.String = "Vendor"                   ;
%     else
%         UI1.Vnd.String = iC.S.Vnd                 ;
%     end % if
    
    UI1.Stp = uicontrol( 'Style', 'pushbutton'  ,   ... % ui Stp
        'ToolTipString' , '[c/d] TR ±1 s'       ,   ...
        'String'        , 'Step TR'             ,   ...
        'FontWeight'    , 'bold'                ,   ...
        'Units'         , 'normalized'          ,   ...
        'Position'  , [ uiX uiY-7*eS eW eH ]    ,   ...
        'Callback'  , {@steptr_callback,1000}   )   ;   % UI to debug by increasing TR in steps

    createUIPanel()                                 ;
    switch UI1.Show
        case true ; F1_Pan.Visible = 'on'           ;
        case false; F1_Pan.Visible = 'off'          ;
    end % switch UI1.Show
end % IR_CreateUI

function createUIPanel()
    % Make UI control buttons for saving images instead of the variable, like with printVars?
    %      Probably not necessary, since you can save the image from the File menu.
    
    boxClr = iC.P.FigClr                            ;   % UI panel background color (same as fig background)
%     uiX = 58.3*fSz; uiY = 54.5*fSz                ;   % Default start position (in px) for UI controls
    eN = 11; eW = 130; eH = 30                      ;   % Num. of UI elements; default element width/height (were 120/30)
    mX = 10; mY =  6                                ;   % x/y margins in UI panel
    eW2 = eW/2 - 3; mX2 = mX + eW/2 + 3             ;   % Width/position of half-width elements
    UI1.Hpx = eN*eH + 2*mY; uiT = UI1.Hpx - mY + 3  ;   % UI panel height in pixels; top position for elements
    UI1.Wpx = eW + 2*mX                             ;   % UI panel width in pixels
    
    F1_Pan = uipanel( iC.P.Fig1                 ,   ... % F1_Pan
        'BackgroundColor'   , boxClr            ,   ... %  [ 0.8 0.8 0.4 ]
        'BorderType'        , 'etchedin'        ,   ... % default 'etchedin'
        'BorderWidth'       , 1.0               ,   ... % default 1
        'Clipping'          , 'off'             ,   ... % default 'on'; 'off' is good for debugging
        'Units'             , 'normalized'      ,   ... % Normalized units facilitate resizing to axes
        'ResizeFcn'         , @placeUiBox       )   ;   % Resize function for the UI panel
%         'Units'             , 'pixels'          ,   ... % With normalized units (def.), got Pos. [ .726 .420 .172 .500 ]
%         'Title'             , 'Controls'        ,   ... % If specifying a title, it'll show up on the panel
%         'FontSize'          , 10                ,   ...
%         'Position'  , [ uiX uiY-uiH uiW uiH ]   ,   ... % UI panel container for the controls (pixel units)
%     uistack( F1_Pan, 'top' )                        ;   % Bring the panel to the top of the graphic stack (not necessary)
    placeUiBox()                                    ;   % Adjust the UI panel position to the figure window
    
    for i = 1:eN; eY = uiT - i*eH                   ;   % Generate each UI element on a new row in the panel
    switch i
    case  1
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui pT1s
            'Style'         , 'pushbutton'      ,   ...
            'ToolTipString' , '[f] Print info'  ,   ... % 'Print T1 times'
            'String'        , 'Info'            ,   ...
            'Position'  , [ mX eY+0 eW2 24 ]    ,   ...
            'Callback'  , @printInfo_callback   )   ;   % UI to print settings to the command window
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'pushbutton'      ,   ...
            'String'        , 'Results'         ,   ...
            'ToolTipString' ,'[p] Print results',   ...
            'Position'  , [ mX2 eY+0 eW2 24 ]   ,   ...
            'Callback'  , @printVars_callback   )   ;   % UI to print results to the command window

    case 2
        uiStr = 'T1s:'                              ;
        ttStr = 'Ref. values for T1s'               ;
        UI1.B0T = uicontrol( 'Parent', F1_Pan   ,   ... % ui B0[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'      , [ mX eY+2 eW 16 ] )   ; 
        UI1.B0S = uicontrol( 'Parent', F1_Pan   ,   ...
            'Style'         , 'popup'           ,   ...
            'String'        , iC.P.B0str        ,   ...
            'Value'         , iC.P.B0sel        ,   ...
            'Position'  ,[ mX+24 eY+0 eW-34 23 ],   ... % Note: Height setting won't change the actual height here
            'Callback'  , @b0_sel_callback      )   ;   % UI to set relaxation times based on B0 and ref.
        
    case  3
        uiStr = 'TR:                          ms'   ;
        ttStr = '[u] Repetition time'               ;
        t2Str = 'TR = TR_ef + T_Ro + T_T2p'         ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui TR[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        UI1.TRS = uicontrol( 'Parent', F1_Pan   ,   ...
            'Style'         , 'edit'            ,   ...
            'String'        , iC.S.TR           ,   ...
            'Value'         , iC.S.TR           ,   ...
            'Min' , 11      , 'Max' , 11        ,   ... % Tip: Max > Min allows multiline edit
            'ToolTipString' , t2Str             ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @tr_set_callback      )   ;   % UI to set TR

    case  4
        uiStr = 'T_Ro:                      ms'     ;
        ttStr = 'TSE Readout time'                  ;
        t2Str = 'ETD = ETL * ES'                    ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui TRo[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , t2Str             ,   ...
            'String'        , iC.S.ETD          ,   ...
            'Value'         , iC.S.ETD          ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @tro_set_callback     )   ;   % UI to set T_Read

    case  5
        uiStr = 'Inv.eff.:                   %'     ;
        ttStr = 'Inversion efficiency'              ;
        t2Str = 'IEf = -cos(FA_inv)'                ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui IEf[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+2 eW 16 ]     )   ;
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , t2Str             ,   ...
            'String'        , iC.S.IEf          ,   ...
            'Value'         , iC.S.IEf          ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @ie_set_callback      )   ;   % UI to set inversion efficiency

    case  6
        uiStr = 'T_Inv.:                     ms'    ;
        ttStr = 'Inv. pulse dur.'                   ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui IPD[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , iC.S.IPD          ,   ...
            'Value'         , iC.S.IPD          ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @ipd_set_callback     )   ;   % UI to set inversion pulse duration

    case  7
        uiStr = 'T2p:                        ms'    ;
        ttStr = 'T2 preparation'                    ;
        t2Str = 'E2p = exp(-T2p/T2)'                ;
        t3Str = '[2] T2 prep. On/Off'               ;
        UI1.T2pT = uicontrol( 'Parent', F1_Pan  ,   ... % ui T2p[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
%             'Enable'        ,        'Inactive' ,   ... % TODO: Make the text clickable, reflecting state?
%             'ButtonDownFcn','toggleT2prep( -1 )')   ;   % Problem: Fcn will be run in main data space so no local fns.
        UI1.T2pS = uicontrol( 'Parent', F1_Pan  ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , t2Str             ,   ...
            'String'        , iC.S.T2p          ,   ...
            'Value'         , iC.S.T2p          ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @t2p_set_callback     )   ;   % UI to set T2prep duration (in ms)
%         UI1.T2pS.ForegroundColor  = iC.P.Gray       ;   % Gray out the T2prep text for now since it's WIP
        UI1.T2pB = uicontrol( 'Parent', F1_Pan  ,   ... % ui T2p[T|S]
            'Style'         , 'radio'           ,   ...
            'ToolTipString' , t3Str             ,   ...
            'Position'  , [ mX+23 eY+3 14 14 ]  ,   ...
            'Callback'  , @t2p_toggle_callback  )   ;   % UI to set T2prep duration (in ms)
        if ( iC.S.T2pOn )
            UI1.T2pS.Enable = 'on'                  ;
            UI1.T2pB.Value  = 1                     ;
        else
            UI1.T2pS.Enable = 'off'                 ;
            UI1.T2pB.Value  = 0                     ;
        end % if

    case  8
        uiStr = 'T1_n1:                    ms'      ;   % Note: UIControl text can't support Tex/LaTeX
        ttStr = 'Longer T1 to null'                 ;
        UI1.Tn1T = uicontrol( 'Parent', F1_Pan  ,   ... % ui Tn1[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'edit'            ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'String'        , iC.T.Tn1          ,   ...
            'Value'         , iC.T.Tn1          ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @t1n1_set_callback    )   ;   % UI to set Tn1 (longer T1, typically CSF)

    case  9
        if ismember( iC.P.Mode , [ 2 3 ] )               % DIR T1 nulling and T1-DIR null two T1
            uiStr = 'T1_n2:                    ms'  ;
            ttStr = 'Shorter T1 to null'            ;
        UI1.Tn2T = uicontrol( 'Parent', F1_Pan  ,   ... % ui Tn2[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        UI1.Tn2S = uicontrol( 'Parent', F1_Pan  ,   ...
            'Style'         , 'edit'            ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'String'        , iC.T.Tn2          ,   ...
            'Value'         , iC.T.Tn2          ,   ...
            'Position'  , [ mX+40 eY+0 50 20 ]  ,   ...
            'Callback'  , @t1n2_set_callback    )   ;   % UI to set Tn2 (shorter T1; GM/WM/etc)
        end % if

    case 10
        uiStr = 'Null T1:'                          ;
        ttStr = 'Ref. tissue to null'               ;
        UI1.STnT = uicontrol( 'Parent', F1_Pan  ,   ... % ui STn[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment',  'Left'      ,   ...
            'Position'  , [ mX eY+1 eW 16 ]     )   ;
        if ismember( iC.P.Mode, [ 1 2 ] )                % 1IR and DIR T1 nulling use tissue selector
        UI1.STnS = uicontrol( 'Parent', F1_Pan  ,   ...
            'Style'         , 'popup'           ,   ...
            'String'        , iC.M.Tag(2,:)     ,   ...
            'Value'         , iC.P.TsSel        ,   ...
            'Position'  ,[ mX+40 eY+0 eW-50 22 ],   ...
            'Callback'  , @t1_sel_callback      )   ;   % UI to set T1_n(2) by tissue
        if iC.T.T1z ~= iC.M.T(1,iC.P.TsSel)
            UI1.STnS.ForegroundColor  = iC.P.Gray   ;   % Gray out the tissue selector if entering values manually
        end % if iC.T.T1z
        else % ismember iC.P.Mode
            UI1.STnT.String = 'Rel. S0:                     %' ;
            if ~isfield(iC.S,'oldS'); iC.S.oldS = iC.S.newS; end % if       % For use in the relSNR uicontrol box
            if ~isfield(iC.S,'oldB'); iC.S.oldB = iC.S.B0  ; end % if       % (w/ normal global var., was if isempty)
            relS0 = uint16(100*(iC.S.newS/iC.S.oldS)*(iC.S.B0/iC.S.oldB));  % Rel.S0 is prop. to CNR, and rel.SNR here
            if ( debugInfo == true )
                fprintf("   "); fprintf("%-10s", ["new S0L","old S0L","ratio","(rel. units)"] ); fprintf("\n") ;
                disp( [ iC.S.newS, iC.S.oldS, iC.S.newS/iC.S.oldS ] )
            end % if debug
            iC.S.oldS = iC.S.newS                   ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % Text box to show rel. S0 (instead of tissue selector)
            'Style'         , 'edit'            ,   ...
            'String'        , relS0             ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'Enable'        , 'Inactive'        ,   ...
            'Position'  , [ mX+40 eY+0 60 20 ]  )   ;
        end % if

    case  11
        uiStr = 'Res.S:                 %'          ;
        ttStr = 'Residual signal'                   ;
        t2Str = 'May improve CNR?'                  ;
        uicontrol( 'Parent' , F1_Pan            ,   ... % ui RS[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX eY+2 eW 16 ]     )   ;
        uicontrol( 'Parent' , F1_Pan            ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , t2Str             ,   ...
            'String'        , 100*iC.P.RCS      ,   ...
            'Value'         , 100*iC.P.RCS      ,   ...
            'Position'  , [ mX+40 eY+0 40 20 ]  ,   ...
            'Callback'  , @rs_set_callback      )   ;   % UI to set desired residual CSF signal

%     case  11
%         uicontrol( 'Parent' , F1_Pan            ,   ... % ui FAex[T|S]
%             'Style'         , 'text'            ,   ...
%             'ToolTipString','Excitation pulse FA',  ...
%             'BackgroundColor', boxClr           ,   ...
%             'HorizontalAlignment', 'Left'       ,   ...
%             'Visible'       , 'off'             ,   ... % Disable this for now (not needed)
%             'Position'  , [ mX eY+1 eW 16 ]     ,   ...
%             'String','FAex:                         °');
%         UI1.FAexS = uicontrol( 'Parent', F1_Pan ,   ...
%             'Style'         , 'edit'            ,   ...
%             'String'        , iC.S.FAx          ,   ...
%             'Value'         , iC.S.FAx          ,   ...
%             'Min' , 11      , 'Max' , 11        ,   ... % Tip: Max > Min allows multiline edit
%             'Visible'       , 'off'             ,   ... % Disable this for now (not needed)
%             'Position'  , [ mX+40 eY+0 60 20 ]  ,   ...
%             'Callback'  , @fa_set_callback      )   ;   % UI to set FA (in degrees)
        

% TODO: Make a button to show a list box for selecting what to plot, in a separate GUI panel?
%     case 11
%         UI1.ShTT = uicontrol( 'Parent', F1_Pan  ,  ... % ui ShT[T|S]
%             'Style'         , 'text'            ,   ...
%             'ToolTipString' , 'Tissues to show' ,   ...
%             'String'        , 'Showing'         ,   ...  
%             'BackgroundColor', boxClr           ,   ...
%             'HorizontalAlignment', 'Left'       ,   ...
%             'Position'  , [ mX eY+1 eW 16 ]     )   ;
%         UI1.ShTS = uicontrol( 'Parent', F1_Pan  ,   ...
%             'Style'         , 'listbox'         ,   ...
%             'String'        , iC.M.Tag(1,:)     ,   ...
%             'Value'         , iC.M.Leg          ,   ...
%             'Min' , 00      , 'Max' , 11        ,   ... % Tip: Max > Min allows multiple listbox selection
%             'Position', [ mX+40 eY-60 eW-40 82 ],   ...
%             'Callback'  , @t_disp_callback      )   ;   % UI to set which tissues to display/calculate

    end % switch i
    end % for i
    
    switch iC.P.Mode
        case 1                                          % 1IR T1 nulling
        UI1.Tn1T.String  = 'T1_n:                      ms'  ;   % There is only one T1 to null in 1IR
        UI1.STnT.Position = UI1.STnT.Position + [ 0 eH 0 0 ] ;  % Move the tissue selector one line up
        UI1.STnS.Position = UI1.STnS.Position + [ 0 eH 0 0 ] ;  % --"--
%       UI1.Tn2S.Visible = 'off'                    ;   % Hide the tissue 2 edit text/box (ui###T/S, if created)
        case 2                                          % DIR T1 nulling
%       UI1.B0S.ForegroundColor = iC.P.Gray         ;   % Gray out the B0 select text/box (ui###T/S, if created)
        case 3                                          % T1W nulling
        UI1.Tn2S.Enable           = 'inactive'      ;   % The tissue 2 selector is not editable
        UI1.Tn2S.ForegroundColor  = iC.P.Gray       ;   % Gray out the tissue 2 edit text (it's set automatically)
%       UI1.FAexS.Enable          = 'off'           ;   % Gray out the FA edit box (FA is irrelevant in this mode)
    end % switch iC.P.Mode

end % createUIPanel

function val = str2val( str )                           % Function to process numeric input for UI control callback fcns
    val = str2double( str )                         ;
    if ( val ~= val )                                   % test for NaN from malformed input
        val = 0                                     ;
    end % if
end % fcn

function masterUI_callback(~,~)                         % UI that shows/hides the UI control panel
    UI1.Show = ~UI1.Show                            ;
    switch UI1.Show
        case true ; F1_Pan.Visible = 'on'           ;
        case false; F1_Pan.Visible = 'off'          ;
    end
end % fcn

function cMode_callback(~,~,toMode)                     % Switch between calculation modes
    switch iC.P.Mode
        case 1                                          % 1IR T1 nulling    ->|
            SetMode( toMode(1) )                    ;
        case 2                                          % DIR T1 nulling    |<-
            SetMode( toMode(2) )                    ;
        case 3                                          % T1W nulling       <<-
            SetMode( toMode(3) )                    ;
    end % switch                                        % Note that the figure title shows the active mode
    main();
end % fcn

function vendor_callback(~,~)                           % Switch between vendor implementations
    switch iC.S.Vnd
        case ""
            iC.S.Vnd = "GE"                         ;   % -> GE
            iC.S.T2pOld = iC.S.T2p                  ;   % Store the current T2 prep time
            iC.S.T2p = 200.0                        ;   % GE uses this value consistently
        case "GE"
            iC.S.Vnd = "Siemens"                    ;   % -> Sie
            iC.S.T2p = 170.0                        ;   % Sie uses this value consistently
%             if FLAIR, use T2P, if DIR turn it off by def.? Sie doesn't yet provide T2P w/ DIR.
        case "Siemens"
            iC.S.Vnd = "Philips"                    ;   % -> Phi
            iC.S.T2p = 125.0                        ;   % Phi uses this value as default
        case "Philips"
            iC.S.Vnd = ""                           ;   % <<-
            iC.S.T2p = iC.S.T2pOld                  ;   % Retrieve the previous T2 prep time
    end % switch                                        % Note that the figure title shows the active vendor
    SetRelaxTimes( iC.P.B0sel )                         % Vendor implementation and mode may affect relaxation time settings
    main();
end % fcn

function plotT2_callback(~,~)                           % Run T2 plot (in separate script)
    iT2s        =   struct(                         ...
        'Fig1'  ,   iC.P.Fig1                   ,   ...
        'S0'    ,   iC.R.S0(1:4)                ,   ... % S0
        'T2'    ,   iC.M.T(2,:)                 ,   ... % T2
        'TsTag' ,   iC.M.Tag(1,1:4)             ,   ... % Tags
        'T2tis' ,   [ 4 1 2 ]                   )   ;   % Tissues to T2 plot: WML, WM, GM
    IR_T2_SigPlot( iT2s )                           ;
end % fcn

function b0_sel_callback(src,~)                         % UI to set rel. times based on B0/ref.
    SetRelaxTimes( src.Value )                      ;   % Was get(src,'Value') pre-R2014b
    main();
end % fcn

function tr_set_callback(src,~)                         % UI to set TR
    iC.S.TR = str2val( src.String )                 ;   % str2double( src.String );
    Lim = iC.S.ETD + 5; if iC.S.TR < Lim, iC.S.TR = Lim; end  % TR > T_Read
    main();
end % fcn

function steptr_callback(~,~,inc)                       % UI to step TR (DEBUG)
    iC.S.TR = iC.S.TR + inc                         ;
    main();
end % fcn

function tro_set_callback(src,~)                        % UI to set T_Read
    iC.S.ETD = str2val( src.String )                ;
    Lim = iC.S.TR - 5; if iC.S.ETD > Lim, iC.S.ETD = Lim; end % TR > T_Read
    main();
end % fcn

function ipd_set_callback(src,~)                        % UI to set Inv. pulse dur.
    Lim = (iC.S.TRef + iC.S.IPD) - 1000             ;
    if ( iC.P.Mode ~= 1 )
        Lim = double( iC.T.Ti2n ) - 50              ;   % TODO: The second TI needs recalculation!
    end
    iC.S.IPD = str2val( src.String )                ;
    if iC.S.IPD >= Lim, iC.S.IPD = Lim; end
    if iC.S.IPD <  0  , iC.S.IPD = 0  ; end
    main();
end % fcn

function t1n1_set_callback(src,~)                       % UI to set Tn1 (longer T1)
    iC.T.Tn1 = str2val( src.String )                ;
    if iC.P.Mode > 1                                     % Using multiple inversions...
        Lim = iC.T.Tn2 + 1                          ;
        if iC.T.Tn1 < Lim, iC.T.Tn1 = Lim; end              % Tn1 is the longest T1 time to null
    end
%   toggleT2prep( false )                           ;   % T2prep is incompatible with setting T1 manually (no known T2)
    main();
end % fcn

function t1n2_set_callback(src,~)                       % UI to set Tn2 (shorter T1)
    iC.T.Tn2 = str2val( src.String )                ;
    Lim = iC.T.Tn1 - 1                              ;
    if iC.T.Tn2 > Lim, iC.T.Tn2 = Lim; end              % Tn1 is the longest T1 time to null
%   toggleT2prep( false )                           ;   % T2prep is incompatible with setting T1 manually (no known T2)
    main();
end % fcn

function t1_sel_callback(src,~)                         % UI to select a tissue T1 to null
    SetT1n( iC.P.Mode, src.Value )                  ;
%   toggleT2prep( true )                            ;   % T2prep is compatible with tissue selection
    main();
end % fcn

% function fa_set_callback(src,~)                         % UI to select the flip angle (in degrees)
%     iC.S.FAx = str2val( src.String )                ;
%     main();
% end % fcn

function t2p_set_callback(src,~)                        % UI to select the T2 preparation time (in ms)
    iC.S.T2p = str2val( src.String )                ;
    main();
end % fcn

function ie_set_callback(src,~)                         % UI to select the inversion efficiency (in percent)
    iC.S.IEf = str2val( src.String )                ;
    if iC.S.IEf <= 0, iC.S.IEf = 1; end                 % Zero/negative Inv.Eff. crashes the program
    main();
end % fcn

function rs_set_callback(src,~)                         % UI to select the desired residual CSF signal (in percent)
    iC.P.RCS = str2val( src.String )/100            ;
    Lim = 0.66                                      ;   % Too high residual signal crashes the program
    if iC.P.RCS >  Lim, iC.P.RCS =  Lim ; end
    if iC.P.RCS < -Lim, iC.P.RCS = -Lim ; end
    main();
end % fcn

function printVars_callback(~,~)                        % UI callback that prints results to the command window
    fprintf( iC.P.pStr, iC.P.pVar )                 ;
    if ( debugInfo ) && ismember( iC.P.Mode, [ 1 2 ] )  % Debug info for plot
%         len = length(iC.R.MzIt)                     ;
%         fprintf('\nMagZ(TR) over the first %i repetitions (starting at 1):\n', len);
%         disp( iC.R.MzIt )                           ;   % Show the iterative Z magnetization over the first TRs
%         disp( iC.R.MzIt(:,2:its+1) )                ;
        fStr = [ '\nMagZ and Signal(0) by T1 at'    ...
            ' readout with flip angle %3.1f°'       ...
            ' [S0=%1.2f*MagZ(TI)]:\n'           ]   ;
        fprintf(fStr, round(iC.S.FAx,1,'significant'), sin(iC.S.FAx*pi/180)) ;  % Heading
        fprintf("\t ") ; fprintf("%-10s", iC.M.Leg) ; fprintf("\n") ;       % Tissue legends
        if ( iC.S.FAx == 90 )                             % Report resulting S0 iff angle less than 90°
            fStr =                iC.R.S0'          ;
        else
            fStr = [ iC.R.MzTI' ; iC.R.S0' ]        ;
        end % if
        disp ( fStr )                               ;   % Mz/S (after FA° pulse) at readout
    end % if
end % fcn

function printInfo_callback(~,~)                        % UI callback that prints settings to the command window
    fStr = [                               '\n'     ...
            '*******************************\n'     ... % fprintf( fStr, iC.P.B0str(iC.P.B0sel), T1set );
            '***    Relaxation times     ***\n'     ...
            '***    for %-15s',       '  ***\n'     ...
            '*******************************\n'     ... % [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ]
            '*  T1 %-8s', '     = %+4s ms  *\n'     ... % Used %4.i when printing T1s directly
            '*  T1 %-8s', '     = %+4s ms  *\n'     ...
            '*  T1 %-8s', '     = %+4s ms  *\n'     ...
            '*  T1 %-8s', '     = %+4s ms  *\n'     ...
            '*  T1 %-8s', '     = %+4s ms  *\n'     ...
            '*******************************\n'     ...
            '*  T2 %-8s', '     = %+4s ms  *\n'     ...
            '*  T2 %-8s', '     = %+4s ms  *\n'     ...
            '*  T2 %-8s', '     = %+4s ms  *\n'     ...
            '*  T2 %-8s', '     = %+4s ms  *\n'     ...
            '*  T2 %-8s', '     = %+4s ms  *\n'     ...
            '*******************************\n'     ...
            '*         Settings:           *\n'     ...
            '*  %-11s',   '     = %+4s ms  *\n'     ...
            '*  %-11s',   '     = %+4s ms  *\n'     ...
            '*  %-11s',   '     = %+4s %%   *\n'    ...
            '*  %-11s',   '     = %+4s ms  *\n'     ...
            '*******************************\n' ]   ; 
    T1info = [ iC.M.Tag(2,1:5); num2cell(iC.M.T(1,1:5)) ] ; % Use this trick to weave strings and numbers into a string array
    T2info = [ iC.M.Tag(2,1:5); num2cell(iC.M.T(2,1:5)) ] ;
    if ( iC.S.T2pOn )
        T2prep = iC.S.T2p                           ;
    else
        T2prep = "--"                               ;
    end % if
    T1pars = [  [  "TR    ", "T_Readout", "Inv.Eff.", "T_T2prep"];          ...
        num2cell([  iC.S.TR  ,  iC.S.ETD    ,  iC.S.IEf   ,  T2prep   ]) ]  ;
    fprintf( fStr, iC.P.B0str(iC.P.B0sel), T1info{:}, T2info{:}, T1pars{:} );
end % fcn

function refresh_callback(~,~)                          % Refresh T1/T2 values
    fprintf( "<>   IR-Calc T1/T2 refreshed\n" )     ;
    SetRelaxTimes( iC.P.B0sel )                     ;
    main();
end % fcn

function keyPress_callback(~,evt)
    switch evt.Key
        case 'r'                                        % Refresh T1/T2 values
            refresh_callback([],[])                 ;
        case 'h'                                        % Show/hide UI panel
            masterUI_callback([],[])                ;   % []: Calling with dummy values
        case 's'                                        % Plot T2 decay
            plotT2_callback([],[])                  ;
        case 'v'                                        % Switch vendor implementation
            vendor_callback([],[])                  ;
        case 'd'                                        % Step TR (DEBUG)
            steptr_callback([],[], 1000)            ;
        case 'c'                                        % Step TR (DEBUG)
            steptr_callback([],[],-1000)            ;
        case 't'                                        % TrueT2-DIR shortcut
            switch iC.P.Mode
                case 1
                    cMode_callback([],[],[ 2 1 0 ]) ;
                case 2                                  % T1-nulling, plot T2 decay, focus back on this fig.
                    cMode_callback([],[],[ 0 3 2 ]) ;
                    cMode_callback([],[],[ 0 3 2 ]) ;
                    plotT2_callback([],[])          ;
                    figure( iC.P.Fig1 )             ;
                case 3
                    cMode_callback([],[],[ 0 3 2 ]) ;
            end % switch Mode
        case 'f'                                        % Print info
            printInfo_callback([],[])               ;
        case 'p'                                        % Print results
            printVars_callback([],[])               ;
        case 'w'                                        % Switch IR mode
            if ismember( iC.P.Mode, [ 1 2 ] )
                cMode_callback([],[],[ 2 1 0 ])     ;
            end % if
        case 'a'                                        % Change active figure
            if isfield( iC.P, 'Fig2' )                  % ishandle( iC.P.Fig2 )
                figure( iC.P.Fig2    )              ;
            end % if
        case 'u'                                        % Focus on TR setting uicontrol
            uicontrol( UI1.TRS )                    ;
        case '2'                                        % Toggle T2prep on/off
            t2p_toggle_callback                     ;
    end % switch
end % fcn

function placeUiBox(~,~)                                % Function to resize the UI panel (automatically called)
    if ishandle( F1_Pan )                               % (The first time this fn is called, the panel doesn't exist yet)
        pxFig = getpixelposition(iC.P.Fig1)         ;   % Figure size in px, for manual normalizing
        myAx = gca                                  ;   % A handle to the current graphics axes/plot
        uiW = UI1.Wpx / pxFig(3)                    ;   % Normalized UI panel width
        uiH = UI1.Hpx / pxFig(4)                    ;   % Normalized UI panel height
        uiX = myAx.Position(3) - uiW + 0.115        ;   % The axis' inner right x pos. (but why the "fudge factor"?)
        uiY = myAx.Position(4) - uiH + 0.090        ;   % The axis' inner upper y pos. (--"--)
        F1_Pan.Position = [ uiX uiY uiW uiH ]       ;
    end % if
end % fcn

function t2p_toggle_callback(~,~)                       % UI callback that toggles T2 prep on/off
    toggleT2prep(-1)                                ;
    main();
end % fcn

function toggleT2prep( set )
    if ( set == -1 )
        iC.S.T2pOn = ~iC.S.T2pOn                    ;
    else
        iC.S.T2pOn = set                            ;
    end % if set                                        % Note: Can't change UI here, as the GUI is redrawn later.
end % fcn

end % script fn