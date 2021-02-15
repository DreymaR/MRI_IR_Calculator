function T1_CSF = IR_SetRelaxTimes( B0opt )         % Set relax. times based on the selected "B0" option, return T1_CSF
%% (D)IR MRI T1/T2 settings
%   - This function is called from the IR_TI_Calculator
%   - Set T1/T2 relaxation times for relevant tissues (WM/GM/CSF, MS/WML and Fat).
%   - You may edit the Custom set of times at your leisure. Reselect the "T1s" option after saving this file.
%
%  NOTE:
%   - Relaxation times vary widely in literature. We have used values from GE's scanner implementation, and other sources.
%       - Lalande et al, MRI 2016: https://www.sciencedirect.com/science/article/pii/S0730725X16301266
%   - 
%
%  TODO:
%   - 
%
%  DONE:
%   - Moved SetRelaxTimes() to a separate file to ease editing and overview.
%


%% INIT
% debugInfo  = 1;                                     % Show extra info in the command window?

global iPr iS Ts    % iUI iT iR                         % Program/system/tissue % UI/times/results data structs

    iPr.B0sel = B0opt;                                  % Ensures right selection in the UI
    iPr.B0str = [                                   ... % Array of B0/sources for relax. times
                "1.5 T (GE)"    , "3.0 T (GE)"      ,   ...
                "1.5 T (Custom)", "3.0 T (Custom)"  ,   ...
                "1.5 T (Alsop)" , "3.0 T (Alsop)"   ,   ...
                "1.5 T (MRIQue)", "3.0 T (Lalande)" ,   ...
                "1.5 T (Visser)", "7.0 T (Visser)"  ,   ...
                "3.0 T (quant.)"                    ];  % Used in the UI
%     iPr.B0tag = char(iPr.B0str(iPr.B0sel));
%     iPr.B0tag = iPr.B0tag(1:5);                         % Char array of only the field strength, for figure titles
    TisStr = [  " WM     ", " GM     ",   ...           % Array of tissue names
                " CSF    ", " WML/MS ", " Fat    "  ];  %   (for formatted output)
%     T1 = split(TsStr); T1 = char(T1(:,:,3));            % Array of only the tissue type (3rd word) for legends,...
%     TisTag = string(T1(:,1:3,:));                       % ...trimmed to length 3 (e.g., "CSF" or "WML")
    TisTag = [  "WM ",  "GM ",                      ... % Explicit declaration for sanity
                "CSF",  "WML",  "Fat"               ];
    
            % 1.5 T (Visser); Visser 2010 MagResMed - used here as fallback T2 values (1.5 & 3.0 T)
            T2_WM  =   74.0;            % 1.5 T Wehrli, Yacoub
            T2_GM  =   87.0;            % 1.5 T Wehrli, Yacoub
            T2_CSF = 2280.0;            % 1.5 T Helms
            T2_MS  =  100.0;            % 1.5 T Alsop
            T2_Fat =   60.0;            % NB: PURE GUESSWORK! BEWARE!
    switch B0opt                        % Set T1 (and T2 if available) times based on B0 and sources
        case 1                  % 1.5 T (GE)
            iS.B0  =    1.5;
            T1_WM  =  660.0;            % 1.5 T value from GE 3dfse source code
            T1_GM  = 1200.0;            % 1.5 T --"--
            T1_CSF = 4270.0;            % 1.5 T --"--
            if ( iS.Vnd == "GE" ) && ( iPr.Mode ~= 1 )
                T1_CSF = 5100.0;        % GE sets a longer T1_CSF for DIR only! It probably compensates for Mz_end.
            end % if
            T1_MS  =  T1_GM;            % Estimate T1_WML ~ T1_GM
            T1_Fat =  192.0;            % 1.5 T value from GE 3dfse source code (250 ms for InHance only?)
            T2_WM  =   80.0;            % 1.5 T --"--
            T2_GM  =   95.0;            % 1.5 T --"--
            T2_CSF = 3500.0;            % 1.5 T --"--
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
        case 2                  % 3.0 T (GE)
            iS.B0  =    3.0;
            T1_WM  =  825.0;            % 3.0 T value from GE 3dfse source code
            T1_GM  = 1400.0;            % 3.0 T --"--
            T1_CSF = 4270.0;            % 3.0 T --"-- (Mathias Engström, GE)
            if ( iS.Vnd == "GE" ) && ( iPr.Mode ~= 1 )
                T1_CSF = 5000.0;        % GE sets a longer T1_CSF for DIR only! It probably compensates for Mz_end.
            end % if
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  250.0;            % 3.0 T value from GE 3dfse source code (300 ms for InHance only?)
            T2_WM  =   60.0;            % 3.0 T --"--
            T2_GM  =   70.0;            % 3.0 T --"--
            T2_CSF = 2400.0;            % 3.0 T --"-- (Mathias Engström, GE)
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
        case 3                  % 1.5 T (Custom): You can add your own times here.
            iS.B0  =    1.5;            % For display only
            T1_WM  =  660.0;            % 1.5 T GE source value = 660 ms
            T1_GM  = 1200.0;            % 1.5 T --"--           = 1200 ms
            T1_CSF = 4308.0;            % 1.5 T --"--           = 4270 ms; Siemens seem to use around 4308 ms?
            T1_MS  = 1150.0;            % Estimate T1_WML ~ T1_GM?
            T1_Fat =  192.0;            % 1.5 T GE source value = 192 ms
            
            T2_WM  =   80.0;            % 1.5 T --"--           =   80 ms
            T2_GM  =   95.0;            % 1.5 T --"--           =   95 ms
            T2_CSF = 3500.0;            % 1.5 T --"--           = 3500 ms
            T2_MS  =  100.0;            % 1.5 T Alsop           =   60 ms (meas. 100 ms @3T)
            T2_Fat =   60.0;            % NB: PURE GUESSWORK! BEWARE!
        case  4                 % 3.0 T (Custom): You can add your own times here.
            iS.B0  =    3.0;
            T1_WM  =  850.0;            % 3.0 T GE source value =  825 ms
            T1_GM  = 1400.0;            % 3.0 T --"--           = 1400 ms
            T1_CSF = 5000.0;            % 3.0 T --"--           = 4270 ms (5000 for DIR?!)
            T1_MS  = 1200.0;            % 3.0 T Alsop           = 1350 ms
            T1_Fat =  250.0;            % 3.0 T GE source value =  250 ms
            
            T2_WM  =   60.0;            % 3.0 T --"--           =   60 ms
            T2_GM  =   70.0;            % 3.0 T --"--           =   70 ms
            T2_CSF = 2400.0;            % 3.0 T --"--           = 2400 ms
            T2_MS  =  100.0;            % 1.5 T Alsop           =   60 ms (meas. 100 ms @3T)
            T2_Fat =   60.0;            % ???
        case 5                  % 1.5 T (Alsop); Madhuranthakam et al, Mag Res Med 2012 67(1):81-88
            iS.B0  =    1.5;
            T1_WM  =  650.0;            % 1.5 T Alsop
            T1_GM  = 1300.0;            % 1.5 T --"--
            T1_CSF = 4200.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ~ T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            T2_WM  =   70.0;            % 1.5 T Alsop
            T2_CSF = 2000.0;            % 1.5 T --"--
            T2_MS  =  100.0;            % 1.5 T --"--
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
            TisStr(5) = "(Fat)   ";     % (This value was cited from elsewhere)
        case 6                  % 3.0 T (Alsop)
            iS.B0  =    3.0;
            T1_WM  =  750.0;            % 3.0 T Alsop (other litt. has 850 ms or more?)
            T1_GM  = 1400.0;            % 3.0 T litteratur
            T1_CSF = 4200.0;            % 3.0 T litteratur
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean); Myrlund states 230 ms?
            TisStr(2) = "(GM)    ";     % (This value was cited from elsewhere)
            TisStr(5) = "(Fat)   ";     % (This value was cited from elsewhere)
        case 7                  % 1.5 T (MRIQ); MRIquestions.com and other online resources
            iS.B0  =    1.5;
            T1_WM  =  580.0;            % 1.5 T MRIQuestions
            T1_GM  =  940.0;            % 1.5 T --"--
            T1_CSF = 3600.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ~ T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
        case 8                  % 3.0 T (Lalande); Bojorquez et al, MRI 35(2017) 69-80
            iS.B0  =    3.0;
            T1_WM  =  940.0; %842.0;    % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR); Stikov (2015) has 940(860-992) ms.
            T1_GM  = 1425.0;            % 3.0 T Lalande (Liberman 2014 w/ VFA SPGR)
            T1_CSF = 4300.0;            % 3.0 T Lalande (trimmed mean of 4(6) reports)
            T1_MS  = 1350.0;            % 3.0 T Alsop
            T1_Fat =  385.0;            % 3.0 T Lalande (trimmed mean; 405.0 in Barral 2010 w/ SE-IR)
            T2_WM  =   75.0;            % 3.0 T (Lu 2005; jMRI 22(1):13-22)
            T2_GM  =   83.0;            % 3.0 T Lalande
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
        case 9                  % 1.5 T (Visser); Visser 2010 MagResMed
            iS.B0  =    1.5;
            T1_WM  =  656.0;            % 1.5 T Rooney
            T1_GM  = 1188.0;            % 1.5 T --"--
            T1_CSF = 4329.0;            % 1.5 T --"--
            T1_MS  =  T1_GM;            % Estimate T1_WML ~ T1_GM
            T1_Fat =  260.0;            % 1.5 T MRIQuestions
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
            TisStr(5) = "(Fat)   ";     % (This value was cited from elsewhere)
        case 10                 % 7.0 T (Visser); Visser 2010 MagResMed
            iS.B0  =    7.0;
            T1_WM  = 1220.0;            % 7.0 T Rooney
            T1_GM  = 2132.0;            % 7.0 T --"--
            T1_CSF = 4400.0;            % Estimated from the 1.5 T value...
            T1_MS  =  T1_GM;            % Estimate T1_WML ~ T1_GM
            T1_Fat =  600.0;            % 7.0 T estimated(?!?)
            T2_WM  =   46.0;            % 7.0 T Wehrli, Yacoub (= a factor 0.63 from 1.5 T)
            T2_GM  =   55.0;            % 7.0 T Wehrli, Yacoub (= a factor 0.63 from 1.5 T)
            T2_CSF = 1800.0;            % Estimated from the 1.5 T value (> factor 0.63...)
            T2_MS  =   70.0;            % Estimated from the 1.5 T value (> factor 0.63...)
            TisStr(3) = "(CSF)   ";     % (This value was cited from elsewhere)
            TisStr(4) = "(MS/WML)";     % (This value was cited from elsewhere)
            TisStr(5) = "(Fat)   ";     % (This value was cited from elsewhere)
        case 11                 % 3.0 T (meas.); 2019-02-06 T1/T2 brain measurements from OUS UUS NMR3 (by WibeN)
            iS.B0  =    3.0;
            T1_WM  =  850.0;            % 3.0 T value measured on brain at UUS NMR3 w/ multi-IR FSE (GE: 825 ms; meas. 840-900)
            T1_GM  = 1420.0;            % 3.0 T --"-- (GE: 1400 ms)
            T1_CSF = 4270.0;            % 3.0 T value from GE 3dfse source code
            if ( iS.Vnd == "GE" ) && ( iPr.Mode ~= 1 )
                T1_CSF = 5000.0;        % GE sets a longer T1_CSF for DIR only! It probably compensates for Mz_end.
            end % if
            T1_MS  = 1050.0;            % 3.0 T Measured as above (Alsop: 1350 ms; measured 920 & 1050 in two lesions)
            T1_Fat =  250.0;            % 3.0 T value from GE 3dfse source code (300 ms for InHance only?)
            T2_WM  =   55.0;    60.0;   % 3.0 T value measured on brain at UUS NMR3 w/ multi-TE SE (GE: 60 ms; meas. 55 ms)
            T2_GM  =   65.0;    70.0;   % 3.0 T --"-- (GE: 70 ms; meas. 58 ms but everyone has +10 from WM? But, PV from CSF?)
            T2_MS  =  102.0;            % 3.0 T --"-- (measured in one lesion in one healthy volunteer only)
            T2_CSF = 2400.0;            % 3.0 T value from GE 3dfse source code
            TisStr(3) = "(CSF)   ";     % (This value was cited from elsewhere)
            TisStr(5) = "(Fat)   ";     % (This value was cited from elsewhere)
        otherwise
            error('ERROR: Undefined B0/source!');
    end % switch B0opt
    
    iPr.B0tag = sprintf( '%1.1f T', iS.B0 );        % Char array of the field strength in T, for figure titles
    Ts.T  = [ [ T1_WM, T1_GM, T1_CSF, T1_MS, T1_Fat ];  ...
              [ T2_WM, T2_GM, T2_CSF, T2_MS, T2_Fat ] ];
    Ts.Tag  = [ TisTag; TisStr ];

end % script fn
