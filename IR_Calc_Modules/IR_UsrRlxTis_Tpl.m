function IR_UserRelaxTimes()                    % Set user defined relax. times for the IR_SetRelaxTimes module.
%% User defined (D)IR MRI T1/T2 settings for IR_Calculator. WM/GM/CSF, MS(WML) and Fat.
%   - These values are read at startup by the IR_SetRelaxTimes() IR-Calc module.                     -->
%   - The user values can be selected in the IR-Calculator as "T1s" option 3 or 4.                   -->
%   - If you change them while IR-Calc is running, press 'r' to refresh the values.                  -->
%
%  TODO:
%   - Migrate to .xml format for the MatLab R2020b+ readstruct() fn.
%

%% USERDATA

global iC                                               % IR-Calc globally used data structs

% User1 settings:                               % 1.5*T (Custom): You can add your own times here.
    iC.S.U1.B0sys   =    1.5                        ;   % #.# T (Used for display only)

    iC.S.U1.T1_WM   =  660.0                        ;   % 1.5 T GE source value = 660 ms
    iC.S.U1.T1_GM   = 1200.0                        ;   % 1.5 T --"--           = 1200 ms
    iC.S.U1.T1_CSF  = 4308.0                        ;   % 1.5 T --"--           = 4270 ms; Siemens seem to use around 4308 ms?
    iC.S.U1.T1_MS   = 1200.0                        ;   % Estimate T1_WML ~ T1_GM?
    iC.S.U1.T1_Fat  =  192.0                        ;   % 1.5 T GE source value = 192 ms

    iC.S.U1.T2_WM   =   80.0                        ;   % 1.5 T --"--           =   80 ms
    iC.S.U1.T2_GM   =   95.0                        ;   % 1.5 T --"--           =   95 ms% 
    iC.S.U1.T2_CSF  = 3500.0                        ;   % 1.5 T --"--           = 3500 ms
    iC.S.U1.T2_MS   =  100.0                        ;   % 1.5 T Alsop           =   60 ms (meas. 100 ms @3T)
    iC.S.U1.T2_Fat  =   60.0                        ;   % ???

% User2 settings:                               % 3.0*T (Custom): You can add your own times here.
    iC.S.U2.B0sys   =    3.0                        ;   % #.# T (Used for display only)

    iC.S.U2.T1_WM   =  850.0                        ;   % 3.0 T GE source value =  825 ms
    iC.S.U2.T1_GM   = 1400.0                        ;   % 3.0 T --"--           = 1400 ms
    iC.S.U2.T1_CSF  = 5000.0                        ;   % 3.0 T --"--           = 4270 ms (5000 for DIR?!)
    iC.S.U2.T1_MS   = 1400.0                        ;   % 3.0 T Alsop           = 1350 ms
    iC.S.U2.T1_Fat  =  250.0                        ;   % 3.0 T GE source value =  250 ms

    iC.S.U2.T2_WM   =   60.0                        ;   % 3.0 T --"--           =   60 ms
    iC.S.U2.T2_GM   =   70.0                        ;   % 3.0 T --"--           =   70 ms
    iC.S.U2.T2_CSF  = 2400.0                        ;   % 3.0 T --"--           = 2400 ms
    iC.S.U2.T2_MS   =  100.0                        ;   % 1.5 T Alsop           =   60 ms (meas. 100 ms @3T)
    iC.S.U2.T2_Fat  =   60.0                        ;   % ???

end % script fn
