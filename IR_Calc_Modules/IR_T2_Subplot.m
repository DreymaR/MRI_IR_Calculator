function Fig2 = IR_T2_Subplot( iT, iS2 )    % S0/T2/tags for tissues to plot; find TEopt for two ( S0set, T2set, TsTag, TEset )
%% (D)IR MRI TE calculation and T2 decay plot for IR_Calc
%   - This function is called from the IR_TI_Calculator
%   - Plot the T2 decay at readout, to determine the optimal TE.
%
%  NOTE:
%   - Use symbolic solution instead of just the graphical one? Not critical.
%
%  TODO:
%   - Add max TE for GM/WML contrast (for WM-DIR)
%   - Find TE range for 90% (80%?) DelS? Find the point under the threshold that has a neighbor above? Or solve?
%   - Make the TE_opt calculation more robust for different cases?
%   - Make TE_set more flexible. Allow unsetting it, updating for each replot. Maybe a "Freeze TE" button instead?
%       - Best to start it @ 0 and ignore it then? But then we can't compare T=0. Set it to '--'?
%       - Make it so a click on the text resets it? Secret mode, or in mouse-over text?
%   - Make the T2 plot a subplot instead of a new figure!? Send this script our main figure handle then.
%       - Specify subplot properly, like [ ax = subplot(1, 2, 1, 'Parent', fig); ]?
%       - Must then be specific about which axes I'm drawing to in each case!? Set current etc.
%       - But it seems that subplots are intended to be in panels in the same window? So, less attractive?
%
%  DONE:
%   - Make a T2-TE figure that can be updated from the GUI panel for the current signal/plot settings
%   - Calculate relative SNR (abs. S(TE) for WML) and CNR (abs. deltaS(TE) for WML-WM)
%   - CNR/t (vs TR)! Could use something like seconds(duration([0,5,57])) (but would like duration as "mm:ss")?
%       - Simply though, just use (SNR1/SNR2)/sqrt(TR1/TR2) as a measure of relative efficiencies, as in literature
%   - A freeze button for RelS/C so you can do multiple comparisons to the same standard setting
%   - Use less globals! Conflict with other script when names are the same. Structs are economic.
%   - Report contrast between WML and both secondary tissues (WM and GM). Always relevant?
%


%% INIT
debugInfo  = 1;                                     % Show extra info in the command window?
Tis     = iT.T2tis;                                 % Tissues to T2 plot: Normally WML, WM, GM
% nTis    = size( Tis, 2 );                           % # of tissues to compare to the first one, including itself
DUR     = 499   ; % iS2.ETD;                        % Use the sequence Echo Train Duration as plot duration? No, less.

global iC  iPr                                      % IR-Calc Program data struct (common global with Fig1)
global UI2 iR2                                      % UI/tissue/times data structs (sequence struct is iS2)
global iC_Ini2 UI2Bx UI2on frzS2                    % Non-struct globals
global oldS2 oldC2 oldTR  TEopt TEset TEmin TEstr oTs1 oTs2     % Result globals (TODO: Put some in structs)

fSz = 10.80; % sSz(4)/100.00;                       % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])    % Default fig. position & size
figName = 'MR T2 contrast calculator';
iPr.Fig2 = findobj( 'Type', 'Figure',           ... % Find this script's figure (possibly among others)
        'Name', figName                         );
if isempty(iPr.Fig2)                                % Iff my figure does not yet exist... (was ~ishghandle)
    iPr.Fig2 = figure(                          ... % ...make a new figure!
        'Name', figName,                        ...
        'Color', iPr.FigClr,                    ...
        'KeyPressFcn', @keyPress_callback,      ...
        'Position', [ 14 16 80 60 ]*fSz         );  % At pos. (x,y)% of ScreenSize, etc.
else
    figure(iPr.Fig2);                               % If my figure exists, make it the CurrentFigure
end % if
Fig2 = iPr.Fig2;                                    % Return the handle to this figure

t = sym('t');                                       % Set the symbolic explicitly as it's used in nested functions

if isempty( iC_Ini2 )
    iC_Ini2  = true ;                               % The program has been initialized
    UI2on   = true  ;                               % Show the UI control panel
    frzS2   = 0     ;                               % Don't freeze the RelSNR at first
    TEmin   = 14    ; % ms                          % Minimum allowed TE
    TEset   = -1.0  ;                               % The manually set TE (will be set to TEopt)
end % if iC_Ini2

%% MAIN
main();
function main()
    clf                                             % Clear the figure so we can change it
    mainTEcalc();
    createUI();
end % main

% -) Calculate & plot MRI S(t) over a duration, and TE giving optimal contrast
function mainTEcalc()
    oTs1  = Tis(1);    oTs2  = Tis(2);              % Usually WML, WM, GM
    if ( abs(03.5*iT.S0(oTs2)) < abs(iT.S0(oTs1)) ) % If 2nd signal (WM) is very low (nulled)...
        oTs2 = Tis(3);                              % ...use the 3rd tissue (GM) instead (TODO: Use both and compare?)
    end % if
    NTs = length(iT.S0);
    St = zeros(DUR+1,NTs);                          % Array of signals
    
    St(1,:) = iT.S0;
    for tpt = 1:DUR
        for i = 1:NTs
            St(tpt+1,i) = iT.S0(i).*exp(-tpt./iT.T2(i));
        end % for i
    end % for tpt
    
    % TODO: Calculate both for Ts1 vs Ts2 (WML vs WM) and Ts1 vs Ts3 (WML vs GM)?!
    S_dif = ( St(:,oTs1) - St(:,oTs2) );            % Calculate signal difference between two T2s; don't use Abs() here
    [ dif_max, maxind ] = max(S_dif)        ;       % Find the max difference
    [ abs_max, absind ] = max(abs(S_dif))   ;       % Find the max difference of the absolute function
    [ neg_max, negind ] = max(-S_dif)       ;       % Find the max difference of the negative function
    if ( debugInfo )
        fprintf("   "); fprintf("%-6s", ["pos","abs","neg","t (ms)"] ); fprintf("\n");     % Legends
        disp( [ maxind, absind, negind ]-1 )        % -1 converts indices to times (ms)
    end % if debug
    S0_1 = iT.S0(Tis(1)); S0_2 = iT.S0(Tis(2));
    if ( S0_1 < 0 ) && ( S0_2 < 0 )                 % TODO: Is this robust?
        S_dif = -S_dif;                             %       (abs_max > dif_max) || (maxind == 1) && (absind ~= 1) ???
        dif_max = neg_max; maxind = negind;         % abs_max; absind;
    end % if
    TEopt = maxind -1;
    [       ~, S_atTE ] = max(abs(St(maxind,:) ));  % The strongest signal at the time of max difference
%     TI2_n = uint16(FindOptTE(                       % TODO: Solve symbolically?
    
    if ( TEset < 0 )
        TEset = TEopt;
    end % if
    
    hold on
    plot( St,                                   ... % Plot the signals calculated above
        'LineWidth',    1.0                     )
    plot( S_dif,        'b-' ,                  ...
        'LineWidth',    2.0                     )
    plot( zeros(1,DUR), ':',                    ... % Add a zero line
        'Color',        'black',                ...
        'LineWidth',    1.0                     )
    hold off
    
%     whitebg( 'white' );
%     whitebg( [ 0.70 1.00 0.70 ] );                  % Debug: See white plot/GUI elements
    Tag1 = regexp(iT.TsTag(oTs1),' ','split');      % Trick to trim space chars off single-word strings
    Tag2 = regexp(iT.TsTag(oTs2),' ','split');      % TODO: Use 'tokens' instead for sanity? Trickyish though.
    Tag3 = num2str(iS2.TRef);
%     txt = join([ "T2 decay and TE_{opt} (" Tag1 " vs " Tag2 ")" ],"");
    txt = join([ "T2 decay and TE_{opt} (" Tag1 " vs " Tag2 "), for TR_{eff} = " Tag3 " ms"],"");
    title(  txt,                                ...
            'FontSize', 16                      );
    ylabel( 'Rel. MR signal', 'FontSize', 12 )
    xlabel( 'time (ms)'     , 'FontSize', 12 )
%     xlm = [ 0 200*ceil( TI2_max/200 ) ]; xlim(xlm); % Round up to 200 on the x axis
%     ylm = [ -1 1 ];                      ylim(ylm); % Rel. S is in [-1,1]
%     ylm = ylim();
    
    if ( TEopt == DUR )
        TEstr = "N/A";
    elseif (TEopt == 0 )
        TEopt = TEmin;
        TEstr = string([ '(~' num2str(TEopt) ' ms)' ]);
    else
        TEstr = string([ num2str(TEopt) ' ms' ]);
    end % if
    
    lxlm = [ TEopt            , TEopt   ];          % Time of max signal difference
    lylm = [ St(maxind,S_atTE), dif_max ];          % Lower/upper y of line (= S_dif curve)
    line( lxlm, lylm,                           ... % Vertical line @ crossing
        'Color', 'black',                       ...
        'LineWidth', 1.5,                       ...
        'DisplayName', 'TE_{opt}',              ... % (Tip: If declared after legend, the line appears in it)
        'LineStyle', ':'                        );
    
    LegS = strings(1,NTs);                          % Make a legend...
    for i = 1:NTs                                   % ...showing the tissue short names and T2 times
        LegS(i) = sprintf('%s (T2 = %i ms)', iT.TsTag(i), iT.T2(i) );
    end % for i
    legend( LegS,                               ...
            'Location', 'South'                 );  % 'Best' may conflict with UI and TextBox
    
    iR2.magS    = abs( St(TEopt+1, oTs1) );
    iR2.magC    = abs( dif_max           );         % TODO: How does this work? It should be negative if WML<WM
    iR2.magCW   = magC( St(TEopt+1, Tis(1)) , St(TEopt+1, Tis(2)) );    % WML-WM contrast
    iR2.magCG   = magC( St(TEopt+1, Tis(1)) , St(TEopt+1, Tis(3)) );    % WML-GM contrast
    if isempty(oldS2)   ; oldS2 = iR2.magS  ;   end % Old signal strength   for comparison; as float
    if isempty(oldC2)   ; oldC2 = iR2.magC  ;   end % Old signal difference for comparison; as float
    if isempty(oldTR)   ; oldTR = iS2.TR    ;   end % Old rep. time         for comparison; as float
    relS    = uint16(100*(iR2.magS/oldS2));         % Rel. Signal for WML
    relC    = uint16(100*(iR2.magC/oldC2));         % Rel. Contrast between WML and tissue
    relS_t  = uint16(100*( (iR2.magS/oldS2)/sqrt(iS2.TR/oldTR) ));  % Rel. SNR per time (ratio of S over TR^0.5)
    relC_t  = uint16(100*( (iR2.magC/oldC2)/sqrt(iS2.TR/oldTR) ));  % Rel. CNR per time (ratio of S over TR^0.5)
    fStr = ['\\itOptimal TE:\n\\rm\\bf'         ... % Tex \bf\it\rm = bold/italic/normal; escape \ for sprintf!
            ' TE_{opt}    = %3s \n'             ... % %3.i ms
            ' S_{WML}    = %4.3f \n'            ... % %4.3f
            ' \\DeltaS_{max}  = %4.3f \n'       ... % %4.3f
            ' rS_{WML}   = %3.f%%\n'            ... % %4.3f
            ' rC(L-?)  = %3.f%%\n'              ... % %4.3f
            ' rSNR/t   = %3.f%%\n'              ... % %4.3f
            ' rCNR/t   = %3.f%%'                ];  % %4.3f
    fVar = [ TEstr, iR2.magS, iR2.magC, relS, relC, relS_t, relC_t ]; % 
    figTextBox( fStr, fVar );                       % Display an info text box on the figure
    if ( frzS2 == 0 )                               % If RelZ isn't frozen (by the UI button)
        oldC2 = iR2.magC;
        oldS2 = iR2.magS;
        oldTR = iS2.TR;
    end % if frzS2
    
    if ( TEset ~= TEopt )
        iR2.magSset = abs( St(TEset+1, oTs1) );
        iR2.dirCset = magC( St(TEset+1, oTs1  ) , St(TEset+1, oTs2  ) );
        iR2.C_WMset = magC( St(TEset+1, Tis(1)) , St(TEset+1, Tis(2)) );    % WML-WM contrast
        iR2.C_GMset = magC( St(TEset+1, Tis(1)) , St(TEset+1, Tis(3)) );    % WML-GM contrast
    end % if

end % mainTEcalc

function figTextBox( fStr, fVar )
    txt = sprintf( fStr, fVar );                    % Formatted string "print"
%     lin = count( fStr, '\n' ) + 1;                  % Count lines in the string
    xlm = xlim; xt =  0.780*double(xlm(2));         % Max value of axis; (x,y) is a plot point!
    ylm = ylim; yt = double( (ylm(2)+ylm(1))/2. );  % Place the text box at the middle of the y axis
%     yt = yt + 0.035*lin;                            % Compensate for number of lines
    text( xt, yt, txt,                          ... % Note: annotation('textbox' etc doesn't scale with text
        'EdgeColor', 'black',                   ...
        'Margin', 4,                            ...
        'LineWidth', 1,                         ...
        'LineStyle', '-'                        );
end % fcn

%% FUNCTIONS
function FltOrSym = FoS( InNr )                     % Cast as double for, e.g., exp() unless it's used as a sym
    if      isa( InNr, 'sym' ), FltOrSym = InNr;
    else,                       FltOrSym = double(InNr);
    end % if
end % fcn

function St = CalcSt( S0, T2, Ts )                  % Calculate generic S(t) = S0*exp(-t/T2)
    S0 = FoS(S0); T2 = FoS(T2); Ts = FoS(Ts);       % Cast input as float unless it's a sym
    St = S0*exp(-Ts./T2);                           % One signal value per given T2
end % fcn

function TEm = FindOptTE( S0_1, S0_2, T2_1, T2_2 )  % Find TE to give max contrast between two decaying tissues
    T2_1 = FoS(T2_1); T2_2 = FoS(T2_2);             % Cast input as float unless it's a sym
    % TODO: How to use vpasolve to find the max?
    TEm = vpasolve( ( CalcSt(S0_1,T2_1) ...         % Solve |S1(t) - S2(t)| == max for t = TE_opt
                    - CalcSt(S0_2,T2_2) ), t);      % Return TE_opt
end % fcn

function mC = magC( SL, SB )                        % Calculate contrast between "lesion" and "background" magnitude S
    mC = abs( SL ) - abs( SB );                     % Use the difference between abs. signals as each voxel is mag.
end % fcn

%% PLOT

% Plot S(t) under T2 decay
% function PlotSt( T2plt, DUR )                       % DUR = Total duration for the plot
% end % fcn

%% USER INTERFACE

function createUI()                                 % Display a button to show/hide the UI controls
    eW = 0.065 ; eH = 0.065; eS = eH + 0.025;       % Button width/height/spc for masterUI controls
    myAx = gca;                                     % A handle to the current graphics axes (plot)
    uiX = myAx.OuterPosition(3) - 0.015 - eW;       % Use the axis' outer x position (def. [ 0 0 1 1 ]) for UI pos.
    uiY = myAx.Position(4)      + 0.010     ;       % Use the axis' inner y position to determine UI pos.
    
    uicontrol( 'Style'      ,   'pushbutton',   ... % For all modes:
        'String'            ,   'Panel'     ,   ...
        'FontWeight'        ,   'bold'      ,   ...
        'ToolTipString', '[h] Show/hide UI',    ...
        'Units'             ,   'normalized',   ...
        'Position', [ uiX uiY-0*eS eW eH ]  ,   ...
        'Callback',         @masterUI_callback  );  % UI to show/hide the UI control panel
    
    createUIPanel();
    switch UI2on
        case true ; UI2Bx.Visible = 'on';
        case false; UI2Bx.Visible = 'off';
    end % switch UI2on
end % createUI

function createUIPanel()
    boxClr = iPr.FigClr;                            % UI panel background color (same as fig background)
    mX = 10; mY =  6;                               % x/y margins in UI panel
    eN = 04; eW = 130; eH = 30;                     % Number of UI elements; default element width/height
    UI2.Hpx = eN*eH + 2*mY; uiT = UI2.Hpx - mY + 3; % UI panel height in pixels; top position for elements
    UI2.Wpx = eW + 2*mX;                            % UI panel width in pixels
    
    UI2Bx = uipanel(            iPr.Fig2    ,   ... % UI2Bx
        'BackgroundColor'   ,   boxClr      ,   ... %  [ 0.8 0.8 0.4 ]
        'BorderType'        ,   'etchedin'  ,   ... % default 'etchedin'
        'BorderWidth'       ,   1.0         ,   ... % default 1
        'Clipping'          ,   'off'       ,   ... % default 'on'; 'off' is good for debugging
        'Units'             ,   'normalized',   ... % Normalized units facilitate resizing to axes
        'ResizeFcn'         ,   @placeUiBox     );  % Resize function for the UI panel
    placeUiBox();                                   % Adjust the UI panel position to the figure window
    
    for i = 1:eN; eY = uiT - i*eH;                  % Generate each UI element on a new row in the panel
        switch i
        case  1
        uicontrol(  'Parent',   UI2Bx           ,   ...
           'Position'   , [ mX+00 eY+00 eW 24 ] ,   ...
            'Style'         ,   'pushbutton'    ,   ...
            'String'        ,   'Results'       ,   ...
            'ToolTipString' ,   '[p] Print results',...
            'Callback'      ,   @printVars_callback );  % UI to print info/results to the command window
%             'Units'         ,   'normalized'    ,   ...
        
        case  2
        uicontrol(              'Parent', UI2Bx ,   ... % ui TE[T|S]
            'Position',   [ mX+00 eY+01 eW 16 ] ,   ...
            'Style'         ,   'text'          ,   ...
            'BackgroundColor',  boxClr          ,   ...
            'HorizontalAlignment', 'Left',          ...
            'ToolTipString' ,   '[q] Echo time' ,   ...
            'String', 'TE_set:                      ms');
        UI2.TES = uicontrol(    'Parent', UI2Bx ,   ...
            'Position',   [ mX+40 eY+00 60 20 ] ,   ...
            'Style'         ,   'edit'          ,   ...
            'Min', 50000    ,   'Max', 50000    ,   ...
            'ToolTipString', 'Can be set manually', ...
            'String'        ,   TEset           ,   ...
            'Value'         ,   TEset           ,   ...
            'Callback', @te_set_callback            );  % UI to set TR

        case  3
        uicontrol(              'Parent', UI2Bx ,   ... % uiSetS
           'Position',    [ mX+10 eY+00 eW 24 ] ,   ...
            'Style'         ,   'radiobutton'   ,   ... % 'togglebutton'/'checkbox'
            'BackgroundColor',  boxClr          ,   ...
            'Value'         ,   frzS2           ,   ...
            'String'        ,   'Freeze relSNR' ,   ...
            'ToolTipString', '[r] Set comparison std.', ...
            'Callback'      ,   @setSignal_callback );  % UI to set the relative signal for comparison
%             'Min',  0,          'Max',  1,          ... % 0/1 are default for Matlab toggles

        case  4
        uicontrol(              'Parent', UI2Bx ,   ... % ui DUR[T|S]
            'Position',   [ mX+00 eY+01 eW 16 ] ,   ...
            'Style'         ,   'text'          ,   ...
            'BackgroundColor',  boxClr          ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'ToolTipString', '[d] Plot duration',   ...
            'String', 'Dur.:                          ms');
        UI2.DURS = uicontrol(   'Parent', UI2Bx ,   ...
            'Position',   [ mX+40 eY+00 60 20 ] ,   ...
            'Style'         ,   'edit'          ,   ...
            'Min', 50000    ,   'Max', 50000    ,   ...
            'ToolTipString', 'Can be set manually', ...
            'String'        ,   DUR             ,   ...
            'Value'         ,   DUR             ,   ...
            'Callback'      ,   @dur_set_callback   );  % UI to set DUR
        
        end % switch i
    end % for i
end % createUIPanel

function masterUI_callback(~,~)                     % UI that shows/hides the UI control panel
    UI2on = ~UI2on;
    switch UI2on
        case true ; UI2Bx.Visible = 'on';
        case false; UI2Bx.Visible = 'off';
    end
end % fcn

function printVars_callback(~,~)                    % UI that prints info to the command window
    fStr = ['\n'                                ... % TODO: Isn't fprintf obeying the field width operator for strings?!
            '*******************************\n' ...
            '***    T2 plot results:     ***\n' ...
            '*******************************\n' ...
            '*  S0_%3s         =  % 1.3f   *\n' ... % TODO: Align when positive?
            '*  S0_%3s         =  % 1.3f   *\n' ...
            '*  S0_%3s         =  % 1.3f   *\n' ...
            '*  S0_%3s         =  % 1.3f   *\n' ...
            '*******************************\n' ...
            '*  TE_opt         =  %+7s','  *\n' ...
            '*  S_WML @TE_opt  =  % 1.3f   *\n' ...
            '*  C_WM   --"--   =  % 1.3f   *\n' ...
            '*  C_GM   --"--   =  % 1.3f   *\n' ];  % TODO: Proper right alignment of positive/negative numbers
    fVar = { iT.TsTag(1),iT.S0(1), iT.TsTag(2),iT.S0(2),    ...
             iT.TsTag(3),iT.S0(3), iT.TsTag(4),iT.S0(4),    ...
             TEstr, iR2.magS, iR2.magCW, iR2.magCG          };
    fprintf( fStr, fVar{:} );
    if ( TEset ~= TEopt )
        fStr = [                                    ...
                '*******************************\n' ...
                '*  TE_set         =  %4.i ms  *\n' ...
                '*  S_WML @TE_set  =  % 1.3f   *\n' ...
                '*  C_WM   --"--   =  % 1.3f   *\n' ...
                '*  C_GM   --"--   =  % 1.3f   *\n' ...
                '*  relS(set/opt)  =  %4.1i %%   *\n' ...
                '*  relC(set/opt)  =  %4.1i %%   *\n' ];
        fVar = { uint16(TEset), iR2.magSset, iR2.C_WMset, iR2.C_GMset,              ...
                 uint16(100*iR2.magSset/iR2.magS), uint16(100*iR2.dirCset/iR2.magC) };
        fprintf( fStr, fVar{:} );
%                 '*******************************\n' ...
%                 '***     Manual echo time:   ***\n' ...
    end % if
    fprintf(    '*******************************\n' );
end % fcn

function te_set_callback(src,~)                     % UI to set TE manually
    val = str2val( src.String );
    if ( ( val >= 0 ) && ( val < iS2.ETD ) )        % Allowed range
        TEset = val;
    else
        TEset = TEopt;                              % You can reset with a neg. value
    end % if
    main();
end % fcn

function setSignal_callback(~,~)                    % UI callback to toggle Freeze RelSNR
    frzS2 = ~frzS2;                                 % was src.Value;
    main();
end % fcn

function dur_set_callback(src,~)                    % UI to set DUR manually
    val = str2val( src.String );
    if ( ( val > TEopt ) && ( val < iS2.ETD ) )     % Allowed range
        DUR = val;
    else
        DUR = iS2.ETD;                              % You can reset with a neg. value
    end % if
    main();
end % fcn

function keyPress_callback(~,evt)
    switch evt.Key
        case 'h'                                    % Show/hide UI panel
            masterUI_callback([],[])                % Calling with dummy values
        case 'p'                                    % Print results
            printVars_callback([],[])
        case 'q'                                    % Focus on TE setting uicontrol
            uicontrol( UI2.TES );
        case 'd'                                    % Focus on DUR setting uicontrol
            uicontrol( UI2.DURS );
        case 'a'                                    % Change active figure
            if isfield( iPr, 'Fig1' )
                figure( iPr.Fig1    );
            end % if
    end % switch
end % fcn

function placeUiBox(~,~)                            % Function to resize the UI panel (automatically called)
    if ishandle(UI2Bx)                              % (The first time this fn is called, the panel doesn't exist yet)
        pxFig = getpixelposition(iPr.Fig2);         % Figure size in px, for manual normalizing
        myAx = gca;                                 % A handle to the current graphics axes/plot
        uiW = UI2.Wpx / pxFig(3);                   % Normalized UI panel width
        uiH = UI2.Hpx / pxFig(4);                   % Normalized UI panel height
        uiX = myAx.Position(3) - uiW + 0.115;       % The axis' inner right x pos. (but why the "fudge factor"?)
        uiY = myAx.Position(4) - uiH + 0.090;       % The axis' inner upper y pos. (--"--)
        UI2Bx.Position = [ uiX uiY uiW uiH ];
    end % if
end % fcn

function val = str2val( str )                       % Function to process numeric input for UI controls
    val = str2double( str );
    if ( val ~= val )                               % test for NaN from malformed input
        val = -1;
    end % if
end % fcn

end % script fn