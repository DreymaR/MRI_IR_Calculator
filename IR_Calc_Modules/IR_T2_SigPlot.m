function IR_T2_SigPlot( iT )                    % S0/T2/tags for tissues to plot; find TEopt for two ( S0set, T2set, TsTag, TEset )
%% (D)IR MRI TE calculation and T2 decay plot for IR_Calc
%   - This function is called from the IR_TI_Calculator
%   - Plot the T2 decay at readout, to determine the optimal TE.
%
%  NOTE:
%   - 
%
%  TODO:
%   - Simulate echo train by the Hennig phase evolution method?!?
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
%  ONHOLD:
%   - Use symbolic solution instead of just the graphical one? Not critical.
%


%% INIT
debugInfo  = 1                                      ;   % Show extra info in the command window?
Tis     = iT.T2tis                                  ;   % Tissues to T2 plot: Normally WML, WM, GM
% nTis    = size( Tis , 2 )                         ;   % # of tissues to compare to the first one, including itself
DUR     = 499   ;   % iC.S.ETD                      ;   % Use the sequence Echo Train Duration as plot duration? No, less.

global iC                                               % IR-Calc globally used data structs
global iC_Ini2 F2_Pan UI2                               % Non-struct globals; UI for Fig 2

fSz = 10.80 ;   % sSz(4)/100.00                     ;   % Scaling factor for figures, ?1.0% of screen size
% set(0,'defaultfigureposition',[3*fsz 1*fsz 4*fsz 3*fsz])  % Default fig. position & size
figName = 'MR T2 contrast calculator'               ;
iC.P.Fig2 = findobj( 'Type' , 'Figure'          ,   ... % Find this script's figure (possibly among others)
                    'Name'  , figName           )   ;
if isempty(iC.P.Fig2)                                   % Iff my figure does not yet exist... (was ~ishghandle)
    iC.P.Fig2 = figure(                             ... % ...make a new figure!
        'Name'          , figName               ,   ...
        'Color'         , iC.P.FigClr           ,   ...
        'Position'      , [ 14 16 80 60 ]*fSz   ,   ...  % At pos. (x,y)% of ScreenSize, etc.
        'KeyPressFcn'   , @keyPress_callback    )   ;
else
    figure(iC.P.Fig2)                               ;   % If my figure exists, make it the CurrentFigure
end % if

t = sym('t')                                        ;   % Set the symbolic explicitly as it's used in nested functions

if isempty( iC_Ini2 )                                   % Parameter presets for the first run only
    UI2.PanVis  = true      ;                           % Show the UI control panel
    UI2.FrzSNR  = false     ;                           % Don't freeze the RelSNR at first
    iC.TEmin    = 14        ; % ms                      % Minimum allowed TE
    iC_Ini2     = true      ;                           % The program has been initialized
end % if iC_Ini2

%% MAIN
main();
function main()
    clf                                                 % Clear the figure so we can change it
    mainTEcalc()                                    ;
    createUI()                                      ;
end % main

% -) Calculate & plot MRI S(t) over a duration, and TE giving optimal contrast
function mainTEcalc()
    iC.M.o1 = Tis(1) ; iC.M.o2 = Tis(2)             ;   % Usually WML, WM, GM
    if ( abs(03.5*iT.S0(iC.M.o2)) < abs(iT.S0(iC.M.o1)) )   % If 2nd signal (WM) is very low (nulled)...
        iC.M.o2 = Tis(3)                            ;       % ...use the 3rd tissue (GM) instead (TODO: Use both and compare?)
    end % if
    NTs = length(iT.S0)                             ;
    St = zeros(DUR+1,NTs)                           ;   % Array of signals
    
    St(1,:) = iT.S0                                 ;
    for tpt = 1:DUR
        for i = 1:NTs
            St(tpt+1,i) = iT.S0(i).*exp(-tpt./iT.T2(i)) ;
        end % for i
    end % for tpt
    
    % TODO: Calculate both for Ts1 vs Ts2 (WML vs WM) and Ts1 vs Ts3 (WML vs GM)?!
    S_dif = ( St(:,iC.M.o1) - St(:,iC.M.o2) )       ;   % Calculate signal difference between two T2s; don't use Abs() here
    [ dif_max , maxind ] = max(S_dif)               ;   % Find the max difference
    [ abs_max , absind ] = max(abs(S_dif))          ;   % Find the max difference of the absolute function
    [ neg_max , negind ] = max(-S_dif)              ;   % Find the max difference of the negative function
    if ( debugInfo )
        fprintf("   ") ; fprintf("%-6s", ["pos","abs","neg","t (ms)"] ) ; fprintf("\n") ;   % Legends
        disp( [ maxind , absind , negind ] -1 )         % -1 converts indices to times (ms)
    end % if debug
    S0_1 = iT.S0(Tis(1)); S0_2 = iT.S0(Tis(2))      ;
    if ( S0_1 < 0 ) && ( S0_2 < 0 )                     % TODO: Is this robust?
        S_dif = -S_dif                              ;   %       (abs_max > dif_max) || (maxind == 1) && (absind ~= 1) ???
        dif_max = neg_max; maxind = negind          ;   % abs_max; absind;
    end % if
    iC.R.TEopt = maxind -1                          ;
    [ ~ , S_atTE ] = max(abs( St(maxind,:) ))       ;   % The strongest signal at the time of max difference
    
    if ~isfield( UI2, 'TEset' )
        UI2.TEset = iC.R.TEopt                      ;
    end % if
    
    hold on
    plot( St            ,                           ... % Plot the signals calculated above
        'LineWidth'     , 1.0                   )
    plot( S_dif         , 'b-'                  ,   ...
        'LineWidth'     , 2.0                   )
    plot( zeros(1,DUR)  , ':'                   ,   ... % Add a zero line
        'Color'         , 'black'               ,   ...
        'LineWidth'     , 1.0                   )
    hold off
    
%   whitebg( 'white' )                              ;   % DEBUG
%   whitebg( [ 0.70 1.00 0.70 ] )                   ;   % DEBUG: See white plot/GUI elements
    Tag1 = regexp(iT.TsTag(iC.M.o1),' ','split')    ;   % Trick to trim space chars off single-word strings
    Tag2 = regexp(iT.TsTag(iC.M.o2),' ','split')    ;   % TODO: Use 'tokens' instead for sanity? Trickyish though.
    Tag3 = num2str(iC.S.TRef)                       ;
    txt = join([ "T2 decay and TE_{opt} (" Tag1 " vs " Tag2 "), for TR_{eff} = " Tag3 " ms"],"");
    title(  txt              , 'FontSize' , 16 )    ;
    ylabel( 'Rel. MR signal' , 'FontSize' , 12 )
    xlabel( 'time (ms)'      , 'FontSize' , 12 )
%     xlm = [ 0 200*ceil( TI2_max/200 ) ] ; xlim(xlm) ;   % Round up to 200 on the x axis
%     ylm = [ -1 1 ]    % ylim()          ; ylim(ylm) ;   % Rel. S is in [-1,1]
    
    if ( iC.R.TEopt == DUR )
        UI2.TEstr = "N/A"                           ;
    elseif ( iC.R.TEopt == 0 )
        iC.R.TEopt = iC.TEmin                       ;
        UI2.TEstr = string([ '(~' num2str(iC.R.TEopt) ' ms)' ]) ;
    else
        UI2.TEstr = string([      num2str(iC.R.TEopt) ' ms'  ]) ;
    end % if
    
    lxlm = [ iC.R.TEopt , iC.R.TEopt ]              ;   % Time of max signal difference
    lylm = [ St(maxind,S_atTE) , dif_max ]          ;   % Lower/upper y of line (= S_dif curve)
    line( lxlm , lylm   ,                           ... % Vertical line @ crossing
        'Color'         , 'black'               ,   ...
        'LineWidth'     , 1.5                   ,   ...
        'DisplayName'   , 'TE_{opt}'            ,   ... % (Tip: If declared after legend, the line appears in it)
        'LineStyle'     , ':'                   )   ;
    
    LegS = strings(1,NTs)                           ;   % Make a legend...
    for i = 1:NTs                                       % ...showing the tissue short names and T2 times
        LegS(i) = sprintf('%s (T2 = %i ms)', iT.TsTag(i), iT.T2(i) ) ;
    end % for i
    legend( LegS        ,                           ...
        'Location'      , 'South'               )   ;   % 'Best' may conflict with UI and TextBox
    TEo = iC.R.TEopt + 1                            ;
    iC.R.magS  =  abs( St(TEo,iC.M.o1 ) )           ;
    iC.R.magC  =  abs( dif_max        )             ;   % TODO: How does this work? It should be negative if WML<WM
    iC.R.magCW = magC( St(TEo,Tis(1)) , St(TEo,Tis(2)) );   % WML-WM contrast
    iC.R.magCG = magC( St(TEo,Tis(1)) , St(TEo,Tis(3)) );   % WML-GM contrast
    if ~isfield(iC.R,'olS2') ; iC.R.olS2 = iC.R.magS ; end  % Old signal strength   for comparison; as float
    if ~isfield(iC.R,'olC2') ; iC.R.olC2 = iC.R.magC ; end  % Old signal difference for comparison; as float
    if ~isfield(iC.R,'olTR') ; iC.R.olTR = iC.S.TR   ; end  % Old rep. time         for comparison; as float
    relS    = uint16(100*(  iC.R.magS/iC.R.olS2  )) ;   % Rel. Signal for WML
    relC    = uint16(100*(  iC.R.magC/iC.R.olC2  )) ;   % Rel. Contrast between WML and tissue
    relS_t  = uint16(100*( (iC.R.magS/iC.R.olS2)/sqrt(iC.S.TR/iC.R.olTR) )) ;   % Rel. SNR per time (ratio of S over TR^0.5)
    relC_t  = uint16(100*( (iC.R.magC/iC.R.olC2)/sqrt(iC.S.TR/iC.R.olTR) )) ;   % Rel. CNR per time (ratio of S over TR^0.5)
    fStr = ['\\itOptimal TE:\n\\rm\\bf'             ... % Tex \bf\it\rm = bold/italic/normal; escape \ for sprintf!
            ' TE_{opt}    = %3s \n'                 ... % %3.i ms
            ' S_{WML}    = %4.3f \n'                ... % %4.3f
            ' \\DeltaS_{max}  = %4.3f \n'           ... % %4.3f
            ' rS_{WML}   = %3.f%%\n'                ... % %4.3f
            ' rC(L-?)  = %3.f%%\n'                  ... % %4.3f
            ' rSNR/t   = %3.f%%\n'                  ... % %4.3f
            ' rCNR/t   = %3.f%%'                ]   ;   % %4.3f
    fVar = [ UI2.TEstr, iC.R.magS, iC.R.magC, relS, relC, relS_t, relC_t ]  ; 
    figTextBox( fStr , fVar )                       ;   % Display an info text box on the figure
    if ( ~UI2.FrzSNR )                                  % If RelZ isn't frozen (by the UI button)
        iC.R.olC2 = iC.R.magC                       ;
        iC.R.olS2 = iC.R.magS                       ;
        iC.R.olTR = iC.S.TR                         ;
    end % if FrzSNR
    
    if ( UI2.TEset ~= iC.R.TEopt )
        TEs = UI2.TEset + 1                         ;
        iC.R.magSset =  abs( St(TEs,iC.M.o1  ) )    ;
        iC.R.dirCset = magC( St(TEs,iC.M.o1  ) , St(TEs,iC.M.o2 ) ) ;
        iC.R.C_WMset = magC( St(TEs,Tis(1)   ) , St(TEs,Tis(2)  ) ) ;   % WML-WM contrast
        iC.R.C_GMset = magC( St(TEs,Tis(1)   ) , St(TEs,Tis(3)  ) ) ;   % WML-GM contrast
    end % if

end % mainTEcalc

function figTextBox( fStr, fVar )
    txt = sprintf( fStr , fVar )                    ;   % Formatted string "print"
%     lin = count( fStr , '\n' ) + 1                  ;   % Count lines in the string
    xlm = xlim  ; xt = double( xlm(2) ) * 0.780     ;   % Max value of axis; (x,y) is a plot point!
    ylm = ylim  ; yt = double( (ylm(2)+ylm(1))/2. ) ;   % Place the text box at the middle of the y axis
%     yt = yt + 0.035*lin                             ;   % Compensate for number of lines
    text( xt , yt , txt ,                           ... % Note: annotation('textbox' etc doesn't scale with text
        'EdgeColor'     , 'black'               ,   ...
        'Margin'        , 4                     ,   ...
        'LineWidth'     , 1                     ,   ...
        'LineStyle'     , '-'                   )   ;
end % fcn

%% FUNCTIONS
function FltSym = FoS( inNr )                           % Cast as double for, e.g., exp() unless it's used as a sym
    if  isa( inNr , 'sym' ) , FltSym = inNr         ;   % TODO: This exists both in the main script and T2_Subplot now.
    else,                     FltSym = double(inNr) ;
    end % if
end % fcn

function St = CalcSt( S0 , T2 , Ts )                    % Calculate generic S(t) = S0*exp(-t/T2)
    S0 = FoS(S0) ; T2 = FoS(T2) ; Ts = FoS(Ts)      ;   % Cast input as float unless it's a sym
    St = S0 * exp( -Ts./T2 )                        ;   % One signal value per given T2
end % fcn

function TEm = FindOptTE( S0_1, S0_2, T2_1, T2_2 )      % Find TE to give max contrast between two decaying tissues
    T2_1 = FoS(T2_1); T2_2 = FoS(T2_2)              ;   % Cast input as float unless it's a sym
    % TODO: How to use vpasolve to find the max?
    TEm = vpasolve( ( CalcSt(S0_1,T2_1)             ... % Solve |S1(t) - S2(t)| == max for t = TE_opt
                    - CalcSt(S0_2,T2_2) ), t)       ;   % Return TE_opt
end % fcn

function mC = magC( SL, SB )                            % Calculate contrast between "lesion" and "background" magnitude S
    mC = abs( SL ) - abs( SB )                      ;   % Use the difference between abs. signals as each voxel is mag.
end % fcn

%% PLOT

% Plot S(t) under T2 decay
% function PlotSt( T2plt, DUR )                           % DUR = Total duration for the plot
% end % fcn

%% USER INTERFACE

function createUI()                                     % Display a button to show/hide the UI controls
    eW = 0.065 ; eH = 0.065; eS = eH + 0.025        ;   % Button width/height/spc for masterUI controls
    myAx = gca                                      ;   % A handle to the current graphics axes (plot)
    uiX = myAx.OuterPosition(3) - 0.015 - eW        ;   % Use the axis' outer x position (def. [ 0 0 1 1 ]) for UI pos.
    uiY = myAx.Position(4)      + 0.010             ;   % Use the axis' inner y position to determine UI pos.
    
    uicontrol( 'Style'      , 'pushbutton'      ,   ... % For all modes:
        'String'            , 'Panel'           ,   ...
        'FontWeight'        , 'bold'            ,   ...
        'ToolTipString'     , '[h] Show/hide UI',   ...
        'Units'             , 'normalized'      ,   ...
        'Position'  , [ uiX uiY-0*eS eW eH ]    ,   ...
        'Callback'  , @masterUI_callback        )   ;   % UI to show/hide the UI control panel
    
    createUIPanel()                                 ;
    switch UI2.PanVis
        case true ; F2_Pan.Visible = 'on'           ;
        case false; F2_Pan.Visible = 'off'          ;
    end % switch UI2 PanVis
end % createUI

function createUIPanel()
    boxClr = iC.P.FigClr                            ;   % UI panel background color (same as fig background)
    mX = 10; mY =  6                                ;   % x/y margins in UI panel
    eN = 04; eW = 130; eH = 30                      ;   % Number of UI elements; default element width/height
    UI2.Hpx = eN*eH + 2*mY; uiT = UI2.Hpx - mY + 3  ;   % UI panel height in pixels; top position for elements
    UI2.Wpx = eW + 2*mX                             ;   % UI panel width in pixels
    
    F2_Pan = uipanel(         iC.P.Fig2         ,   ... % F2_Pan
        'BackgroundColor'   , boxClr            ,   ... %  [ 0.8 0.8 0.4 ]
        'BorderType'        , 'etchedin'        ,   ... % default 'etchedin'
        'BorderWidth'       , 1.0               ,   ... % default 1
        'Clipping'          , 'off'             ,   ... % default 'on'; 'off' is good for debugging
        'Units'             , 'normalized'      ,   ... % Normalized units facilitate resizing to axes
        'ResizeFcn'         , @placeUiBox       )   ;   % Resize function for the UI panel
    placeUiBox();                                       % Adjust the UI panel position to the figure window
    
    for i = 1:eN; eY = uiT - i*eH                   ;   % Generate each UI element on a new row in the panel
        switch i
        case  1
        uicontrol( 'Parent' , F2_Pan            ,   ... % ui pT2s
            'Style'         , 'pushbutton'      ,   ...
            'String'        , 'Results'         ,   ...
            'ToolTipString' ,'[p] Print results',   ...
            'Position'  , [ mX+00 eY+00 eW 24 ] ,   ...
            'Callback'  , @printVars_callback   )   ;   % UI to print info/results to the command window
%             'Units'         , 'normalized'      ,   ...
        
        case  2
        uiStr = 'TE_set:                      ms'   ;
        ttStr = '[q] Echo time'                     ;
        t2Str = 'Can be set manually'               ;
        uicontrol( 'Parent' , F2_Pan            ,   ... % ui TE[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX+00 eY+01 eW 16 ] )   ;
        UI2.TES = uicontrol( 'Parent' , F2_Pan  ,   ...
            'Style'         , 'edit'            ,   ...
            'Min' , 11      , 'Max' , 11        ,   ... % Tip: Max > Min allows multiline edit
            'ToolTipString' , t2Str             ,   ...
            'String'        , UI2.TEset         ,   ...
            'Value'         , UI2.TEset         ,   ...
            'Position'  , [ mX+40 eY+00 60 20 ] ,   ...
            'Callback'  , @te_set_callback      )   ;   % UI to set TR

        case  3
        uiStr = 'Freeze relSNR'                     ;
        ttStr = '[r] Set comparison std.'           ;
        uicontrol( 'Parent' , F2_Pan            ,   ... % uiSetS
            'Style'         , 'radiobutton'     ,   ... % 'togglebutton'/'checkbox'
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor',  boxClr          ,   ...
            'Value'         , UI2.FrzSNR        ,   ...
            'Position'  , [ mX+10 eY+00 eW 24 ] ,   ...
            'Callback'  , @setSignal_callback   )   ;   % UI to set the relative signal for comparison
%             'Min' , 0       , 'Max' , 1         ,   ... % Max/Min 0/1 are default for Matlab toggles

        case  4
        uiStr = 'Dur.:                          ms' ;
        ttStr = '[d] Plot duration'                 ;
        t2Str = 'Can be set manually'               ;
        uicontrol( 'Parent' , F2_Pan            ,   ... % ui DUR[T|S]
            'Style'         , 'text'            ,   ...
            'ToolTipString' , ttStr             ,   ...
            'String'        , uiStr             ,   ...
            'BackgroundColor', boxClr           ,   ...
            'HorizontalAlignment', 'Left'       ,   ...
            'Position'  , [ mX+00 eY+01 eW 16 ] )   ;
        UI2.DURS = uicontrol( 'Parent' , F2_Pan ,   ...
            'Style'         , 'edit'            ,   ...
            'ToolTipString' , t2Str             ,   ...
            'String'        , DUR               ,   ...
            'Value'         , DUR               ,   ...
            'Min' , 11      , 'Max' , 11        ,   ...
            'Position'  , [ mX+40 eY+00 60 20 ] ,   ...
            'Callback'  , @dur_set_callback     )   ;   % UI to set DUR
        
        end % switch i
    end % for i
end % createUIPanel

function masterUI_callback(~,~)                         % UI that shows/hides the UI control panel
    UI2.PanVis = ~UI2.PanVis                        ;
    switch UI2.PanVis
        case true ; F2_Pan.Visible = 'on'           ;
        case false; F2_Pan.Visible = 'off'          ;
    end
end % fcn

function printVars_callback(~,~)                        % UI that prints info to the command window
    fStr = ['\n'                                    ... % TODO: Isn't fprintf obeying the field width operator for strings?!
            '*******************************\n'     ...
            '***    T2 plot results:     ***\n'     ...
            '*******************************\n'     ...
            '*  S0_%3s         =  % 1.3f   *\n'     ... % TODO: Align when positive?
            '*  S0_%3s         =  % 1.3f   *\n'     ...
            '*  S0_%3s         =  % 1.3f   *\n'     ...
            '*  S0_%3s         =  % 1.3f   *\n'     ...
            '*******************************\n'     ...
            '*  TE_opt         =  %+7s','  *\n'     ...
            '*  S_WML @TE_opt  =  % 1.3f   *\n'     ...
            '*  C_WM   --"--   =  % 1.3f   *\n'     ...
            '*  C_GM   --"--   =  % 1.3f   *\n' ]   ;   % TODO: Proper right alignment of positive/negative numbers
    fVar = { iT.TsTag(1),iT.S0(1), iT.TsTag(2),iT.S0(2),    ...
             iT.TsTag(3),iT.S0(3), iT.TsTag(4),iT.S0(4),    ...
             UI2.TEstr, iC.R.magS, iC.R.magCW, iC.R.magCG } ;
    fprintf( fStr, fVar{:} );
    if ( UI2.TEset ~= iC.R.TEopt )
        fStr = [                                    ...
            '*******************************\n'     ...
            '*  TE_set         =  %4.i ms  *\n'     ...
            '*  S_WML @TE_set  =  % 1.3f   *\n'     ...
            '*  C_WM   --"--   =  % 1.3f   *\n'     ...
            '*  C_GM   --"--   =  % 1.3f   *\n'     ...
            '*  relS(set/opt)  =  %4.1i %%   *\n'   ...
            '*  relC(set/opt)  =  %4.1i %%   *\n' ] ;
        fVar = { uint16(UI2.TEset), iC.R.magSset, iC.R.C_WMset, iC.R.C_GMset,               ...
                 uint16(100*iC.R.magSset/iC.R.magS), uint16(100*iC.R.dirCset/iC.R.magC) }   ;
        fprintf( fStr, fVar{:} )                    ;
%                 '*******************************\n' ...
%                 '***     Manual echo time:   ***\n' ...
    end % if
    fprintf(    '*******************************\n' );
end % fcn

function te_set_callback(src,~)                         % UI to set TE manually
    val = str2val( src.String )                     ;
    if ( ( val >= 0 ) && ( val < iC.S.ETD ) )             % Allowed range
        UI2.TEset = val                             ;
    else
        UI2.TEset = iC.R.TEopt                      ;   % You can reset with a neg. value
    end % if
    main();
end % fcn

function setSignal_callback(~,~)                        % UI callback to toggle Freeze RelSNR
    UI2.FrzSNR = ~UI2.FrzSNR                        ;   % was src.Value;
    main();
end % fcn

function dur_set_callback(src,~)                        % UI to set DUR manually
    val = str2val( src.String )                     ;
    if ( ( val > iC.R.TEopt ) && ( val < iC.S.ETD ) )     % Allowed range
        DUR = val                                   ;
    else
        DUR = iC.S.ETD                                ;   % You can reset with a neg. value
    end % if
    main();
end % fcn

function keyPress_callback(~,evt)
    switch evt.Key
        case 'h'                                        % Show/hide UI panel
            masterUI_callback([],[])                ;   % Calling with dummy values
        case 'p'                                        % Print results
            printVars_callback([],[])               ;
        case 'q'                                        % Focus on TE setting uicontrol
            uicontrol( UI2.TES )                    ;
        case 'd'                                        % Focus on DUR setting uicontrol
            uicontrol( UI2.DURS )                   ;
        case 'r'                                        % Freeze SNR
            setSignal_callback([],[])               ;
        case 'a'                                        % Change active figure
            if isfield( iC.P, 'Fig1' )
                figure( iC.P.Fig1    )              ;
            end % if
    end % switch
end % fcn

function placeUiBox(~,~)                                % Function to resize the UI panel (automatically called)
    if ishandle( F2_Pan )                               % (The first time this fn is called, the panel doesn't exist yet)
        pxFig = getpixelposition(iC.P.Fig2)         ;   % Figure size in px, for manual normalizing
        myAx = gca                                  ;   % A handle to the current graphics axes/plot
        uiW = UI2.Wpx / pxFig(3)                    ;   % Normalized UI panel width
        uiH = UI2.Hpx / pxFig(4)                    ;   % Normalized UI panel height
        uiX = myAx.Position(3) - uiW + 0.115        ;   % The axis' inner right x pos. (but why the "fudge factor"?)
        uiY = myAx.Position(4) - uiH + 0.090        ;   % The axis' inner upper y pos. (--"--)
        F2_Pan.Position = [ uiX uiY uiW uiH ]       ;
    end % if
end % fcn

function val = str2val( str )                           % Function to process numeric input for UI controls
    val = str2double( str )                         ;
    if ( val ~= val )                                   % test for NaN from malformed input
        val = -1                                    ;
    end % if
end % fcn

end % script fn