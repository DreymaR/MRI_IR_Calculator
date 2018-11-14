function IR_T1_Magplot( TIabs )
%% (D)IR MRI Z magnetization plot
%   - This function is called from the IR_TI_Calculator
%   - Plot the longitudinal magnetization Mz of several tissue T1s under multiple inversion time points.
%   - (NB: Absolute TI times are used here, not relative to readout!)

%  NOTE:
%   - 
%
%  TODO:
%   - Add the Mz_end correction for FLAIR if iS.Vnd == "GE"? Must sort out which tissues to correct for then
%
%  DONE:
%   - 
%

%% INIT
% debugInfo  = 1;                                     % Show extra info in the command window?
global iS Ts iR                                       % Globally used data structs
global IE                                             % Globals that didn't thrive in structs or need short names

%% MAIN
    NIT = length(TIabs)-2;                          % # of inversion time points, without readout time and duration
    if ( iS.T2pOn )
        TIabs(2:NIT+1) = TIabs(2:NIT+1) + iS.T2p;   % We simulate T2 prep first in the sequence, so add it to TIs
    end % if
    TIn = TIabs(1:NIT); DUR = TIabs(NIT+2);         % Actual inversion time points & total duration for the plot
    TIR = TIabs(NIT+1); xRo = double([ TIR TIR ]);  % Readout time (=TI)
%     sZRo = 0.02*iS.ETD;                             % NOTE: sZRo is irrelevant as spin lock dictates readout from TIR?
%     TRo = xRo + [ -sZRo, (iS.ETD-sZRo) ];           % Readout time range. Linear(sZRo=0.5) or center-out (sZRo small)?
    TRo = xRo + [ 0, iS.ETD ];                      % Readout time range. Due to spin lock, this should start at TIR.
    dT = 1.0; TPv = (0:dT:DUR-dT);                  % Sampling interval dT (ms, Int), TPv is the time point vector
    
    T1Plt = Ts.Plt(1,:); T2Plt = Ts.Plt(2,:);
    Ntis = length( Ts.Plt(1,:) );                   % # of relaxation times for which to calculate/plot
%     if ( ~iS.T2pOn )                                % WIP: Plotting extra tissues is hard if T2 in T2p is unknown!
%         if ~ismember( iT.T1z, T1Plt )               % Add T1z if not in the preselected T1 ensemble
%             Ntis  = Ntis + 1;
%             T1Plt = double( [ T1Plt, iT.T1z ] );
%             T2Plt = double( [ T2Plt, iT.T1z ] );    % WIP: Estimate T2 (or look it up)? No.
%             Ts.Leg = [ Ts.Leg, sprintf('T1#%s', num2str(Ntis)) ];
%         end % if ~ismember iT.T1z
%     end % if iS.T2pOn
    Mz = ones( Ntis, DUR );                         % MagZ for all tissues and time points (was zeros(NT1s+1...)
    
    % WIP: Starting Z magnetization calculation.
%     Mz0 = ones(1,length(T1Plt));                    % Test/debug: Start with relaxed signal (TR >> T1)
%     Mz0 = -IE*Mz0;                                  % --"--: Start with inversion
%     Mz0 = 0*Mz0;                                    % --"--: Start with saturation
%     Mz0 = -CalcS0( [ iS.TRef, 0 ], [ T1Plt ] )      % --"--: CalcS0 uses TI differences.
%     Mz0 = CalcMagZ( [ 0 ], iS.TRef, T1Plt )         % --"--: [1 - exp(-iS.TRef./T1Plt)]
    % TODO: Use derivation instead of iteration to arrive at SS? Maybe not so critical. Also, nice to simulate it.
    RTI = zeros( 1, iS.ETL);
    for i = 1:iS.ETL
        RTI(i) = TRo(1) + round(i*iS.ESP);          % Array of readout refocus times. TODO: Avoid rounding?
    end % for
    
    its=3; iR.MzIt = ones(Ntis,its);                % Iterate to get to a steady state. Only need to run twice?
    for i = 1:its                                   % Iteration loop (Note: Needs only two TR if Mz(TRo)=0)
        Mz(:,1) = Mz(:,DUR);                            % The starting magnetization (before first inv.) is MagZ(TR)
        for j = 2:DUR                                   % Step through all time points
            tpt = double( (j-2)*dT );                   % Was (j-1)*dT, but I want the time point 0 to be defined
            Mz1 = Mz(:,j-1);                            % MagZ at the previous time point
            if ismember( tpt, TIn )                     % Inversion time (was any(Ts == TIn(:)) earlier)
                IPBeg = tpt; IPEnd = tpt + iS.IPD;      % Pulse time (For now, I use TI as the start of the pulse...)
%                 if ( IPBeg < 0 )
%                     IPBeg2 = IPBeg + iS.TR;             %   WIP: If some pulse time is < 0, move it to the end of TR?!?
%                     IPEnd2 = IPEnd + iS.TR;             %   Or, should I just use positive times (GE style)?!?
%                 else % ( IPEnd > iS.TR )
%                     IPBeg2 = IPBeg - iS.TR;
%                     IPEnd2 = IPEnd - iS.TR;
%                 end % if
                MzAtIP = Mz1;
            end % if tpt                                % (Used a new if here because that seems more robust)
            if ( tpt >= IPBeg ) && ( tpt < IPEnd )
                rPt = ( tpt*2 - (IPBeg+IPEnd) )/( IPEnd-IPBeg );   % Scan the "relative pulse time" [-1,1]
                rMz = -atan( 6.0*rPt )/atan( 6.0 );     % A generic scaled sigmoid. WIP: Not sure if it's right for Mz!
                Mz(:,j) = MzAtIP.*(rMz*(1+IE)/2 + (1-IE)/2);  % Account for Inv.Eff.
            elseif ( tpt == IPEnd ) % || ( tpt == IPEnd2 )
                Mz(:,j) = -IE*MzAtIP;                   % Invert Mz, accounting for Inv.Eff. (not trusting the sigmoid)
            elseif ( tpt == TRo(1) )                    % Readout time
                iR.MzTI = Mz1;                          % MagZ at readout
                iR.S0   = sin(iS.FAx*pi/180)*Mz1;       % Signal strength at readout
                Mz(:,j) = cos(iS.FAx*pi/180)*Mz1;
                MzRo    = iR.MzTI - Mz(:,j);
%             elseif ( Ts > TRo(1) ) && ( Ts < TRo(2) )   % Readout with CPMG train, so signal doesn't evolve (so much)
            elseif ismember( tpt, RTI )                 % CPMG refocus time (rounded to the nearest ms...)
                Mz(:,j) = -iS.IET*Mz1;
%                 for k = 1:NT1s                        % TODO: Evolve slowly during readout with adjusted T1s? Or?!?
%                     Mz(k,j) = 1 - (1 - Mz(k,j-1))*exp(-dT*T1rhF/T1Plt(k));
%                 end % for k
            elseif ( tpt == TRo(2) )                    % Readout end
                Mz(:,j) = Mz1 + iS.Rew*MzRo;            % A drive pulse reclaims some excitation signal as Mz? Or can it?!?
%             elseif ( Ts > iS.TR - iS.T2p )              % T2 preparation (effectively decays Mz with T2)
            elseif ( tpt > iS.IPD ) && ( tpt < iS.T2p + iS.IPD ) && ( iS.T2pOn )     % T2 preparation (effectively decays Mz with T2)
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
        iR.MzIt(:,i) = Mz(:,DUR);                   % The end Z magnetization for all T1s, for this iteration
    end % for i
    iR.MagZ = Mz(:,DUR)';                           % Return the end Z magnetization for all T1s as a vector
    
%     Mz = [ Mz; zeros(1,DUR) ];                    % Add a zero line to the plot
    plot( TPv, Mz, '-',                         ... % Plot the (steady state) Z magnetization
        'LineWidth', 1.5                        )
    
    xlabel( 't (ms)'      , 'FontSize', 12 )
    ylabel( 'Rel. Z magn.', 'FontSize', 12 )
    xlm = [  0  DUR ]; xlim(xlm);
    ylm = [ -1   1  ]; ylim(ylm);
%     set( gca, 'XTick', [] );                        % Remove x axis ticks
%     set( gca, 'YTick', [] );                        % Remove y axis ticks

    line( xlm, [ 0 0 ],                             ... % Zero line
        'Color', 'black',                           ...
        'LineWidth', 0.7,                           ...
        'LineStyle', ':'                            );
    lineY = [ -0.975 -0.975 ]; lineW  = 10;             % Y position and width for horz. marker lines
    line( TRo, lineY,                               ... % Horizontal line representing readout
        'DisplayName', 'Echo train',                ...
        'Color', [ 0.3 0.6 0.3 ],                   ...
        'LineWidth', lineW,                         ...
        'LineStyle', '-'                            );
    
    if ( iS.T2pOn )
        line( [ 0 iS.T2p ]+iS.IPD, lineY,           ... % Horz. line representing T2prep
            'DisplayName', 'T2 prep.',              ...
            'Color', [ 0.9 0.6 0.3 ],               ...
            'LineWidth', lineW,                     ...
            'LineStyle', '-'                        );
    end % if T2pOn

    for i = 1:NIT
    line( [ TIn(i) TIn(i)+iS.IPD+10 ], lineY,       ... % Horz. line(s) representing inv. (a little extra for display)
        'DisplayName', 'Inversion',                 ...
        'Color', [ 0.5 0.3 0.7 ],                   ...
        'LineWidth', lineW,                         ...
        'LineStyle', '-'                            );
%         'LineJoin', 'chamfer',                      ... % Note: This doesn't round corners unless they're joined.
    end % for i
    
%     line( xRo, ylm,                             ... % Vertical line at TEeff (=TI)
%         'Color', 'black',                       ... % (Tip: If declared after legend, the line appears in it)
%         'DisplayName', 'TI_{1}',                ... % 'TE_{eff}'
%         'LineStyle', ':'                        );
    
    if ( iR.RCS ~= 0 )
%         fprintf( "NIT:\n" );      % DEBUG
%         disp( NIT );              % DEBUG
        y = (-1)^(NIT+1)*iR.RCS;                        % Sensible imperfect nulling depends on number of inversions
        line( [ 0 TRo(1) ], [ y y ],                ... % Horz. line for desired residual CSF signal (imperfect nulling)
            'Color', [ 0.929 0.694 0.125 ],         ... % (Tip: lines(#) gives the color scheme for # plot lines)
        'LineWidth', 1.0,                           ...
        'LineStyle', '-.'                            );
    end % if iR.RCS
    
    LegS = strings(1,Ntis);                             % Figure legends
    for i = 1:Ntis
        LegS(i) = sprintf('%4s (T1: %4i ms)', Ts.Leg(i), Ts.Plt(1,i) );   %legend('WM', 'GM', 'CSF'...
    end % for i
    legend( LegS, 'Location', 'NorthWest');

end % script fn