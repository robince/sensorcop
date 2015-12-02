function splitcolorplot(x,y,val,fmt1,fmt2)
% SPLITCOLORPLOT Plots graph split into two colors above and below a given threshold value
%   SPLITCOLORPLOT(X,Y,VAL) plots vector Y against vector X, using one
%   color when Y > VAL and another when Y < VAL.  If VAL is not specified,
%   a dividing value of 0 is used.
%   As with PLOT, SPLITCOLORPLOT(Y,VAL) plots Y against index numbers,
%   instead of X values.
%   SPLITCOLORPLOT(X,Y,VAL,FMT1,FMT2) uses standard PLOT formatting strings
%   to define the color, linestyle, and marker of the two lines.  FMT1
%   defines the formatting for the line when Y > VAL.
%
%   Acceptable calling syntaxes:
%     SPLITCOLORPLOT(Y)                     [index used for X & VAL = 0]
%     SPLITCOLORPLOT(X,Y)                   [VAL = 0]
%     SPLITCOLORPLOT(Y,VAL)                 [index used for X]
%     SPLITCOLORPLOT(X,Y,VAL)
%     SPLITCOLORPLOT(X,Y,FMT1,FMT2)         [VAL = 0]
%     SPLITCOLORPLOT(Y,VAL,FMT1,FMT2)       [index used for X]
%     SPLITCOLORPLOT(X,Y,VAL,FMT1,FMT2)
%
%   Example:
%   >> x = 0:0.1:10;
%   >> y = sin(4*x);
%   >> splitcolorplot(x,y,0.4)
%   >> splitcolorplot(y,-0.2,'rx--','k*:')

%  Parse inputs (nested function used to set x, y, val, and all necessary
%  format options)
checkinputs(nargin)

%  Find zero crossings ("zero" = specified level)
if isrow(y)
    idx = [0,find(sign(y(2:end)-val)~=sign(y(1:end-1)-val)),lngth];
else
    idx = [0;find(sign(y(2:end)-val)~=sign(y(1:end-1)-val));lngth];
end
nidx = length(idx);

%  Make an invisible plot -- this means that splitcolorplot has the same
%  overwriting behavior as plot
plot(x,y,'visible','off')

%  Loop over each segment
for k=2:nidx
    %  Get first and last indices
    k1 = idx(k-1)+1;
    k2 = idx(k);
    %  Choose format, depending on whether y > level or not
    fmt = (y(k1)<=val) + 1;
    %  Plot line segment
    line(x(k1:k2),y(k1:k2),'color',clr{fmt},'marker',mrkr{fmt},'linestyle',lntp{fmt})
    %  Plot linear interpolation to beginning of next line segment
    if k<nidx
        %  Start point is last point of previous line segment
        %  End point is beginning of next segment
        x1 = x(k2);
        x2 = x(k2+1);
        %  Get associated y values
        y1 = y(k2);
        y2 = y(k2+1);
        %  Slope
        m = (y2-y1)/(x2-x1);
        %  Find point on line where y = val
        x0 = x1 + (val-y1)/m;
        %  Plot from end of first segment to y = val
        line([x1,x0],[y1,val],'color',clr{fmt},'linestyle',lntp{fmt})
        %  Change formats and plot from y = val to start of 2nd segment
        fmt = mod(fmt,2) + 1;
        line([x0,x2],[val,y2],'color',clr{fmt},'linestyle',lntp{fmt})
    end
end

    %  Parse all the inputs
    function checkinputs(n)
        %  As long as we have some inputs, get the length of the first
        %  input (sort out other possible errors later)
        if n
            lngth = length(x);
        end
        %  Assign values and defaults as best we can (and trust the later
        %  error checking to sort out any problems)
        %  How many Romans?!
        switch n
            %  i
            case 1
                % Single input is vector of y values
                y = x;
                % Everything else is default
                x = 1:lngth;
                val = 0;
                fmt1 = 'o-';
                fmt2 = 'o-';
            %  ite!
            case 2
                % Two inputs could be x & y or y & level
                % Two equal-sized vectors => x & y
                if isequal(size(x),size(y))
                    val = 0;
                % Scalar y => y is level (and x is actually y)
                elseif numel(y)==1
                    val = y;
                    y = x;
                    x = 1:lngth;
                else
                    error('x & y must be same size, and level must be a scalar')
                end
                fmt1 = 'o-';
                fmt2 = 'o-';
            %  if you don't understand these comments, you seriously need
            %  to watch "Life of Brian"
            case 3
                % x, y, and level are all set, so set default formats
                fmt1 = 'o-';
                fmt2 = 'o-';
            case 4
                % Four inputs => last two are formats
                fmt2 = fmt1;
                fmt1 = val;
                % First two inputs could be x & y or y & level
                % Two equal-sized vectors => x & y
                if isequal(size(x),size(y))
                    val = 0;
                % Scalar y => y is level (and x is actually y)
                elseif numel(y)==1
                    val = y;
                    y = x;
                    x = 1:lngth;
                else
                    error('x & y must be same size, and level must be a scalar')
                end
            case 5
                % Nothing to do here
            otherwise
                % No.
                error('Wrong number of inputs')
        end
        %  Now error-check the hell out of the x, y, val, fmt1, & fmt2 that
        %  we ended up with.
        if ~isequal(size(x),size(y))
            error('x & y must be same size')
        elseif ~isvector(x)
            error('x & y must be vectors')
        elseif ~isscalar(val)
            error('Level must be a scalar')
        elseif ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(val)
            error('x, y, and level must be numeric')
        elseif ~ischar(fmt1) || ~ischar(fmt2)
            error('Formats must be strings')
        end
        %  Parse and deconstruct the format strings
        %  Look for color specifiers first
        [clr{1},fmt1] = getcolor(fmt1,'b');
        [clr{2},fmt2] = getcolor(fmt2,[0 0.5 0]);
        %  And now for markers
        [mrkr{1},fmt1] = getmarker(fmt1);
        [mrkr{2},fmt2] = getmarker(fmt2);
        %  And what's left is a linestyle
        fmt1 = getlnstyle(fmt1,mrkr{1});
        fmt2 = getlnstyle(fmt2,mrkr{2});
        lntp = {fmt1,fmt2};
    end

    function [c,fmt] = getcolor(fmt,def)
        %  Possible values
        clrs = '[bgrcymkw]';
        %  See if fmt1 includes any of the possible values
        idx = regexp(fmt,clrs);
        if isempty(idx)
            % No?  Set default
            c = def;
        else
            % Yes?  Cool, that's the color.  Remove all color specifiers
            % from the string.
            if numel(idx)~=1
                warning('SPLITCOLORPLOT:multipleColors',...
                    'More than one color detected')
            end
            c = fmt(idx(1));
            fmt(idx) = [];
        end
    end

    function [m,fmt] = getmarker(fmt)
        marks = '[.ox+*sdv^<>ph]';
        idx = regexp(fmt,marks);
        %  Special case: if there's a dot (.), we need to make sure it
        %  isn't part of the linestyle specifier -.
        if ~isempty(regexp(fmt,'.','once'))
            %  Go through all possibilities (.s) in reverse order --
            %  reverse because we will remove any that shouldn't be there,
            %  and that will mess up the indexing
            for j=length(idx):-1:1
                %  If there's a - in front of the ., remove it from the
                %  list of possibilities (it's a line specifier, not a
                %  marker)
                if idx(j)>1 && strcmp(fmt(idx(j)-1:idx(j)),'-.')
                    idx(j) = [];
                end
            end
        end
        % OK, now go through and work out marker
        if isempty(idx)
            % Nothing given?  Set default
            m = 'none';
        else
            % Otherwise, that's the marker.  Remove all marker specifiers
            % from the string.
            if numel(idx)~=1
                warning('SPLITCOLORPLOT:multipleMarkers',...
                    'More than one marker type detected')
            end
            m = fmt(idx(1));
            fmt(idx) = [];
        end
    end

    function fmt = getlnstyle(fmt,m)
        %  If there's nothing left, no linestyle was specified.  However,
        %  the depends on whether a marker was given
        if isempty(fmt)
            % No marker (and no linestyle) => default solid line
            if strcmp(m,'none')
                fmt = '-';
            % Marker given but no linestyle => no line (just marker)
            else
                fmt = 'none';
            end
        %  Something was specified.  Was it valid?
        elseif ~ismember(fmt,{'-','--',':','-.'})
            error('Unknown format specification')
        end
    end

end