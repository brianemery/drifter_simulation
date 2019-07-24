function c1 = arrowplot(arg1,arg2,arg3,arg4,arg5,arg6)
%ARROWPLOT  draws vectors at every (x,y) pair in vectors X and Y.	
%	ARROWPLOT(X,Y,DX,DY) draws little arrows at every (x,y) pair in
%	the position vectors X and Y. The (dx,dy) pairs in vectors DX and DY
%	determine the direction and magnitude of the arrows.  This program
%	expects at least 4 arguments (X,Y,DX,DY) and at most 6.
%
%	ARROWPLOT(X,Y,DX,DY,S) applies the scalar S as a scale
%	factor to the lengths of the arrow. For example, S = 2 doubles
%	their relative length and S = 0.5 halves them.
%
%	A final trailing string argument specifies linetype and color using
%	any legal line specification as described under the plot command.
%
%   C1 = ARROWPLOT(X,Y,DX,DY,S,'R')
%	The output argument C1 is the handle to the little arrows so the
%   user can change the color of the arrows.
%
%   NOTE:  The arrows will be distorted (i.e. the direction may not
%          be true) if the length of a unit on the x and y axes are
%          not the same.  This can be taken care of 2 ways:
%
%          1)  Set the x and y units to be the same length with the
%              axis('equal') command.
%          2)  If 1) isn't possible, scale the vector to match the
%              scaling of the axes.  See help on the 'aspectratio'
%              of the axes command in the Matlab Reference Manual
%              to see how the axis scaling factors can be accessed.
	
%	Modified by Mike Cook from the Matlab ver. 4.0 quiver.m - NOV 93
%   Modified by Mike Cook to pass vectors handle back to calling
%	program  -  FEB 94.
%   Modified by Mike Cook to fix when axis is drawn when minx=maxx or 
%               miny=maxy.  OCT 96.
%
%	Charles R. Denham, MathWorks 3-20-89
%	Modified 12-19-91, LS.
%	Copyright (c) 1984-92 by The MathWorks, Inc.

x = [0 1 .8 1 .8].';
y = [0 0 .08 0 -.08].';
arrow = x + y.*sqrt(-1);

eval(['last = arg' int2str(nargin) ';']);
if isstr(last)
	eval(['style = arg' int2str(nargin), ';']);
	narg = nargin - 1;
else
	style = '''Color'',''k'''; % the default
	narg = nargin;
	eval(['last = arg' int2str(narg-1) ';']);
	if isstr(last)
		error('Only the final argument can be a string.');
	end
end
if narg == 0
	error('First 4 arguments must be numeric.')
end
eval(['lastdim = max(size(arg' int2str(narg) '));']);
if lastdim == 1
	if isstr(eval(['arg' int2str(narg)]))
		error('Scalar scale argument expected.')
	end
	eval(['scale = arg' int2str(narg) ';']);
	narg = narg - 1;
else
	scale = 1;
end
if narg == 0
	error('First 4 arguments must be numeric.')
end

if narg < 4
        error('MUST supply the following 4 arguments: x,y,dx, and dy.')
else
	if isstr(arg1) | isstr(arg2) | isstr(arg3) | isstr(arg4)
		error('First 4 arguments must be numeric.')
	end
	xx = arg1;
	yy = arg2;
	px = arg3;
	py = arg4;
end

% figure out delx and dely so spacing is accounted for in z
%%mc delx = xx(1,2)-xx(1,1);
%%mc dely = yy(2,1)-yy(1,1);

grid = xx + yy.*sqrt(-1); grid = grid(:);
px = px(:); py = py(:);
maxlen = max(sqrt(px.^2 + py.^2));
%%mc maxlen = max(sqrt((px/delx).^2+(py/dely).^2));

z = (px + py.*sqrt(-1)).';

%   This scale factor is scaled in part by the maximum length 
%   of the vectors, so legend vectors from vector field to 
%   the next vector field would fluctuate.  Change to make the
%   scale factor absolute, i.e. don't scale by a variable other
%   than the scale factor, use 1 as default or arg5 if supplied.
%%scale = scale*0.90 ./ maxlen;

a = scale * arrow * z + ones(5,1) * grid.';

% append nan's so we get one handle
a = [a; nan*ones(1,size(a,2))];
a = a(:);

cax = newplot;
%c1=plot(real(a), imag(a), style); 	% Plots the vectors

if length(style)<5, style=['''Color'',''' style ''''];, end
eval(['c1=plot(real(a), imag(a), ' style ');']); 
next = lower(get(cax,'NextPlot'));

if ~ishold
	minx = min(min(xx));
	miny = min(min(yy));
	maxx = max(max(xx));
   maxy = max(max(yy)); 
   
	if minx == maxx
	    if miny ~= maxy
	       set(gca,'ylim',[miny, maxy]);
	    end
	elseif miny == maxy
	    if minx ~= maxx
	       set(gca,'xlim',[minx, maxx]);
	    end
	else
	    set(gca,'xlim',[minx, maxx],'ylim',[miny, maxy]);
	end


	view(0,90);
end
