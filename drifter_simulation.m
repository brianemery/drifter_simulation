function DRFT = drifter_simulation(TUV, CFG, TURB)
% DRIFTER SIMULATION.M - simulate drifter trajectories using CODAR data
% DRFT = drifter_simulation(TUV, CFG, TURB)
%
% Simulate lagrangian drifters using codar HF radar data. Optionally, add a 
% markovian turbulent component. 
%
% OUTPUTS
% Matricies of drifter lat lons, with each row a drifter, and each column 
% a time. The matrix is padded with NaN's, so drifters can start and end 
% at times other than time=1 and time=end times ... all bundled into a
% structure (eg SIM)
% 
% REQUIRED INPUTS
%
% TUV standard total data structure. The TUV.TimeStamp variable defines the
%       length and timestep of the simulation run
%
% CFG (Configurations):
%     CFG.deploy_locations  = deployment location [lon lat]*
%     CFG.deploy_times = deployment times (generally same # rows as locn)
%
% * Note that the best input is to specifty a lonlat and a deploy time for
% every desired simulated drifter explicitly.
%
% OPTIONAL CFG fields:
%     CFG.ang_diff = [20 160];
%     CFG.maxDist = 50; ?? does?
%     CFG.check_plots = 1; 
%     CFG.drive =  '/Data/totals/';
%     % CFG.eof_filled = 1;                  % Deprecated ...
%     % CFG.eof_data_path = mainPath;
%     CFG.randn_state = sum(100*clock);      % allows state to be reused
%
%
% TURB structure, containing settings for turbulence model
% eg:
% TURB.n_clusters = 100;  % This is the # of simulation clusters OBSOLETE?
% TURB.add = 1;           % Switch 'noise' on or off
% TURB.stdev_u = 2.83;    % cm/sec - Standard deviation of turbulent velocities
% TURB.stdev_v = 2.83;    % cm/sec - these are 2.3238  3.4496 in paper (near/off shore)
% TURB.Tu = 179;          % min. 'Lagrangian Integral Time Scale', aka decorrelation time
% TURB.Tv = 183;          % min.
% TURB.dt = 60;           % min. Model delta time in minutes (must be less than Tu, Tv)
% TURB.sd = false; ** need info about this ....
% 
% EXAMPLE: (OLD)
% locn=[-120.723333 34.61033333];
% st=datenum(2002,5,31,15,0,0):4/24: datenum(2002,6,30,23,0,0);
% [Lon,Lat,T,iStart,iEnd]=drifter_simulation(locn,st);

% Copyright (C) 2009-2011 Brian Emery
% Version 1.0, Using Krisada's calcposdrifter.m and supporting RungaKutta code
% Version 8.0 14jul2009 Major Major revisions, reorganizing, turbulence... etc
% Version 8.1 15nov2009 Added variance as function of cluster separation
% Version 8.2 4Jan2010 newer/better calcposdrifter using interp3
% Version 9.0 14Jan2011 generalized time setp, check configurations
%                       moved the total data clean up and loading to 
%                       calling function, made TUV an input

% TODO
% Generalize timestep. Notes:
% The km grid is a sticking point. It needs to be regular for the 3-d gridding
% to work 
%    
% Update example
% 
% TEST CASE IDEA
% use real drifter data to create radar 'totals' - then use this mfile to
% recreate the drifter data
%
% Parallelize loops?
%
% Cull unneeded code that is sticking around since this was way less
% generalized (eof_fill field for example, a lot of plotting stuff)

% --------------------------------------------------------- 
%  DEFAULT SETTINGS / INPUT CHECKS
%---------------------------------------------------------- 

% If CFG not defined, set as empty, check fields later
if nargin < 2, disp('CFG not defined'), keyboard, end


% If TURB not defined, set switch
if nargin < 3
    TURB.add = 0;
else
    % otherwise set the delta time in minutes
    TURB.dt = mode(diff(TUV.TimeStamp))*1440;
end


% Check CFG fields, set to default if missing
CFG = cfg_check(CFG);


% TOTAL GRID IN KM
% Total grids loaded from CODAR OS total files will have the grid already
% present. The grid needs to be uniform in order for the grid_codar_data
% subfunction to work right (I think )

% Try to get center of grid as defined in TUV data
try
    loc.central = TUV.LonLat(TUV.OtherSpatialVars.X ==0 & TUV.OtherSpatialVars.Y ==0,:);
    HFx = TUV.OtherSpatialVars.X;
    HFy = TUV.OtherSpatialVars.Y;
    

catch
    % valid for SBC only
    loc = codar_sites({'central'});

    % convert gridd in (lon,lat)  to griddXY in km
    % [HFx,HFy] = lonlat2km_dev(central_loc(1), central_loc(2),TUV.LonLat(:,1), TUV.LonLat(:,2));
    [HFx,HFy] = lonlat2km(loc.central(1), loc.central(2),TUV.LonLat(:,1), TUV.LonLat(:,2));

end

% --------------------------------------------------------- 
%  INITIALIZE 
%---------------------------------------------------------- 

% Initialize randn to a different state
randn('state',CFG.randn_state);                                               


% CREATE MATRICIES FOR SIMULATED DRIFTERS
% Create X,Y, T matricies for the sim. drifters, a row for each drifter and
% a column for each time. This will be all NaN except where a drifter is
% being deployed. The code computes the postitions in time and replaces the
% NaN's with X,Y values
[drftX,drftY]= init_simulation_matrix(CFG.deploy_locations,CFG.deploy_times,loc.central,TUV.TimeStamp);  


% set initial turbulent velocities to zero (for markov model)
[u0,v0] = deal(zeros(size(drftX,1),1));


% OPTIONALLY, INITIATE PLOTTING
% Option to plot based on CFG.check_plots setting 
[cstr,hh,hd] = init_plot(size(drftX,1),TURB,CFG,loc.central,drftX,drftY,'c.');


% PUT CODAR DATA INTO 3-D MATRICIES
% in km space, creating the grid from HFx, HFy. Round grid to 1m. I assume
% this is ok to do, errors in the HF data are likely much larger. 
[X,Y,U,V] = grid_codar_data(round(HFx.*1e3)./1e3,round(HFy.*1e3)./1e3,TUV.U,TUV.V);
    


% ---------------------------------------------------------
%  RUN SIMULATION                                             
%----------------------------------------------------------
% Loop over HFR data and run the simulation calculation
%
% This works one drifter, one time step at a time
% Everything below here is done in km space ...


% Verbosity
disp(['Running drifter_simulation.m from ' ...
             datestr(TUV.TimeStamp(1)) ' to ' datestr(TUV.TimeStamp(end))])
disp(['Begin: ' datestr(now)])
         

% LOOP OVER TIMES
% Technically, this is looping over the codar data.                    
for NT = 1:size(U,3)-1   
    

    % Get indecies (rows) of valid drifters for this time step
    validRows = find(~isnan(drftX(:,NT)+drftY(:,NT)));
    
    disp(['... running for ' datestr(TUV.TimeStamp(NT)) ', ' num2str(length(validRows)) ' drifters'])

    
    % need to exit loop if validRows is empty!
    if isempty(validRows), break, end

    % Optional plotting 
    hh = codar_plot(hh,TUV,NT,validRows,CFG);

    
    % POSITION FROM 'MEAN' FLOW COMPONENT
    % Calculate drifter trajectories using fourth-order runge-kutta method.
    % Note that 'for', below, uses each column of the input variable
    % Loop over each drifter
    for drftIdx = validRows'; 

        [nextXpos,nextYpos] = calc_drifter_position(drftX(drftIdx,NT), ...
                                         drftY(drftIdx,NT),X,Y,U,V,NT,TUV.TimeStamp);

        % place the computed postitions in the X,Y matricies
        drftX(drftIdx,NT+1) = nextXpos;
        drftY(drftIdx,NT+1) = nextYpos; 
    end
    
    

    % ADD POSITION FROM TURBULENT COMPONENT
    % Use Joint Markovian Model (x and u) to model turbulence
    % Do this outside the loop over each drifter, since this does
    % not depend on the Codar HF flow field. Works on one colum of the
    % drifter position matricies at a time.
    if ~isempty(validRows) && TURB.add 
        
        if  lfnn(drftX(validRows,NT+1)+drftY(validRows,NT+1)) > 0
            
        [drftX(:,NT+1),drftY(:,NT+1),u0,v0] = add_joint_markovian(drftX(:,NT+1),drftY(:,NT+1),u0,v0,TURB,validRows);
        
        else
            keyboard
        end
        
    end


    % Optional plot of EACH drifter
    hd = drifter_plot(loc.central,drftX(validRows,:),drftY(validRows,:),cstr(validRows,:),hd,CFG);

    
end 


% Optional final check plot
% last_check_plot(loc.central,drftX,drftY,size(drftX,1),CFG)



% --------------------------------------------------------- 
%  FINALIZE OUTPUTS 
%---------------------------------------------------------- 
% Init output structure
DRFT = DRFTstruct(size(drftX));

% Create the output lon lats 
[DRFT.Lon,DRFT.Lat] = km2lonlat(loc.central(1),loc.central(2),drftX,drftY);

% finds the start and end matrix addresses of each drifter 
[DRFT.iStart,DRFT.iEnd] = find_first_last(DRFT.Lon+DRFT.Lat);

% other stuff
DRFT.FixTimes = TUV.TimeStamp;
DRFT.TimeStamp = TUV.TimeStamp;
DRFT.X = drftX;
DRFT.Y = drftY;
DRFT.GridOrigin = loc.central;
DRFT.ProcessingSteps = {mfilename};
DRFT.Type = 'Simulated Drifter Data';

% Store Configurations
DRFT.OtherMetadata.CFG = CFG;
DRFT.OtherMetadata.TURB = TURB;

% Final check plot
% spaghetti_plot(outLon,outLat,runTimes,size(drftX,1),TURB,CFG)


disp(['Done at ' datestr(now)])

 
end



% -----------------------------------------------------------------------
% LOCAL FUNCTIONS
function [x,y] = init_simulation_matrix(lonlat,deployTime,central,runTimes)
% INIT_SIMULATION_MATRIX
% [x,y]= init_simulation_matrix(lonlat,deployTime,central,runTimes)
% Create X,Y matricies for the sim. drifters, a row for each drifter and
% a column for each time. This will be all NaN except where a drifter is
% being deployed. The code computes the postitions in time and replaces the
% NaN's with X,Y values. Output columns correspond to times in runTimes.
%
% INPUTS:
% size(lonlat,1)=length(deployTime) or length(deployTime)==1 ...

% Brian Emery

% to do:
%  one location, every tstep
%  multi location, every tstep
%  one location, only one tstep

% SORT INTPUTS
% Note that the default input is to specifty a lonlat and a deploy time for
% every drifter explicitly. This option covers most every deployment
% scenario
if size(lonlat,1)~=length(deployTime) && length(deployTime)==1
    % Case of multiple deploy locations, one deploy time
    deployTime(1:size(lonlat,1))=deployTime;

elseif size(lonlat,1)~=length(deployTime)
    % we dont know what to do in this case
    disp('Expecting size(lonlat,1)=length(deployTime) or length(deployTime)==1')
    keyboard

end


% get unique deployment times to use as column indecies
st=unique((round(deployTime*1440))/1440);

% Using max here will give the outputs a column for each run time, but st
% is still used to get the column index.
lon=NaN(length(lonlat),length(runTimes));
lat=NaN(length(lonlat),length(runTimes));

for i = 1:size(lonlat,1)
    lon( i, (round(deployTime(i)*1440)/1440) == (round(runTimes*1440)/1440) ) =  lonlat(i,1);
    lat( i, (round(deployTime(i)*1440)/1440) == (round(runTimes*1440)/1440) ) =  lonlat(i,2);
end


% get rid of all NaN rows:
ll=lon+lat;
lon=lon(logical(sum(~isnan(ll),2)),:);
lat=lat(logical(sum(~isnan(ll),2)),:); clear ll

% convert to km grid
[x,y] = lonlat2km(central(1), central(2),lon,lat);

if size(lonlat,1)~=size(x,1)
    disp(['Only using ' num2str(size(x,1)) ' of ' num2str(size(lonlat,1)) ' input locations'])
    %keyboard
end

% % ** add this here?? **
% 
% keyboard 
% 
% % VERIFY INPUT DATA
% % The time column indecies of codar data and drifter data should correspond
% if size(U,3) ~= size(drftX,2) || size(U,3) ~= length(runTimes)
%     disp('drifter/codar data size discrepancy - assumptions violated'), keyboard  %<---- ** TIMESTEP? **
% end




end
% -----------------------------------------------------------------------
function [nextXpos,nextYpos,ut,vt] = add_joint_markovian(nextXpos,nextYpos,u0,v0,TURB,validRows)
% ADD JOINT MARKOVIAN
% [nextXpos,nextYpos,ut,vt] = add_joint_markovian(nextXpos,nextYpos,u0,v0,TURB)
% 
% To the existing new position (just one column), add a turbulent 'memory' 
% term and a random impulse term. X and Y positions are in km. 
% 
% TURB structure must contain:
%
% TURB.dt
% TURB.Tu
% TURB.Tv
% TURB.SIGMAu <--- these used when sigma is fixed
% TURB.SIGMAv
%
% This model currently uses a function for variable sigma derived from
% drifter data as a function of separation distance.
%
% REFERENCE:
% Stochastic Modelling in Physical Oceanography, pg 113-120 (equations 6
% and 7 in chapter 4 of the section).

% Copyright (C) 2009-2010 Brian Emery
% Version 1.0 10Aug2009
% standard deviation based on separation dist added Nov 2009

% Check inputs
field_check(TURB,{'dt','Tu','Tv','sd'})

% Turbulent u is added to only positions which were valid at the previous
% time step
if nargin<6
    validRows = 1:size(nextXpos,1);
end

% make output u's same size as u0's but ZERO
ut=u0.*0;
vt=v0.*0;

% simplfy equation code by defining these:
dt=TURB.dt;        % delta time (in minutes)
Tu=TURB.Tu;        % turbulent decorrelation time (Set in Minutes)
Tv=TURB.Tv;

% --------------------------------
% DETERMINE SIGMA
% --------------------------------
% Also know as the standard deviation. See compute_variance_vs_separation.m
% on how these numbers are computed. 

if TURB.sd

    % SIGMA FUNCTION OF SEPERATION DISTANCE

    % Compute the separation distance of all the valid drifters
    sd = compute_sd(nextXpos(validRows,1),nextYpos(validRows,1));

    % compute the variance, and standard deviation
    % sd in km. 'min' part is to set the Max at 2km separation distance
    % The ./sqrt(2) comes from the fact that the equation gives us total
    % variance and we want the stdev of each component.
    % Convert [cm/s] to [km/min]
    std_u = ( min([TURB.max sqrt(TURB.A* (mean(sd).^TURB.x))]) ./ sqrt(2) ).*(60./1e5); 
    std_v = std_u;

else
    % FIXED SIGMA
    
    % Convert [cm/s] to [km/min]
    std_u=TURB.stdev_u.*(60./1e5); 
    std_v=TURB.stdev_v.*(60./1e5);

end

% --------------------------------

% Get random numbers
% Note, randn makes normally distibuted random numbers with zero mean and
% unit standard deviation. 
RNu=randn(size(nextXpos(validRows,1)));
RNv=randn(size(nextYpos(validRows,1)));

% Version of Carter's equation 3 
% note that at randn=1 and the current employed settings, the second term 
% results in about 140 m, with a max around 500 m for 10k random numbers
% ut and vt have units [km/min]
ut(validRows) = u0(validRows) .*(1-(dt./Tu)) +  std_u.*sqrt(2.*dt./Tu).*RNu;
vt(validRows) = v0(validRows) .*(1-(dt./Tv)) +  std_v.*sqrt(2.*dt./Tv).*RNv;

% add turb component in km units
nextXpos(validRows,1) = nextXpos(validRows,1) + ut(validRows).*dt;  
nextYpos(validRows,1) = nextYpos(validRows,1) + vt(validRows).*dt;

% disp(['Mean Total Added = ' num2str(sqrt( (mean_noNaN(abs(ut(:).*dt)).^2)  ...
%                + (mean_noNaN(abs(vt(:).*dt)).^2)  ).*1000) 'm'])

end
% -----------------------------------------------------------------------
function [X,Y,Ug,Vg] = grid_codar_data(HFx,HFy,U,V)
% GRID CODAR DATA
% [X,Y,Ug,Vg]=grid_codar_data(HFx,HFy,U_,V_)
% Turn arrays in to matricies, gridded in km space, creating the grid 
% from the inputs HFx, HFy. This makes the code for calc pos drifter much simpler.
% 
% WHy do this? For inputs to interp2 and interp3:
% "X and Y must be monotonic vectors or matrices produced by MESHGRID."

% Brian Emery 4 Jan 2010

% The 'round' here was used to fix the PWS codar grid. Shouldn't need it
% for the SBC grid or the ROMS PWS grid. SHouldn't need it for the PWS 
% grid based on the Codar software either.
%
% 
% keyboard
% 
% disp('** MINOR ROUNDING OF HF GRID **'), %keyboard
% [X,Y,Ug] = my_griddata(round(HFx),round(HFy),U);
% [X,Y,Vg] = my_griddata(round(HFx),round(HFy),V);
% 

[X,Y,Ug] = my_griddata(HFx,HFy,U);
[X,Y,Vg] = my_griddata(HFx,HFy,V);

% % some code to check 
% 
% % CHECK PLOTS
% close all
% figure
% for i = 55:size(U,2)
%     
%     hh=arrowplot(HFx, HFy,U(:,i),V(:,i),1/10);
% 
%     hold on
%     hn=arrowplot(X,Y,Ug(:,:,i),Vg(:,:,i),1/10,'r');
% 
%     pause(1)
%     delete([hh(:); hn(:)])
%     title(num2str(i))
% 
% end
% keyboard

end
% -----------------------------------------------------------------------
function CFG = cfg_check(CFG)
% CFG CHECK
% check for existence of fields and/or set defaults for b/w compatibility

% check for required fields
field_check(CFG,{'deploy_locations','deploy_times'})

% Field names to check, and set to defaults if not present
fn = {'maxDist'
      'drive'
      'eof_filled'
      'eof_data_path'
      'check_plots'
      'map_fxn'
      'map_axis'
      'ang_diff'
      'randn_state'
      'run_str'};


%  DEFAULT SETTINGS (USER DEFINABLE)
C.check_plots = 0;
C.ang_diff = [20 160];
C.maxDist = 20; 
C.drive = '/Data/totals/'; % <- usually unused. See get_total_data.m
C.eof_filled = 0;
%C.time_step = 1/24; % ?? this is set by the total data?
C.randn_state = sum(100*clock);
C.map_fxn = 'sbchan_map';
C.map_axis = [-120.4274 -119.1082   33.7805   34.6137];
C.eof_data_path =[];

% this is a string to add to the saved simulation data file name, as
% default use a time stamp
C.run_str = datestr(now,'_yyyymmdd_HHMM');

% Make sure CFG has all fields, use defaults for missing
for i = 1:numel(fn)
   if ~isfield(CFG,fn{i})
       CFG.(fn{i}) = C.(fn{i});
   end
end
end


% -----------------------------------------------------------------------
% RUNGA-KUTTA CODE
%
function [nextXpos,nextYpos] = calc_drifter_position(PDx,PDy,HFx,HFy,U,V,NT,serialtime)
% CALC DRIFTER POSITION.M - increment drifter position using runga-kutta
% [nextXpos,nextYpos] = calcPosDrifter(PDx,PDy,HFx,HFy,maxDist,U,V,NT,serialtime);
% Calculate drifter trajectories using fourth-order runge-kutta method
%
% INPUTS:
%   PDx, PDy = drifter positions (km) from central location
%    HFx,HFy = grid positions (km) from central location
%        U,V = 3-D CODAR velocity arrays (x,y,t)
%         NT = time index
% serialtime = decimal day

% Based on calcPosDrifter by Krisada Letcherononguong
% Modifications by Brian Emery 26jun03 
% Modified Aug 2007 by Kirk Ireson
%   - provide a vector in cases of low data availability, that is if any 
%     of the 9 boxes by 3 time frames (past, present, future) have a NaN 
%     value, then the calculation of a vector returns a NaN. It'd be better
%     in the case of low data availability to provide a value anyway
% Major Streamlining July 2009 by Brian Emery
% Version 8: Rewritten Jan 2010 to use gridded codar data. Now uses full 
% grid to do the interpolation instead of an arbitrary selection of grid 
% points

% % check plot code:
% hh=arrowplot(HFx,HFy,U(:,:,NT),V(:,:,NT),1/10,'k');
% hold on
% plot(PDx,PDy,'r.')

% Define constants
cm2km = 1/1e05;
dt = (serialtime(NT+1)-serialtime(NT)) * 24 * 3600;     % in seconds


% ---------------------------------------------------------
%  RUNGE-KUTTA 
% ---------------------------------------------------------
% use rugnge-kutta method to calculate the next position

% Use interp to get the best velocity value at the initial drifter position
u0 = interp2(HFx,HFy,U(:,:,NT),PDx,PDy);
v0 = interp2(HFx,HFy,V(:,:,NT),PDx,PDy);

% Try nearest neighbor if u0 or v0 is NaN. This appears to work if the
% drifter is within 1 grid point from the edge of coverage
if isnan(u0), u0 = interp2(HFx,HFy,U(:,:,NT),PDx,PDy,'nearest'); end
if isnan(v0), v0 = interp2(HFx,HFy,V(:,:,NT),PDx,PDy,'nearest'); end


%  STEP 1 ----------------
k1 = dt * u0*cm2km;         
l1 = dt * v0*cm2km;          

t1 = serialtime(NT)+0.5*(dt/(24*3600));
PDx1 = PDx+0.5*k1;
PDy1 = PDy+0.5*l1;

[u1,v1]=interp_velocities(PDx1,PDy1,t1,HFx,HFy,U,V,serialtime);

%  STEP 2 ----------------
k2 = dt*u1*cm2km;                   
l2 = dt*v1*cm2km;                    

t2 = serialtime(NT)+0.5*(dt/(24*3600));
PDx2 = PDx + 0.5*k2;
PDy2 = PDy + 0.5*l2;

[u2,v2]=interp_velocities(PDx2,PDy2,t2,HFx,HFy,U,V,serialtime);


%  STEP 3 ----------------
k3 = dt*u2*cm2km;                 
l3 = dt*v2*cm2km;                    

t3 = serialtime(NT)+(dt/(24*3600));
PDx3 = PDx + k3;
PDy3 = PDy + l3;

[u3,v3]=interp_velocities(PDx3,PDy3,t3,HFx,HFy,U,V,serialtime);

% STEP 4 ----------------
k4 = dt*u3*cm2km;             
l4 = dt*v3*cm2km;                

% PUT IT ALL TOGETHER -- 
nextXpos = PDx + 1/6*( k1 + 2*k2 + 2*k3 + k4);
nextYpos = PDy + 1/6*( l1 + 2*l2 + 2*l3 + l4);

end
% -----------------------------------------------------------------------
function [u1,v1]=interp_velocities(PDx1,PDy1,t1,x,y,U,V,t)
% INTERP VELOCITIES
% use interp3 to interpolate gridded 3-d codar data to a given position and
% time.
%
% INPUTS
% PDx1,PDy1,t1: the x,y and time (scalars) of a drifter
% x,y,U,V,t: the gridded (x,y 2-d, U,V 3-d with time t) data making the
% table.

% Brian Emery, Jan 4th 2010

% initialize outputs
u1 = NaN;
v1 = NaN; 

% TILE INPUTS
% time is a bit more tricky. The key is to expand it using repmat, then
% reshape it 3-d
x = repmat(x,[1 1 length(t)]);
y = repmat(y,[1 1 length(t)]);
t = reshape( repmat(t,[size(x,1)*size(x,2) 1]) ,[size(x,1) size(x,2) length(t)]);

% USE INTERP3
% x,y,t,u are the table, all must be same size
% PDx1,PDy1,t1 are the position to look up
% Define everythign explicitly ... this assumes something about the grid
% and number of neighbors ... the values in x and y must correspond to the
% values in u, and they must be gridded
if ~isnan(PDx1 + PDy1)
    try
        u1 = interp3(x,y,t,U,PDx1,PDy1,t1);
        v1 = interp3(x,y,t,V,PDx1,PDy1,t1);

        % use nearest neighbor if either is just outside the data coverage
        if isnan(u1), u1 = interp3(x,y,t,U,PDx1,PDy1,t1,'nearest'); end
        if isnan(v1), v1 = interp3(x,y,t,V,PDx1,PDy1,t1,'nearest'); end

    catch
        % Interp will choke on NaN position inputs, but there is otherwise
        % no catching going on. If we get here, u1,v1 are NaN and that just
        % meanst that the simulated drifter path will end
        fprintf('.')

    end
end

% % CHECK PLOTTING CODE
% plot(x(:,:,1),y(:,:,1),'.')
% hold on
% close all
% hh=arrowplot(x(:,:,1),y(:,:,1),U(:,:,NT),V(:,:,NT),1/10,'k');
% [PDx,PDy] = ginput(1)
% plot(PDx,PDy,'r.')
% hold on
% hh=arrowplot(x(:,:,1),y(:,:,1),U(:,:,NT),V(:,:,NT),1/10,'k');
% plot(x(:,:,1),y(:,:,1),'c.')

end



% -----------------------------------------------------------------------
%  DEBUG PLOTTING
%  Version uses m_map
%
function [cstr,hh,hd]=init_plot(maxDrifters,TURB,CFG,central_loc,drftX,drftY,ln)
% INIT_PLOT
% Used to initialize plots for debugging, etc
%
% TURB input could be optional


% evaluate option to plot internally
if ~CFG.check_plots
    [cstr,hh,hd]=deal(NaN(maxDrifters,1)); return
end

% Make the background figure 
figure
eval(CFG.map_fxn)
axis(CFG.map_axis)


% CREATE LINESTYLE ARRAY
% define symblols, colors
clr = ['b','g','r','c','m','y','k'];
sym = ['.'];%,'o','x','+','*','s','d','v','^','<','>','p','h'];

% create color/symbol arrays for plotting each drifter differently 
% max of 91 possible permutations if all clr and sym are used
cstr=[repmat(clr(:),length(sym),1) repmat(sym(:),length(clr),1)];

if exist('TURB','var')

    % use this option to make all the drifters in a cluster have the same
    % line style (ie same exact start positions)
    if TURB.n_clusters > length(cstr)
        cstr=repmat(cstr,ceil(TURB.n_clusters/size(cstr,1)),1);
    end
    cstr = cstr(1:TURB.n_clusters,:); 

    %keyboard

end

% Now tile cstr until so have enough for all the drifters
if maxDrifters > size(cstr,1)
    cstr=repmat(cstr,ceil(maxDrifters/size(cstr,1)),1);
end
cstr = cstr(1:maxDrifters,:);

% define codar vector handles
hh=[];


% ADD DRIFTER START LOCATIONS
hSt = drifter_plot(central_loc,drftX,drftY,cstr,{},CFG);
hleg = legend(hSt{1}(1),'Sim Start Locations');

hd={};

set(gcf,'units','inches','position',[1.3752    1.5754   11.6556    8.5581])
end
% -----------------------------------------------------------------------
function hd=drifter_plot(central_loc,nextXpos,nextYpos,cstr,hd,CFG)
% DRIFTER PLOT
% Add drifter to existing plot
%
% hd is the cell array of handles to the drifter plot objects. This allows
% older objects to disapear or turn different colors, etc

% evaluate option to plot internally
if ~CFG.check_plots
    hd=[]; return
end


[nextLon,nextLat] = km2lonlat(central_loc(1),central_loc(2),nextXpos,nextYpos);

jj = find(nextLon == 0); if ~isempty(jj), keyboard, end

if size(cstr,1) == 1
    % Simple case
    hd{end+1}=plot(nextLon,nextLat,cstr);
else
    % Case for many linestyles to be used
    h = NaN(size(cstr,1),1);
    for i= 1:size(cstr,1)
        h(i)=plot(nextLon(i,:),nextLat(i,:),['-' cstr(i,:)]);
    end
    hd{end+1}=h;
end

% Code to make the tails disappear ...
if numel(hd)>5
     delete(hd{1})
     hd=hd(2:end);
     if numel(hd)>10, keyboard, end
end

% % Add to title
% tstr=get(get(gca,'Title'),'String');
% title([tstr ]);

end
% -----------------------------------------------------------------------
function hh=codar_plot(hh,TUV,NT,validRows,CFG)
% DEBUG PLOTTING
% Plot the drifters, with HF data as the code runs

% Brian Emery

% evaluate option to plot internally
if ~CFG.check_plots
    hh=[]; return
end


% set scaling factors
sc_factor = 1/750;
limits = axis;
[ax,sc] = mercat(limits(1:2),limits(3:4));

if ~isempty('hh')
    delete(hh)
end


% Plot CODAR data
hh = arrowplot(TUV.LonLat(:,1), TUV.LonLat(:,2),TUV.U(:,NT),TUV.V(:,NT)*sc,sc_factor);
% hh = plotData( TUV, 'arrow',NT,sc_factor);

title(['Time: ' datestr(TUV.TimeStamp(NT),'ddmmmyyyy HH:MM') ' Hr Idx: ' num2str(NT) ' #Drifters: ' num2str(length(validRows))]);

% PAUSE!
pause(0.01)

end
% -----------------------------------------------------------------------
function spaghetti_plot(outLon,outLat,drifterT,maxDrifters,TURB,CFG)
% SPAGHETTI_PLOT
% Plot all of the tracks

% evaluate option to plot internally
if ~CFG.check_plots
    return
end

% Make the map, get color/symbol strings
if exist('TURB','var')
    cstr=init_plot(maxDrifters,TURB);
else
    cstr=init_plot(maxDrifters);
end

% Add drifters to the plot
for i = 1:maxDrifters
    plot(outLon(i,:),outLat(i,:),['-' cstr(i,:)])
end
title(['Drifter Tracks : ' datestr(drifterT(1)) ' to ' datestr(drifterT(end))])

end
% -----------------------------------------------------------------------
function last_check_plot(central_loc,drifterX,drifterY,maxDrifters,CFG)

% evaluate option to plot internally
if ~CFG.check_plots
    return
end

figure, pws_map, hold on
% this checks that the data in drifterLon,Lat is the same as above
[drifterLon,drifterLat] = km2lonlat(central_loc(1),central_loc(2),drifterX,drifterY);
for drifterIdx = 1:maxDrifters
    plot(drifterLon(drifterIdx,:),drifterLat(drifterIdx,:),'-b.')
end
ax=get_axis(drifterLon,drifterLat);
axis(ax)

end
% -----------------------------------------------------------------------
function ax=get_axis(x,y)
% SET_AXIS
% Set the axis based on data 
%
% x and y can be matrix inputs

x=[min(min(x))  max(max(x))];
y=[min(min(y))  max(max(y))];

dx=diff(x)*.3;
dy=diff(y)*.5;

n=3;
dx=n.*dx;
dy=n.*dy;
ax=[x(1)-dx x(2)+dx y(1)-dy y(2)+dy ];


end


