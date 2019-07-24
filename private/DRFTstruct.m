function T = DRFTstruct(DIM)
% DRFT STRUCT - standardized structure for drifter data
% T = DRFTstruct(DIM)
% 
% Basically a customized copy of HFRP TRAJstruct.m, keeping only parts,
% to initialize a standard format for drifter data
% 
% INPUTS
%  
% DIM: Size of trajectory matrices. A two element vector
% with [ NGridPts, NTimeStamps ].  Defaults to [ 0 0 ].
%  
% OUTPUTS:
%
% DRFT:  the empty drifter data structure 
%
% DRFT STRUCTURE FORMAT NOTES
% There is some variation in drifter structures. These can contain:
% - one element for each drifter (lon lat stime all 1 row, cols = time)
% - one element for each cluster with n drifters (size(lon) = n x time)
% - one element for each bunch of drifters associated with the same time
%   array (the key is that the time is assumed to apply to each column), 
%   such as for the output of drifter_simulation.m
%
% ALL OF THESE will contain fields:
% TimeStamp
% Lon
% Lat

% Copyright (C) 2010  Brian Emery
% 20 Feb 2010
% Custom version to add some useful meta data fields and remove some fields
% that I dont (yet) need.

if ~exist( 'DIM' , 'var' )
  DIM = [0 0];
end

% General
T.Type = ['Drifter Data'];
T.ClusterName = '';
T.TrajectoryDomain = '';
T.CreationInfo = creation_info;
T.GridOrigin =[];
T.SiteOrigin =[];

% Time
T.TimeStamp = repmat(NaN,[1,DIM(2)]);
T.FixTimes = repmat(NaN,[1,DIM(2)]);
T.TimeZone = 'GMT';

% T.CreateTimeStamp = datestr(now);

% Space
[T.Lon,T.Lat] = deal( repmat( NaN, DIM ) );
[T.X,T.Y] = deal( repmat( NaN, DIM ) );

% Other
[T.U,T.V] = deal( repmat( NaN, DIM ) );

% Some units
T.LonLatUnits = 'Decimal Degrees';
T.GridUnits = 'km';
T.UVUnits = 'cm/s';
T.XYUnits = 'km';

% % Other
% T.OtherMatrixVars = []; % Should eventually be a structure.
% T.OtherSpatialVars = []; % Should eventually be a structure.
% T.OtherTemporalVars = []; % Should eventually be a structure.
T.Hourly = []; % For hourly averaged data
T.OtherMetadata = []; % Should eventually be a structure.


T.ProcessingSteps = {};

end