% code to setup bipolar environment and data paths

currentPath = pwd;
parentPath = fileparts(currentPath);
%setenv('BIPOLAR_DATA', fullfile(currentPath, 'data'));
setenv('BIPOLAR_DATA', '/Users/davidcaldwell/Box/KLEENLAB/David/Data');
addpath(genpath(fullfile(currentPath)))
