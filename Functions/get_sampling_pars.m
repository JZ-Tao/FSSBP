function pars = get_sampling_pars(ratio, interplator)
if ~exist('interplator','var')
   interplator = 'tap23';
end


% The default interplator is tap23.
if strcmp(interplator, 'tap23')
    pars.down.using_imresize = 1;
    % not used. Just for clarifing that imresize+nearest == downsample+offset(2,2) 
    pars.down.offset = 2; 
    pars.up.interp_type = 0; % not used. Just for clarification
    pars.up.offset = [1,0]; 
    pars.up.tap = 23;

elseif strcmp(interplator, 'tap7')
    pars.down.using_imresize = 1;
    % not used. Just for clarifing that imresize+nearest == downsample+offset(2,2) 
    pars.down.offset = 2; 
    pars.up.interp_type = 0; % not used. Just for clarification
    pars.up.offset = [1,0]; 
    pars.up.tap = 7;

elseif strcmp(interplator, 'tap8')
    pars.down.using_imresize = 0;
    pars.down.offset = 1;
    pars.up.interp_type = 3; % not used. Just for clarification
    pars.up.offset = [1,1]; 
    pars.up.tap = 8;

elseif strcmp(interplator, 'general')
    pars.down.using_imresize = 0;
    pars.down.offset = floor(ratio/2); %[1,1]; 
    pars.up.interp_type = 1; % not used. Just for clarification
    pars.up.offset = floor(ratio/2); 
    pars.up.tap = 33; 

elseif strcmp(interplator, 'bicubic')
    pars.down.using_imresize = 1;
    pars.down.offset = floor(ratio/2); %[1,1]; 
    pars.up.offset = floor(ratio/2); 
    pars.up.tap = 7; 
    pars.up.interp_type = 1;

end

pars.up.interplator = interplator;
