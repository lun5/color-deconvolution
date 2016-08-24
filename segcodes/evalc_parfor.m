%% call evalc in a parfor loop
% Luong Nguyen 09/21/2015
% INPUT: command: string containing the system command
% OUTPUT: printout of the system command
function s = evalc_parfor(command)
    s = evalc(['system(''' command ''')']);
end