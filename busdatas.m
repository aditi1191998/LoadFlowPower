
function [bus_dat,nbs] = busdatas
% 101 - PQ Bus..
% 102 - PV Bus..
% 103 - Slack Bus..

%          |No | Type | V    | theta | Pg | Qg |   Pl  | Ql  | 
bus_dat = [ 1    103   1.050   0       0     0     0    00.0 ;
            2    102   1.020   0      25     0   05.0   03.0 ;
            3    101   1.000   0       0     0   27.5   06.5 ;
            4    101   1.000   0       0     0   00.0   00.0 ;
            5    101   1.000   0       0     0   15.0   06.0 ;
            6    101   1.000   0       0     0   20.0   02.5 ; ];

        [nbs, x] = size(bus_dat);
end