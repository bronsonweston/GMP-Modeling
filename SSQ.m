function ssq = SSQ(parm)
%This function takes the concentration of [C/EBP]f and [PU.1] and returns
%the norm distance of the rate of change of these proteins from zero.
Cf=parm(1); P=parm(2);
[ dCdt, dPdt ] = dcfdp( Cf, P,0.01 );
ssq = norm([dCdt,dPdt]);
