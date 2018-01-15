function ssq = SSQ2(parm)
global Var constantP
%This function serves as an alternative solving function when the
%bifurcation solution makes a turn. It serves the same role as SSQ, however,
%it allows for the variation in the  bifurcation parameter rather than
%[PU.1]. It returns the norm distance of the rate of change of these 
%proteins from zero.
Cf=parm(1); ChangeVar( Var, parm(2) )
[ dCdt, dPdt ] = dcfdp( Cf, constantP,0.01 );
ssq = norm([dCdt,dPdt]);
