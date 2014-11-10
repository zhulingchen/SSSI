function Qave=qint2qave(Qint,t)
%
% 
% 
% Qint ... vector of interval Q's
% t    ... vector of times the same length as Qint
%NOTE: Qint(k) is the Q for the interval t(k)->t(k+1). This means the last
%Q applies to an undefined interval. We assume this interval to be the mean
%of diff(t).
% 
% Qave ... vector of average Q's the same length at t. Qave(k) is the
% average Q from t(1) to t(k). Qave(1) is identical to Qint(1).


%columnize
Qint=Qint(:);
t=t(:);
%
dt=diff(t);
dtm=mean(dt);
dtq=[dt;dtm];%tack on the last interval
tq=cumsum(dtq);
iQave=cumsum(dtq./Qint)./tq;
Qave=1./iQave;