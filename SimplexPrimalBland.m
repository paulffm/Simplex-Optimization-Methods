function [ xopt,B,message, iter, Zielfktnswert] = SimplexPrimalBland( A,b,c,Binit,xB)
%function [ xopt,B ,message, iter] = primalSimplex( A,b,c,B,xB )
%
% Primales Simplexverfahren mit Bland-Auswahlegel
%
% Input:  A, b, c - Daten für LP in primaler Standardform
%                      min c'x s.t. Ax=b, x>=0
%         B, xB   - primal zulaessige Basis, zugehörige Basislösung
%                   (optional)
% Output: xopt    - optimale Lsg
%         B       - zugehörige Basis
%         message - Information über Optimallösung oder Unbeschraenktheit
%         iter    - Anzahl der Iterationen
%
% Michel Reiffert, TU Darmstadt, Dezember 2017
tol=1e-6;
% Eingabefehler abfangen
if nargin < 4
    error('Nicht genug Eingabeargumente.');
end

[m, n] = size(A); %Einlesen der Dimensionen der Matrix
assert( m <= n );
if length(c)~=n || length(b)~=m
    error('Matrixdimensionen stimmen nicht ueberein.');
end

if length(Binit) ~= m
    error('Startbasis hat falsche Groesse.');
end

if rank(A) < m
    error('Matrix hat keinen vollen Zeilenrang.');
end

if rank(A(:,Binit)) < m
    error('B ist keine zulässige Startbasis.');
end

if nargin < 5
    xB=A(:,Binit)\b;
elseif norm(A(:,Binit)*xB-b,inf)>tol %%Maximumsnorm von (Ax-b) > Rechenungenauigkeit 
    size(xB,1)
    warning('xB ist keine Basisloesung! Bestimme Startloesung neu.');
    xB=A(:,Binit)\b;   
end
if min(xB)<-tol 
         error('B ist keine zulässige Startbasis.');
end

%Initialisierung
N=setdiff(1:n,Binit); %nimmt alle Indizes bis n, die nicht in B sind.
xopt=zeros(n,1); %Initialisiert mit n Nullen
iter=0;
message='';
Zielfktnswert=0;
B=Binit;

while 1 
    iter=iter+1; %in jedem Durchlauf iter hochzählen
    B;
    % (1) BTRAN:
    y=A(:,B)'\c(B);
    % (2) Pricing:
    zN=c(N)-A(:,N)'*y;
    
    % Falls zN nichtnegativ (innerhalb der Toleranz): Optimum gefunden
    if min(zN)>-tol
        xopt(B)=xB;
        B;
        Zielfktnswert = dot(c(B(1:m)),xB(1:m));
        message='Optimum gefunden';
        return %Aus der "while 1" Schleife kommt man nur mit return oder break raus
    end
    %Auswahlregel   
    %Bland-Regel:
    %wähle kleinsten Index j aus der Nichtbasis mit zj<0
    %(wir merken uns auch an welcher Stelle "pos j" in der Nichtbasis j steht
    %um später im Update an diese Stelle die verlassende Variable eintragen
    %zu können)
    [Nj,j]=min(N./(zN<=-tol)'); %(zN<=-tol) gibt Vektor mit 1en wenn der 
                               %zugehörige Eintrag aus zN negativ war, sonst 0.
                               %N./(zN<=-tol) damit entsprechend N. bei
                               %neg. zN oder inf bei positiven.
   
    % (3)FTRAN:
    w=A(:,B)\A(:,Nj);

    % (4)Ratiotest:
    %Bland-Regel
    %wähle kleinsten Index Bi aus der Basis mit wi>0
    
    %Suche positives wi
    l=1;
    while w(l)<tol %Suche den ersten positiven Eintrag von wi
        l=l+1;
        % Falls alle wi nichtpositiv (innerhalb der Toleranz): LP unbeschraenkt
        if l>m
            xopt(B)=xB;
            B;
            Zielfktnswert=inf;
            message='LP unbeschraenkt'     
            return
        end
    end
    
    %Suche unter verbleibenden Einträgen von wi:
    gamma=xB(l)/w(l);
    i=l;
    for l=i+1:m
        % einen Index i, sodass wi>0, xBi/wi minimal und Bi minimal
        if w(l)>tol && xB(l)/w(l)<gamma+tol &&( xB(l)/w(l)<gamma-tol || B(l)<B(i)) %1. Bedingung: w an der Stelle pos. UND
                                                                                   %2. Bedingung: Ratio kleiner als das bisher kleinste UND
                                                                                   %3. Entweder Ratio kleiner als bisher kleinste (dann muss ich es sowieso nehmen)
                                                                                   %   ODER für den Fall, dass der Ratio identisch ist (i.e. zwischen
                                                                                   %   gamma + - tol) nehme ich es falls der Index in der Basis kleiner ist. 
            gamma=xB(l)/w(l); %gamma updaten
            i=l; 
        end
    end
    
   %(5) Update:
    xB=xB-gamma*w;
    N(j)=B(i);
    B(i)=Nj;
    xB(i)=gamma;
    if iter > 200
        return
    end
end
