function [ xopt,B,message, iter, Zielfktnswert] = SimplexDantzigPricing2019( A,b,c,Binit,xB)
%function [ xopt,B ,message, iter] = primalSimplex( A,b,c,B,xB )
%
% Primales Simplexverfahren mit Dantzig Auswahlregel im Pricing
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
    %Auswahlregel nach Dantzig:   
    
    [~,j]=min(zN)
    Nj=N(j)
   
    % (3)FTRAN:
    w=A(:,B)\A(:,Nj)

    % (4)Ratiotest:
     %---- Alternativer Ratiotest ----
     %Bestimme Indizes mit positiven Einträgen in w
    wpos=find(w>=tol);
    
   %Falls es keine gibt: LP unbeschränkt
     if numel(wpos)==0
             xopt(B)=xB;
             message='LP unbeschraenkt';
             return
     end  
     
     %gamma ist das Minimum der Quotienten für positive w-Einträge
    gamma=min(xB(wpos)./w(wpos));
     
     %Indizes i, sodass xBi/wi=gamma
     isgamma=wpos((xB(wpos)./w(wpos))<=gamma+tol);
     
     %kleinstes Bi, sodass xBi/wi=gamma
     [~,index]=min(B(isgamma));
     i=isgamma(index);
     gamma=xB(i)/w(i);
  
    %Update
    xB=xB-gamma*w
    N(j)=B(i)
    B(i)=Nj
    xB(i)=gamma
    if iter > 10
        return
    end
    
     
end
