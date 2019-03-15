% Program for Newton-Raphson Load Flow Analysis..

[bus_dat,nbs] = busdatas;      % Calling busdatas..
Y = ybus;          
BMva = 50;                     % Base MVA..
bus = bus_dat(:,1);            % Bus Number..
type = bus_dat(:,2);           % Type of Bus
V = bus_dat(:,3);              % Specified Voltage..
del = bus_dat(:,4);            % Voltage Angle..
Pg = bus_dat(:,5)/BMva;        
Qg = bus_dat(:,6)/BMva;        
Pl = bus_dat(:,7)/BMva;       
Ql = bus_dat(:,8)/BMva;       
P = Pg - Pl;                
Q = Qg - Ql;                
Psp = P;                    
Qsp = Q;                    
G = real(Y);                % Conductance matrix..
B = imag(Y);                % Susceptance matrix..

pv = find(type == 102 | type == 103);   % PV Buses..
pq = find(type == 101);                 % PQ Buses..
npv = length(pv);                       % No. of PV buses..
npq = length(pq);                       % No. of PQ buses..

Tol = 1; 
Iter = 1;
while (Tol > 1e-3)   % Iteration starting..
    
    P = zeros(nbs,1);
    Q = zeros(nbs,1);
    % Calculate P and Q
    for i = 1:nbs
        for k = 1:nbs
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbs
        if type(i) == 101
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbs);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbs-1,nbs-1); %1 slack bus removed
    for i = 1:(nbs-1)
        m = i+1;
        for k = 1:(nbs-1)
            n = k+1;
            if n == m
                for n = 1:nbs
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    % Vm not multiplied 
    J2 = zeros(nbs-1,npq);
    for i = 1:(nbs-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbs
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbs-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbs-1)
            n = k+1;
            if n == m
                for n = 1:nbs
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbs
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];    % Jacobian Matrix..
    X = gauss_elimination(J,M);
    %X = inv(J)*M           % Correction Vector
    dTh = X(1:nbs-1);      % Change in Voltage Angle..
    dV = X(nbs:end);       % Change in Voltage Magnitude..
    
    % Updating State Vectors..
    del(2:nbs) = dTh + del(2:nbs);    % Voltage Angle..
    k = 1;
    for i = 2:nbs
        if type(i) == 101
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end
    
    Iter = Iter + 1
    Tol = max(abs(M));                  % Tolerance..
    
end
loadflow(V,del,BMva);              % Calling Loadflow.m..