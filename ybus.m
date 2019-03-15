
function Y = ybus  % Returns Y

line_dat = linedatas;           
fb = line_dat(:,1);             
tb = line_dat(:,2);             
r = line_dat(:,3);              
x = line_dat(:,4);              
b = line_dat(:,5);              
a = line_dat(:,6);           
z = r + 1i*x;
y = 1./z; 
b = 1i*b./2;                     
[c, nbs] = busdatas; 
nl = length(fb);                  % No. of branches...
Y = zeros(nbs,nbs);               % Initialise YBus...
 
 % Formation of the Off Diagonal Elements...
 for k = 1:nl
     Y(fb(k),tb(k)) = - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end
 % Formation of Diagonal Elements....
 for m = 1:nbs
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end