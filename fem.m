% solve -u’’ + u =1 on (0,1) with u(0)=0 and u’(1)=0 using FEM
N=1000; % number of elements
x=linspace(0,1,N+1); % set up mesh
h=1/N; % h=constant mesh size
K=zeros(N+1,N+1); % set up empty stiffness matrix
F=zeros(N+1,1); % set up global load vector
for i=1:N % loop over elements
Ktemp = [1 -1;-1 1]/h + [2 1;1 2]*(h/6); % compute entries in element stiffness matrix
Ftemp = [h/2;h/2]; % compute entries in element load vector
K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+Ktemp; % add to global stiffness matrix
F(i:i+1)=F(i:i+1)+Ftemp; % add to global load vector
end
K(1,:)=[]; % apply Dirichlet boundary condition by
K(:,1)=[]; % removing first row and column of stiffness matrix
F(1)=[]; % and first row of load vector
U=K\F; % compute solution to matrix problem
U=[0;U]; % add value at x=0 to solution vector
uExact = dsolve('-D2u + u = 1', 'u(0)=0', 'D1u(1)=0', 'x'); %find the exact solution 
P = zeros(N+1,1);
for j=1:N+1
P(j,1) = -(U(j)-(1 - (exp(2)*exp(-x(j))/(exp(2) + 1) - exp(x(j))/(exp(2) + 1))));
end
disp(P);
%disp(x)
plot(x,P,'rx-')
%plot(x,U,'rx-') % plot solution
%hold
%hold on
%syms x
%di = char(x-0.5*(x.^2));
%di = char(subs(uExact));
%H = zeros(N,1);
%E = zeros(N,1);
%plot(x, subs(uExact), 'b', 'LineWidth', 2)
%for n=1:N
%H(n,1) = 1/n;
%E(n,1) = (H(n,1).^2)/(pi);
%end
%plot(H,E,'r');
%disp(E)
%plot(x, subs(uExact-U), 'b', 'LineWidth', 2)
%fplot(di,[0,1],'b') % plot exact solution
%fplot(di,[0,2],'b')
%disp(U)
%disp(K)