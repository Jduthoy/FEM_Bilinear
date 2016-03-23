function [ e ] = Bilinear2D( H,g,he,ve,ln,f )
%function [ F ] = Bilinear2D( 0,0,3,4,4,@(x,y) 2*(x.^2-x+y.^2-2*y) )
%%%%%%% Julius Duthoy%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Math 520     %%%%%%%%%%%%%%%%%%%%%%%%%%
% Uxx+Uyy=2*(x^2-x+y^2-2y)
%BCs: u(x,y)=H=0    d/dn(u(x,y))=g=0
%%%%%%%%%%List of Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%ln  =# Local Nodes
%gn  =# Global Nodes 
%ne  =# elements- input into the function
%n_eq=# Nonzero Global Basis Fncts
%he=horizontal elements
%ve=vertical elements

ne=ve*he;
gn=(ve+1)*(he+1);

%%%%%%%%%% Set up Global Nodes & Constant Interval h %%%%%%%%%%%%%%%%%%
[y,x] = meshgrid(0:2/ve:2,0:1/he:1);% for the y mesh used 1/2 manually instead of 1/ve 'cuz david
h_y=y(:,2)-y(:,1);
h_x=x(2,:)-x(1,:);
[m,n]=size(x);


i=1:gn;
grph=zeros(he+1,ve+1);
grph(:)=i;
grph=grph';

%%%%%%%%%% ID(A)=P Matrix%%%%%%%%%%%%%%%%%%%%%%%%%% 

n_eq=0;%Determined from ID Matrix below

for A=1:gn
    if ismember(A,grph(:,1))==1 || ismember(A,grph(1,:))==1 || ismember(A,grph(:,he+1))==1 || ismember(A,grph(ve+1,:))==1
        ID(A)=0; 
    else
        n_eq=n_eq+1;
        ID(A)=n_eq;
    end
end
%ID;

%%%%%% IEN(a,e)=A Matrix %%%%%%%%%%%%%%%%%%%%%%%%
B1=grph(1:ve,1:he);
B1=B1';
B2=grph(1:ve,2:he+1);
B2=B2';
B3=grph(2:ve+1,2:he+1);
B3=B3';
B4=grph(2:ve+1,1:he);
B4=B4';
C1=reshape(B1,1,ne);
C2=reshape(B2,1,ne);
C3=reshape(B3,1,ne);
C4=reshape(B4,1,ne);
IEN=zeros(ln,ne);
IEN(:)=[C1; C2; C3; C4];

%%%%%%%%% LM(a,e)=ID(IEN) Matrix %%%%%%%%%%%%%%%%%%%%%%%%
LM=ID(IEN);


%%%%%%% Initialize Local Basis fncts %%%%%%%%%%%%%%%%%%%%%%%%% 

N1=@(x,y) 1-y-x+x*y;
dN1_x=@(x,y) -1+y+x-x; %had to add in weird x-x stuff cuz of how quad works
dN1_y=@(x,y) -1+x+y-y;

N2=@(x,y) x-x*y;
dN2_x=@(x,y) 1-y+x-x;
dN2_y=@(x,y) -x+y-y;

N3=@(x,y) y-x*y;
dN3_x=@(x,y) -y+x-x;
dN3_y=@(x,y) 1-x+y-y;

N4=@(x,y) x*y;
dN4_x=@(x,y) y+x-x;
dN4_y=@(x,y) x+y-y;

N={'1-t-z+z.*t','z-z.*t','z.*t','t-z.*t'};


dN_x={'t.^2-2*t+1+z-z','-t.^2+2*t-1+z-z','t.^2-t+z-z','-t.^2+t+z-z';
      '-t.^2+2*t-1+z-z','t.^2-2*t+1+z-z','-t.^2+t+z-z','t.^2-t+z-z';
      't.^2-t+z-z','t.^2+t+z-z','t.^2+z-z','-t.^2+z-z';
      '-t.^2+t+z-z','t.^2-t+z-z','-t.^2+z-z','t.^2+z-z'};% x-derivatives multiplied by each other
dN_y={'z.^2-2*z+1+t-t','-z.^2+z+t-t','z.^2-z+t-t','-z.^2+2*z-1+t-t';
      '-z.^2+z+t-t','z.^2+t-t','-z.^2+t-t','z.^2-z+t-t';
      'z.^2-z+t-t','-z.^2+t-t','z.^2+t-t','-z.^2+z+t-t'
      '-z.^2+2*z-1+t-t','z.^2-z+t-t','-z.^2+z+t-t','z.^2-2*z+1+t-t'};


%%%%%%%%%Stiffness Matrix and F Vector%%%%%%%%%%%%%%%%%%%%%%%%%
K=zeros(n_eq,n_eq); 
F=zeros(n_eq,1);

%fnct=inline(N{1,2});
%F(LM(2,1))=F(LM(2,1))+h(1).*quad(@(z) f(x(LM(2,1))+h(1).*z).*fnct(z), 0, 1)

%e=0.1048 for 12 elements
%e=0.2438 for 1200 elements
k =[-0.6667   -0.1667    0.3333    0.1667
   -0.1667   -0.6667    0.1667    0.3333
    0.3333    0.1667   -0.6667    0.1667
    0.1667    0.3333    0.1667   -0.6667];


for e=1:ne
       %Set Up Big K
       for a=1:ln
            if LM(a,e)~=0
                    fnct=inline(N{1,a});
                    F(LM(a,e))=F(LM(a,e))+(h_x(a))*(h_y(a)).*dblquad(@(z,t) f(x(IEN(a,e))+h_x(a)*z,y(IEN(a,e))+h_y(a)).*fnct(z,t), 0, 1,0,1);
                for b=1:ln
                     der_x=inline(dN_x{a,b});
                     der_y=inline(dN_y{a,b});
                     %k(a,b)= -1.*dblquad(@(z,t) der_x(z,t) + der_y(z,t),0,1,0,1);
                    if LM(b,e)~=0
                        K(LM(a,e),LM(b,e))=K(LM(a,e),LM(b,e))+k(a,b);
                    end
                end
            end
       end         
end



% Solve for d

d = K\F;

% Get u from d
for A=1:gn
    if ID(A)~=0
        u(A)=d(ID(A));
    else
        u(A)=g;
    end
end
      
%%%%%%%%%%Comparison to Exact Soln%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% zero BC's Solution%%%%%%%%%%%%%%%%%%%%%%%
%{
sol= @(x,y) (1-x)*x*(2-y)*y;
%%%%%%%% Actual BC's of the problem%%%%%%%%%%%%%%%%%%%%%%%%
%sol=@(z)-1*(2*exp(-z/2)*(-34*exp(z/2)*z*sin(8)-28732*sin(8)*sin(2*z)+1666*exp(2)*sin(2*z)+42*exp(z/2)*sin(8)...
%    +136*exp(z/2)*z*cos(8)-28732*cos(8)*cos(2*z)-168*exp(z/2)*cos(8)- 7183*cos(8)*sin(2*z)+7183*sin(8)*cos(2*z)))/...
%    (289*(4*cos(8)-sin(8)));
UU=arrayfun(sol, x, y);%UU=Actual u calculated by hand
%arrayfun(function,inputs) uses x array as inputs to sol to form the actual
UU=reshape(UU,[1,gn]);
x=reshape(x,[1,gn]);
y=reshape(y,[1,gn]);
%%%%%%%%%% Plots the Actual U and my FEM approximation
plot3(x,y, u, 'm*', x,y, UU, 'k')
legend('Approx u', 'Actual u')
str=sprintf('Uxx+Uyy, BC:u(x,y)=%f, du/dn(x,y)=%f',g,H);
% Couldn't get f to change to whatever the user inputs
title(str);
xlabel('x')
ylabel('y')
zlabel('u')
%}

figure(1);
U = reshape(u,[m,n]);
surf(x,y,U)
figure(2);
surf(x,y,(1-x).*x.*(2-y).*y)



%format long
%e=max(abs(UU-u)); %Calculates the error
e=max(max(abs(U-((1-x).*x.*(2-y).*y))));

%}

end

