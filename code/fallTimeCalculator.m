
Rr = 72*10^(-1);    % Rotor length [cm]
num_clip = 2
Rw = 32.4*10^(-1);  % A Rotor width [cm]
Bw = 22.8*10^(-1);  % Body width [cm]
Tl = 42*10^(-1);    % Tail length [cm]
Tw = 30.6*10^(-1);  % Tail width [cm]
rho= 0.001225;        % Air density [g/cm^3]
clipmass = 0.595;  % [g]
h=10.67*100
area_total = Rr*2*Rw+Bw*2*Rw+Tl*Tw;
paperrho = 0.00768425; %g/cm^2
g = 980;       % gravity acceleration [cm/s^2]
w = (area_total*paperrho+clipmass*num_clip)*g;   % weight (g)
a = pi*(Rr)^2;       % Ref. area [cm^2] 
T = readtable('myfile.csv');
function falltime_data = fallTimeCalculator(eq_order,w,rho,cd,a,h,g)
%%% The function to calculate the fall time from the two falling equations
% w    : weight
% rho  : Air density
% c    : Drag coeff
% a    : Ref. area
% h    : Falling height
% g    : gravity acceleration
% t    : falling time
%%% The length unit is cm based

falltime_data = zeros(length(cd),1);
for i=1:length(cd)
    switch eq_order    
        case 'linear'
            %%% Using fsolve to solve nonlinear equation with x0
            c = cd(i);           
            Vss = 2*w/(rho*a*c*91.44);
            f = @(tt) Vss*tt+((Vss)^2)/g*(exp(-g/Vss*tt)-1)-h; 
            x0 = 1;           
            options = optimset('Display','off');  % Turn off display
            falltime_data(i) = fsolve(f,x0,options);
            while ~isreal(falltime_data(i))
                falltime_data(i) = fsolve(f,x0,options);
                x0 = x0+0.2;
            end
 
   
        case 'quadratic'
            %%% Exact solution
            Vss=sqrt(2*w/(rho*a*cd(i)));    % Terminal velocity
            falltime_data(i) = Vss/g*acosh(exp(g*h/(Vss^2)));  % Falling time
            
            %%% Using fsolve to solve nonlinear equation with x0
%             c = cd(i);           
%             Vss = sqrt(2*w/(rho*a*c));
% %             f = @(tt) log((1+exp(2*g*tt/Vss))/2)*(Vss)^2/g-tt*Vss-h;
%             f = @(tt) log((1+exp(-2*g*tt/Vss))/2)*(Vss)^2/g+tt*Vss-h;
%             x0 = 1.5;                        
%             options = optimset('Display','off');  % Turn off display
%             falltime_data(i) = fsolve(f,x0,options);
%             while ~isreal(falltime_data(i))
%                 falltime_data(i) = fsolve(f,x0,options);
%                 x0 = x0+0.2;
%             end
 
                
    end
end