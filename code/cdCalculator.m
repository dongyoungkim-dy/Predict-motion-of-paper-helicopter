function cd = cdCalculator(eq_order,w,rho,a,h,g,t)
%%% The function to calculate the drag coefficient from the two falling
%%% equations
% w    : weight
% rho  : Air density
% c    : Drag coeff
% a    : Ref. area
% h    : Falling height
% g    : gravity acceleration
% t    : falling time
%%% The length unit is cm based
 
cd = zeros(length(t),1);
for i=1:1:length(t)
    switch eq_order
        case 'linear'
            if ~isnan(t(i))
                %%% Linear model equation
                tt = t(i);
                % Using fzero to solve nonlinear equation with x0
                f = @(c) (2*w/(rho*a*c*91.44))*tt+((2*w/(rho*a*c*91.44))^2)/g*(exp(-g/(2*w/(rho*a*c*91.44))*tt)-1)-h; 
                cd(i) = nan;
                x0 = 1;
                options = optimset('Display','off');
                cd(i) = fsolve(f,x0,options);
            else
               cd(i) = nan; 
 
            end
        case 'quadratic'
            if ~isnan(t(i))
                %%% Quadratic model equation        
                tt = t(i);
                % Using fzero to solve nonlinear equation with x0
%                 f = @(c) log((1+exp(2*g*tt/(sqrt(2*w/(rho*a*c)))))/2)*(sqrt(2*w/(rho*a*c)))^2/g-tt*(sqrt(2*w/(rho*a*c)))-h;
                f = @(c) log((1+exp(-2*g*tt/(sqrt(2*w/(rho*a*c)))))/2)*((sqrt(2*w/(rho*a*c)))^2)/g+tt*(sqrt(2*w/(rho*a*c)))-h;
                cd(i) = nan;
                x0 = 1;
                options = optimset('Display','off');
                cd(i) = fsolve(f,x0,options);
            else
               cd(i) = nan; 
            end            
    end
end