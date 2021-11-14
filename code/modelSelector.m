function [rho, g, w, a] = modelSelector(model_name, num_clip)
%%% Function to output the specifications of the selected heli. model
%  Rr : Rotor length [cm]
%  Rw : A Rotor width [cm]
%  Bw : Body width [cm]
%  Tl : Tail length [cm]
%  Tw : Tail width [cm]
%  rho: Air density [g/cm^3]
%  clipmass : Mass of the clip[g]
%  paperrho : Paper density [g/cm^2]
%  g : Gravity acceleration [cm/s^2]
%  w : Paper helicopter weight [g]
%  a : Reference area [cm^2]
 
% All the lenght unit is cm.
switch model_name
    case 'S'
        %%% Small model specifications
        Rr = 72*10^(-1);    % Rotor length [cm]
        Rw = 32.4*10^(-1);  % A Rotor width [cm]
        Bw = 22.8*10^(-1);  % Body width [cm]
        Tl = 42*10^(-1);    % Tail length [cm]
        Tw = 30.6*10^(-1);  % Tail width [cm]
        rho= 0.001225;        % Air density [g/cm^3]
        clipmass = 0.595;  % [g]
        area_total = Rr*2*Rw+Bw*2*Rw+Tl*Tw;
        paperrho = 0.00768425; %g/cm^2
        g = 980;       % gravity acceleration [cm/s^2]
        w = (area_total*paperrho+clipmass*num_clip)*g;   % weight (g)
        a = pi*(Rr)^2;       % Ref. area [cm^2]
    case 'M'
        %%% Medium model specifications
        Rr = 96*10^(-1);    % Rotor length [cm]
        Rw = 43.2*10^(-1);  % A Rotor width [cm]
        Bw = 30.4*10^(-1);  % Body width [cm]
        Tl = 56*10^(-1);    % Tail length [cm]
        Tw = 40.8*10^(-1);  % Tail width [cm]
        rho = 0.001225;        % Air density [g/cm^3]
        clipmass = 0.595;  % [g]
        area_total = Rr*2*Rw+Bw*2*Rw+Tl*Tw;
        paperrho = 0.00768425; %g/cm^2
        g = 980;       % gravity acceleration [cm/s^2]
        w = (area_total*paperrho+clipmass*num_clip)*g;   % weight (g)
        a = pi*(Rr)^2;       % Ref. area [cm^2]
    case 'L'
        %%% Large model specifications
        Rr = 12;   % Rotor length [cm]
        Rw = 5.4;    % A Rotor width [cm]
        Bw = 3.8;    % Body width [cm]
        Tl = 7;    % Tail length [cm]
        Tw = 5.1;    % Tail width [cm]
        rho= 0.001225;        % Air density [g/cm^3]
        clipmass = 0.595;  % [g]
        area_total = Rr*2*Rw+Bw*2*Rw+Tl*Tw;
        paperrho = 0.00768425; %g/cm^2
        g = 980;       % gravity acceleration [cm/s^2]
        w = (area_total*paperrho+clipmass*num_clip)*g;   % weight (dyne)
        a = pi*(Rr)^2;       % Ref. area [cm^2]
           case 'H'
        %%% Helicopter model specifications from the following paper
        %%% "Teaching a Verification and Validation Course ... with Paper Helicopter” 
        %%% Chanyoung Park, Joo-ho Choi and Raphael T. Haftka
        Rr = 11.43;         % Rotor length [cm]
        % Rw = 5.08;          % A Rotor width [cm]
        % Bw = 3.81;          % Body width [cm]
        % Tl = 5.08;          % Tail length [cm]
        % Tw = 2.54;          % Tail width [cm]        
        %%% This two parameters(rho & g) are not exact because there isn't any exact
        %%% value for those on the paper
        rho= 0.0011839;  % Air density in Haftka paper [g/cm^3]    
        g = 980;       % gravity acceleration [cm/s^2]
        if num_clip == 1
            w = 2.05*g; % weight (dyne)
        elseif num_clip == 2
            w = 2.78*g; % weight (dyne)
        end
        a = pi*(Rr)^2;       % Ref. area [cm^2]       
end