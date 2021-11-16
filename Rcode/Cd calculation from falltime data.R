
falltime = read.csv("final_test_data.csv", header = TRUE)
falltime_linear = read.csv("falltime_linear.csv", header = TRUE)
t = falltime[,1]


# w    : weight
# rho  : Air density
# c    : Drag coeff
# a    : Ref. area
# h    : Falling height
# g    : gravity acceleration
# t    : falling time
### The length unit is cm based

Rr = 72*10^(-1);    # Rotor length [cm]
num_clip = 1
Rw = 32.4*10^(-1);  # A Rotor width [cm]
Bw = 22.8*10^(-1);  # Body width [cm]
Tl = 42*10^(-1);    # Tail length [cm]
Tw = 30.6*10^(-1);  # Tail width [cm]
rho= 0.001225;      # Air density [g/cm^3]
clipmass = 0.595;   # [g]
h=10.67*100
area_total = Rr*2*Rw+Bw*2*Rw+Tl*Tw;
paperrho = 0.00768425; #g/cm^2
g = 980;       # gravity acceleration [cm/s^2]
w = (area_total*paperrho+clipmass*num_clip)*g;   # weight (g)
a = pi*(Rr)^2;       # Ref. area [cm^2] 
#time = readtable('test_falltime_data.csv');
t = as.data.frame(t)
cd = numeric(1)
#for i=1:1:length(t)

#if ~isnan(t(i))
## Linear model equation
#tt = t(i);

dslnex <- function(c) {
  for(i in 1:length(t)){
    tt = t[i]
    #c = numeric(length(t))
    #y <- numeric(length(t))
    y <- (2*w/(rho*a*c*91.44))*tt+((2*w/(rho*a*c*91.44))^2)/g*(exp(-g/(2*w/(rho*a*c*91.44))*tt)-1)-h 
  }
}

xstart <- c(2)
fstart <- dslnex(xstart)
nleqslv(xstart, dslnex, control=list(btol=.01))
