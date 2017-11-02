# Automatic PWA Approximation Toolbox #

*Warning:* the code is experimental and subject to frequent changes.

## Prerequisites ##

* MATLAB
* Symbolic toolbox
* Optimization toolbox
* [HYSDEL](http://control.ee.ethz.ch/~hybrid/hysdel/hysdel.php)

## Download ##

The `autoprox` package is available in the [Downloads](/kvasnica/autoprox/downloads) section.

## Bug reports and inquiries ##

Send your comments, questions, and/or bug reports to `michal.kvasnica(at)stuba.sk`

## Further reading ##

See Chapter 3 of [Alex Szucs' PhD thesis](http://www.kirp.chtf.stuba.sk/publication_info.php?id_pub=1522)

## Quick start ##

Download the ZIP package, unpack, and run one of the following demos:

* `demo_1d`
* `demo_2d`
* `demo_nd`
* `demo_nd_sum`
* `demo_deq`

or perform the approximations by the graphical user interface, which can be launched by the following command: `autoprox_gui`

## Usage - command line interface ##

### Approximation of functions of one variable ###

```matlab
syms x
f = x^3;
f_bounds = [-1, 1];
N_segments = 3;
[f_aprx, f_aprx_data] = autoprox_1d(f, f_bounds, N_segments);

% export the approximation to a HYSDEL file
hysdel_1d(f_aprx_data, 'f_aprx.hys');

% plot the original function and its approximation
x = linspace(-1, 1, 100);
plot(x, x.^3, x, f_aprx(x), 'r--');
```

### Approximation of the product of two functions of one variable ###

```matlab
% approximate x1^3 * (abs(x2) + 0.5*x2^2 - sin(x2^3))

syms x1 x2
f1 = x1^3;
f2 = abs(x2) + 0.5*x2^2 - sin(x2^3);

% domains of f1 and f2
f1_bounds = [-1.5, 1.5];
f2_bounds = [-1, 2.5];

% number of approximation segments
f1_segments = 3;
f2_segments = 7;
y1_segments = 3;  % approximation of (f1+f2)^2
y2_segments = 3;  % approximation of (f1-f2)^2

[f_aprx, f_aprx_data] = autoprox_2d(f1, f2, f1_bounds, f2_bounds, ...
    f1_segments, f2_segments, y1_segments, y2_segments);

% export the approximation to a HYSDEL file
hysdel_2d(f_aprx_data, 'f1f2_aprx.hys');

% plot the original function and the approximation
f = f1*f2;
fm = matlabFunction(f);

[x1, x2] = meshgrid(linspace(f1_bounds(1), f1_bounds(2), 100), ...
    linspace(f2_bounds(1), f2_bounds(2), 100));
subplot(1, 2, 1);
surf(x1, x2, fm(x1, x2));
axis tight
title('Original function.');
subplot(1, 2, 2);
surf(x1, x2, f_aprx(x1, x2));
axis tight
title('Approximation.');
```

### Approximation of general N-D functions ###

```matlab
clear
syms z1 z2 z3 z4 z5 z6 

% the symbolic function as a string
fun = 'sin(z1)*z2*z6+z1*log(z3)*z4*sin(z5)+cos(z3)*z4*abs(z5)*cos(z6)';

% vector of symbolic variobles involved in the function
vars = [z1 z2 z3 z4 z5 z6];

% bounds for each variable
xbmin = [-1.5 -1.5 1.5 -1.5 -1.5 -1.5];
xbmax = [1.5 1.5 2.5 1.5 1.5 1.5];

% number of approximation points
Nlin_x1 = [4 2 3];
Nlin_u1 = [3 5 4 2];

Nlin_x2 = [2 4 3 2];
Nlin_u2 = [2 5 2 4 3 3];

Nlin_x3 = [2 4 5 3];
Nlin_u3 = [4 2 2 3 4 3];

Nlin_x = {Nlin_x1 Nlin_x2 Nlin_x3};
Nlin_u = {Nlin_u1 Nlin_u2 Nlin_u3};

% obtain the approximation data
apprx_data = approximation_nd_sum(fun,vars,xbmin,xbmax,Nlin_x,Nlin_u); 
           
% generate the corresponding HYSDEL file
hysdel_nd_sum(apprx_data, 'FNsum_approx.hys'); 

% compile the HYSDEL file to obtain a simulator
hysdel('FNsum_approx.hys', 'FNsum_approx_sim.m');

% obtain values of f_aprx() for a given input data
z_data = [1; 0.1; 2; 0; 0; 0];
[x,d,z,y] = FNsum_approx_sim([], z_data);
y
```

### Approximation of a nonlinear vector field f(x,u) ###

```matlab
% symbolic variables representing the state and input variables
syms ca v vc qc

vars_state = [ca v vc];
vars_input = [qc];

% right-hand sides of the differential equations 

right_sides{1} = '0.2752-0.0652*ca*exp(-v)-v*ca';
right_sides{2} = '20.4130-3.4338*v*ca+0.1337*v+0.0685*vc';
right_sides{3} = '0.25*qc + 0.0685*qc*vc-0.275*vc';

% number of approximation segments for each nonlinear function 

Nlin_x1 = [2 2];
Nlin_u1 = [2 2];

Nlin_x2 = [2 2];
Nlin_u2 = [2 2];

Nlin_x3 = [2 2];
Nlin_u3 = [2 2];

Nlin_x4 = [2 2];
Nlin_u4 = [2 2];

Nlin_x = {Nlin_x1 Nlin_x2 Nlin_x3 Nlin_x4};
Nlin_u = {Nlin_u1 Nlin_u2 Nlin_u3 Nlin_u4};

% lower and upper bounds for the state and input variables, respectively

% states
xbmin = [-3 -3 -3];
xbmax = [3 3 3];

% inputs 
umin = -4;
umax = 4;

% sampling time
Ts = 2;

% obtain the approximation data
apprx_data = autoprox_deq(vars_state,vars_input,xbmin,xbmax,Nlin_x,Nlin_u,umin,umax,right_sides,Ts);

% generate the corresponding HYSDEL file
hysdel_nd_sum(apprx_data, 'diffeq_approx.hys',apprx_data{end}); 

% compile the HYSDEL file to obtain a simulator
hysdel('diffeq_approx.hys', 'diffeq_approx_sim.m','-5');

% obtain values of f_aprx() for a given input data
x_data = [1;0.1;2];
u_data = [2];
[x,d,z,y] = diffeq_approx_sim(x_data,u_data);
```

## Usage - graphical user interface ##

### Brief description ###

The basic window consists of five sections each serving for different tasks:

* *Approximation* - determine the dimension of the approximation
* *Setup* - define the nonlinear function and its bounds and also to plot the original nonlinear term
* *Generate PWA approximation* - to render the approximation or to split the the original function, if approximation of an n-dimensional functions has to be realized
* *Approximation summary* - comprehensive description about the efficiency of the linearization, including several indicators of quality(average error, worst case error ...)
* *Post-processing* - serving for exporting the acquired data to teh HYSDEL language