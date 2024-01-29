# VDX: Automatic Index Tracking for casADi NLPs
VDX is a wrapper for the creation of generic casADi NLPs which provides automatic, transparent, index tracking facilities for named variable classes.
It also provides a direct interface to casADi nlpsol and allows for automatic results extraction for the mentioned varable classes.
## Quickstart Guide for Users
This section provides a quick overview for general users of how to use this tool.
First we create a problem via:
```
prob = vdx.Problem();
```
which creates a problem of the form
```
minimize f(w,p), subject to
lbw <= w <= ubw
lbg <= g(w,p) <= ubg
```
The decision variable `w`, parameters `p`, and the constraint function `g`, are subdivided into individual index tracked "variable classes".
These are useful to track the semantic meanings of variables or contstraints.

### Creating Variables
The following example produces the `x` and `y` variables of hanging masses in a hanging chain problems:
```Matlab
prob = vdx.Problem();
lbx = 0; ubx = 10;
lby = -10; ubx = 10;
for ii=1:n_masses
	x = SX.sym(['x_' num2str(ii)], 1);
	y = SX.sym(['y_' num2str(ii)], 1);
	init_x = ii*10/n_masses;
	init_y = -10 + ii*20/n_masses;
	prob.w.x(ii) = {x, lbx, ubx, init_x};
	prob.w.y(ii) = {y, lby, uby, init_y};
end
```
It is important to note that the variable classes do not need to be pre-defined as the first time a variable class is indexed into it is automatically created.
The syntax for adding a variable is as follows:
```
<variable>.<class>(<index>) = {<symbolic>[,<lower bound>, <upper bound>, <initial value>]}
<variable>.<class>(<index>) = {{<var_name>, <var_length>}[,<lower bound>, <upper bound>, <initial value>]}
```
the first of which takes a casADi symbolic vector as its first member and the second generates the casADi symbolic from a name and a vector size.
You can modify the bounds and initial values of a variable by assigning to the `lb`, `ub`, and `init` fields of a variable class:
```Matlab
prob.w.x(1).lb = 5;
prob.w.y(1).ub = 0;
prob.w.x(1).init = 7;
```
### Creating Constraints
The syntax for adding a constraint is identical to the syntax for adding vectors, however the second syntax is not useful because it creates a new variable.

### Solving VDX Problems and Extracting Results
Currently you can generate a casADi `nlpsol` by calling `Problem.create_solver(casadi_options)` where the argument is a struct containing the casADi nlpsol options.
Here is an example:
```Matlab
prob = vdx.Problem();
...
% define problem
...
casadi_options = struct;
...
% define casadi options
...
prob.create_solver(casadi_options);
stats = prob.solve();
```
Once the solver is created and solve is called the `mult` and `res` fields of `w`, `g`, and `p` are populated with the multipliers and results of the last solver call.
These values can be accessed as in the following example:
```Matlab
x_result = prob.w.x(ii).res;
y_result = prob.w.y(ii).res;
lambda_x = prob.w.x(ii).mult;
lambda_y = prob.w.y(ii).mult;
```
