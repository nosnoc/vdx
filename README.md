# VDX: Automatic index tracking for casADi NLPs
VDX is a wrapper for the creation of generic casADi NLPs which provides automatic, transparent, index tracking facilities for named variable classes.
It also provides a direct interface to casADi nlpsol and allows for automatic results extraction for the mentioned varable classes.
## Quickstart guide for users
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

### Creating variables
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
<variable>.<class>(<index>) = {<symbolic>,[<lower bound>, <upper bound>, <initial value>]}
<variable>.<class>(<index>) = {{<var_name>, <var_length>},[<lower bound>, <upper bound>, <initial value>]}
```
