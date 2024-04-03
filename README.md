# VDX: Automatic Index Tracking for CasADi NLPs
VDX is a Matlab wrapper for the creation of generic CasADi NLPs which provides automatic, transparent, index tracking facilities for named variable classes.
It also provides a direct interface to CasADi nlpsol and allows for automatic results extraction for the mentioned variable classes.
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
lby = -10; uby = 10;
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
the first of which takes a CasADi symbolic vector as its first member and the second generates the CasADi symbolic from a name and a vector size.

If `<index>` is not scalar one can initalize multiple indices at once. THis
This is useful to compactly populate the variable in multiple dimensions.
You can modify the bounds and initial values of a variable by assigning to the `lb`, `ub`, and `init` fields of a variable class:
```Matlab
prob.w.x(1).lb = 5;
prob.w.y(1).ub = 0;
prob.w.x(1).init = 7;
```

### Variable Groups
VDX provides the facility to group variables and index into this group.
```Matlab
var_group = vdx.VariableGroup({<members>}[, {<indexing rules>}])
```
where the first argument is a list of `vdx.Variable` and the optional second parameter allows for a custom indexing function with the signature: `(index,vdx.Variable) -> index`.
Simple indexing rules are implemented in `vdx.indexing` and include rules that map indicies to the (lexicographically) previous or next index.
If variables of different dimension are used in the same variable group then the lower dimensioned variables are indexed by a truncation of the index.
For example the group with `u` 10 by 1, and `x` 10 by 4, indexed at `(4,4)` would yield the vector: `[u(4), x(4,4)]`.

### Creating Constraints
The syntax for adding a constraint is identical to the syntax for adding vectors, however the second syntax is not useful because it creates a new variable.

Constraints for multiple indices can be created via the special syntax:

```Matlab
<variable>.<class>(<index>) = {{<function:casadi.Function>, {<arguments>}[, {<indexing rules>}]}[,<lower bound>, <upper bound>, <initial value>]}
```
with the arguments and indexing rules forming the variable group which makes up the arguments to the function.
This is then applied for all the indices.

### Solving VDX Problems and Extracting Results
Currently you can generate a CasADi `nlpsol` by calling `Problem.create_solver(casadi_options)` where the argument is a struct containing the CasADi nlpsol options.
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
full_result = prob.w.res;

x_result = prob.w.x.res;
y_i_result = prob.w.y(ii).res;
lambda_x = prob.w.x.mult;
lambda_y_ii = prob.w.y(ii).mult;
```

### Other Problem types
We also provide a package of specialized problem formulations.
Currently the only supported formulation is for Mathematical Programs with Complementarity Constraints (MPCC):
```
minimize f(w,p), subject to
lbw <= w <= ubw
lbg <= g(w,p) <= ubg
0 <= G(w,p) perp H(w,p) >= 0
```

## Implementation Details and Notes
There are some details that are important to note:
* Indexing into variable class is zero indexed.
* Assigning to `init`, `lb`, and `ub`, can only be done for scalar indexes. Updating VDX to allow for assigning to non-scalar indices is a soon to be implemented.
* We assume that the length of each index of a variable class is uniformly a single value or empty.
* When accessing a variable class numeric value (`lb`, `ub`, `init`, `res`, `mult`), the values are horizontally concatenated and in column major order, i.e. in "lexicographical" order.
* This project is similar in its goals to the existing (but as far as I know deprecated) `structures` concept in the `casadi.tools` python package. As of now we do not plan to guarantee any sort of common interface however.
