Reduced Order Modelling is a technique used to make the solution of large physical simulations cheaper, by using data from previous simulations and even from experiments in order to express a dynamical system in very efficient coordinates.

In traditional methods, like Finite Volumes, Finite Differences and Finite Elements, we increase the dimensionality of the problem in order to turn the low dimensional PDE into a high dimensional ODE, where each direction in our new coordinate system represents a degree of freedom in a specific node of our spatial discretization.

Since we don't use all the degrees of freedom, meaning that not every state of the state vector could be a valid state (like a random one), we could more efficiently express our problem in other coordinates, which could be Fourier coordinates, Wavelet coordinates, or even tailored coordinates, which use data from other experiments to build a coordinate basis using the Singular Value Decomposition (SVD).

Suppose our problem can be represented in the form of a non-linear dynamical system

![eq1](https://latex.codecogs.com/png.image?\dpi{130}&space;\bg_white&space;\mathbf{u}_t(\mathbf{x},t)&space;=&space;\mathbf{N}(\mathbf{u},&space;\mathbf{u}_x,&space;\mathbf{u}_y,&space;\mathbf{u}_z,&space;\mathbf{u}_{xx},&space;\mathbf{u}_{yy},&space;\dots,&space;\mathbf{x},&space;t))

If we use the separation of variables technique

![eq2](https://latex.codecogs.com/png.image?\dpi{120}&space;\bg_white&space;\mathbf{u}(\mathbf{x},t)&space;=&space;\sum\limits_{k=1}^{r}&space;\psi_k(\mathbf{x})&space;\mathbf{a}_k(t))

And our basis ![eq3](https://latex.codecogs.com/png.image?\dpi{120}&space;\bg_white&space;\Psi&space;=&space;\begin{bmatrix}\psi_1&space;&&space;\psi_2&space;&&space;\dots&space;&&space;\psi_r\end{bmatrix}) is orthonormal, meaning that

![eq4](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cleft%3C%5Cpsi_i,%20%5Cpsi_j%5Cright%3E%20=%20%5Cdelta_%7Bij%7D%20=%20%5Cbegin%7Bcases%7D1%20%5Ctext%7B%20if%20%7D%20i=j%20%5C%5C%200%20%5Ctext%7B%20if%20%7D%20i%5Cneq%20j%5Cend%7Bcases%7D)

Then we can substitute the separation of variables into the dynamical system equation

![eq5](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_k(%5Cmathbf%7Bx%7D)%20%5Cmathbf%7Ba%7D_%7Bk,t%7D(t)%20=%20%5Cmathbf%7BN%7D%5Cleft(%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_k(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_%7Bk,x%7D(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_%7Bk,y%7D(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Cdots,%20%5Cmathbf%7Bx%7D,%20t%5Cright))

And take the internal product with ![eq6](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cpsi_i(%5Cmathbf%7Bx%7D)) in both sides

![eq7](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbf%7Ba%7D_%7Bi,t%7D(t)%20=%20%5Cleft%3C%5Cpsi_i(%5Cmathbf%7Bx%7D),%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_k(%5Cmathbf%7Bx%7D)%20%5Cmathbf%7Ba%7D_%7Bk,t%7D(t)%20%5Cright%3E%20=%20%5Cleft%3C%20%5Cpsi_i(%5Cmathbf%7Bx%7D),%20%5Cmathbf%7BN%7D%5Cleft(%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_k(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_%7Bk,x%7D(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Csum%5Climits_%7Bk=1%7D%5E%7Br%7D%5Cpsi_%7Bk,y%7D(%5Cmathbf%7Bx%7D)%5Cmathbf%7Ba%7D_%7Bk%7D(t),%20%20%5Cdots,%20%5Cmathbf%7Bx%7D,%20t%5Cright)%20%5Cright%3E)

If we discretize the whole problem in space, expressing each basis function as a n-tall vector, where n is the number of discretization points, and ![eq8](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5CPsi) being the matrix whose columns are the basis functions, then the problem can also be expressed as

![eq9](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbf%7Ba%7D_t(t)%20=%20%5CPsi%5ET%20%5Cmathbf%7BN%7D(%5CPsi%5Cmathbf%7Ba%7D,%20%5CPsi_x%5Cmathbf%7Ba%7D,%20%5CPsi_y%5Cmathbf%7Ba%7D,%20%5CPsi_z%5Cmathbf%7Ba%7D,%20%5CPsi_%7Bxx%7D%5Cmathbf%7Ba%7D,%20%5CPsi_%7Byy%7D%5Cmathbf%7Ba%7D,%20%5Cdots,%20%5Cmathbf%7Bx%7D,%20t))

Where ![eq10](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbf%7Ba%7D) is a ![eq11](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20n%5Ctimes%20d) dimensional vector, being d the dimension of ![eq12](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbf%7Bu%7D(%5Cmathbf%7Bx%7D,t)) and ![eq13](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbf%7Ba%7D_k(t)).

If there is as many basis functions as there is discretization points, then the basis functions span ![eq14](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Cmathbb%7BR%7D%5En) and there is no advantage in implementing the reduced order modeling. However, if we choose the basis functions such that the first few basis functions can efficiently represent the behaviour of our problem, then we can use way less basis functions than discretization points, and therefore our solution becomes much less costly than a traditional one.

A good way of choosing ![eq15](https://latex.codecogs.com/png.image?%5Cdpi%7B120%7D%20%5Cbg_white%20%5CPsi) is using the Singular Value Decompostion, an operation that extracts the modes with the bigger variance from an already solved problem, or even from an actual experiment.

In this repository, I'll start solving very easy problems like the Heat Equation problem which is nice and linear, and the goal is to solve more complex and non-linear problems.

All the procedures were inspired by the YouTube channels of the professors Steve Brunton and Nathan Kutz, as well as their book Data Driven Science and Engineering.
