# Effect of Accelerated Global Expansion on the Bending of Light

In this research we consider two types of observers, Static and Comoving which in the limit of constant acceleration (Hubble Constant `H`) correspond to the same definition that is comonly used. McVittie is a Non Vacuume solution to the Einstein equations which corresponds to a spacetime that is rate of expansion accelerates and has a black hole at the center of the coordinates. In that sense it is a good representative of our universe which is known to have an accelerated rate of expansion, however, this solution requires the presense of some type of fluid in the space. In many cases the above problem is worrysome, but for the purpose of light bending it does not cause any problems because it does not interact with the light ray.

We developed a Matlab code that shoots the light ray from the preassumed position of the source star and then follow the path of the light until it reaches to the otherside, where we assumed that the source, lens (black hole) and the observer are all aligned. If the light ray does not reach the location that is required, it iteratively changes the shooting angle until it reaches to the required location with some small tollorance `0.005`.

![Bending](https://github.com/maghili/McVittie/blob/master/Bending.png)

Since the observation is done in a curved spacetime the regular `Euclidean` angles do not work, but in any case the change in this angle is shown here

![Euclidean](https://github.com/maghili/McVittie/blob/master/Euclidean.png)

After making required transformations we can find the change of the angle based on the values of the constant rate `H` and the acceleration `A`.

![Static](https://github.com/maghili/McVittie/blob/master/Static.png)

Comoving is the closest choice to an observer sitting on the Earth and being dragged by the expansion

![Comoving](https://github.com/maghili/McVittie/blob/master/Comoving.png)
