# Group velocity extraction (`GVel`)

We need the group velocity to convert from a transition rate to a cross section. The group velocity is not straightforward to calculate across bands because of possible degeneracies and/or band crossings. However, while not complicated, this task is very tedious and long to do by hand. The purpose of this code is to automate most, if not all of the group velocity calculations. 

Each band at each k-point will have its own group velocity. The group velocity has 3 components in the x/y/z directions. To find group velocity you must:
* Generate 6 slightly-displaced k-points with zero weight around each original k-point (positive and negative in each direction).
* At each k-point and for each direction, match up the bands that fall in a line, taking account of any potential band crossings and/or degeneracies.
* Find the slope of each line, then convert the slope to a group velocity component using $$v\_i = \frac{1}{\hbar} \frac{\partial E}{\partial k\_i},$$ where $i$ indicates each direction x,y,z.
* Combine the components using $$||\mathbf{v}|| = \sqrt{v\_x^2 + v\_y^2 + v\_z^2}$$

This process is straightforward when there are no degeneracies and/or band crossings, but the opposite case is not so easy. There are main problems to solve: matching up bands that fall in a line and combining velocity components (especially with degeneracies). Some things may need to be done by hand, but the goal of this code is to create a tool to assist in that process and minimize the amount of manual work.
