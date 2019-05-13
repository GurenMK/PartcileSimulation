# PartcileSimulation

Serial Optimization:
Using vector to store all the bins and appropriate particles are stored in each bin which is also a vector. Vector of bins is padded to make apply forces loops generic and without having to check for edge bins. The downside is indices may be hard to follow when storing and retrieving particles from bins.

Bins are created and particles are assigned to bins at the beginning of every time step and the bins are destroyed at the end of every time step. I’s not as efficient as adjusting particle pointers but it reduces overhead of keeping track of all the particles and makes it easier to write correct code.

To calculate forces on the particles each particle is checked in order. The bin of each particle is calculated on the fly right before calculating apply forces function and counts as the current bin. For apply force to work another loop was added to only check appropriate neighboring bins. When thinking of the vector in 2d representation the top leftmost neighboring bin is found first in regard to the current bin. Then getting the neighboring bins starting with the first neighbor can be done consecutively and in order from the smaller bin index to the larger. A total of 9 bins are checked for each particle, its own bin and its neighboring bins.

The best performance gains were achieved by switching to vectors and the next jump in performance was by drastically increasing the number of bins. The number of bins is dynamic and is dependent of the total number of particles. First, based on the plane size the size of each bin is calculated which is plane size divided by the square root of the number of particles and divided by 11. The number 11 was chosen arbitrarily and seems to provide the best overall performance.  

OpenMP optimization:
Bins were changed to a shared data structure which fixed the correction issue. At the beginning of every time step a single thread assigns particles to appropriate bins and at the end a signle thread empties the bins from all particles. Apply forces loop was changed to reduce redundancy when calling apply forces. The algorithm still functions the same as in serial.
