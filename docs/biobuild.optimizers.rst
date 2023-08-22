The optimizers package
======================

These are the optimizers of Biobuild. Whenever you want to improve a molecule's conformation or generate new conformers for a molecule, you will want to use the optimizers.



.. .. automodule:: biobuild.optimizers
..    :members:
..    :undoc-members:
..    :show-inheritance:


The Rotatron Environments
-------------------------

Biobuild implements a torsional optimization system. That is, instead of "wiggling" atoms around until a genetically favorable structure is obtained, Biobuild rotates around 
bonds within a structure to find the most favorable conformation. This is done by using the Rotatron Environments. The Rotatron Environments are OpenAI Gym environments that
store an evaluation function to simulate rotating a molecule around a given set of bonds by a given set of angles. There is a base "Rotatron" environment and three subclasses thereof
that can be used for optimization heads-on. These are:


.. tab-set::

   .. tab-item:: The DistanceRotatron

      .. automodule:: biobuild.optimizers.DistanceRotatron
         :members:
         :undoc-members:
         :show-inheritance:

   .. tab-item:: The OverlapRotatron

      .. automodule:: biobuild.optimizers.OverlapRotatron
         :members:
         :undoc-members:
         :show-inheritance:
   
   .. tab-item:: The ForceFieldRotatron

      .. automodule:: biobuild.optimizers.ForceFieldRotatron
         :members:
         :undoc-members:
         :show-inheritance:

   .. tab-item:: The Base Rotatron
         
         .. automodule:: biobuild.optimizers.Rotatron
            :members:
            :undoc-members:
            :show-inheritance:
            

Optimization algorithms
-----------------------

The `Rotatron` environments are used to specify the problems to solve. The optimization algorithms are used to solve them. Biobuild implements a number of classical optimization algorithms
that are tailored to work with the Rotatron environments. These are:

.. dropdown:: Particle Swarm Optimization

   The Particle Swarm Optimization algorithm is a classical optimization algorithm that is based on the behavior of a swarm of particles. Each particle has a position and a velocity. The position
   is the current solution to the problem, and the velocity is the direction in which the particle is moving. The particles are attracted to the best solution found so far, and repelled by the
   worst solution found so far. This way, the particles will move towards the best solution found so far, and will not get stuck in local minima. 

   The algorithm performs well with both small and large inputs, both with AtomGraphs and ResidueGraphs. It is also often the fastest to compute, so it is the `default algorithm` .

   .. autofunction:: biobuild.optimizers.algorithms.swarm_optimize

.. dropdown:: Genetic Algorithm

   The Genetic Algorithm is one of the most iconic optimization algorithms. It is based on the behavior of a population of individuals. Each individual has a "genome", which is the current solution to the problem.
   Each generation (optimization round) individuals are mutated (randomly change their solution), and the best individuals reproduce and make it to the next round. This way, the population will move
   gradually towards good solutions. 

   The algorithm performs well on any scale but gets exceedingly slower the larger the molecules become. Also, it works slightly better with AtomGraphs than with ResidueGraphs.

   .. autofunction:: biobuild.optimizers.algorithms.genetic_optimize

.. dropdown:: Simulated Annealing

   Simulated Annealing is another optimization algorithm that has similarities to both genetic and particle swarm optimization. It explores solutions by randomly changing the current one,
   and accepts or rejects the new solution based on the change in energy. The algorithm is based on the annealing process in metallurgy, where a metal is heated and then slowly cooled down.
   This way, the metal will settle in a more stable state. 

   The algorithm performs well with better with smaller inputs but is suitable for larger ones using both AtomGraphs and ResidueGraphs.

   .. autofunction:: biobuild.optimizers.algorithms.anneal_optimize

.. dropdown:: Gradient-based algorithms

   We implement a direct link to `scipy.optimize.minimize` which provides a number of gradient-based optimization algorithms. These algorithms are usually very fast and perform well on small inputs. However, as evaluation landscapes
   of larger molecules tend to get "rugged" gradient-based methods tend to struggle with larger inputs.

   Any algorithm implemented by `scipy.optimize.minimize` can be used. The default is the `L-BFGS-B` algorithm. For a complete list of available algorithms checkout the `scipy documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`_.

   .. autofunction:: biobuild.optimizers.algorithms.scipy_optimize


.. .. automodule:: biobuild.optimizers.algorithms
..    :members:
..    :undoc-members:
..    :show-inheritance:

Optimization utilities
----------------------

Biobuild also implements a number of utilities that can be used to make the optimization a little easier for the user by automizing certain steps.

.. automodule:: biobuild.optimizers.utils
   :members:
   :undoc-members:
   :show-inheritance:

.. Module contents
.. ---------------

.. .. automodule:: biobuild.optimizers
..    :members:
..    :undoc-members:
..    :show-inheritance:
