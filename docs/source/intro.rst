Intro
=====

Welcome to the docs for The PPOLs Model - a Python-based planet formation model that explores the "pebble snow" mechanism. Improvements to the docs are on-going. Please reach out to S. McCloat via email for questions, code, or links.

The model and results have been peer-reviewed and accepted for publication in the Astrophysical Journal (ApJ), and will appear in an upcoming volume. You can see the accepted manuscript on arXiv [here].

Plain Language Summary
----------------------
We can look at the Solar System in terms of the overall arrangement of the planets: relatively smaller, rocky planets (like Venus, Earth, and Mars) are found closer to the Sun, and the larger gas giants (like Jupiter and Neptune) are found farther away. Up until the mid-2010s, this was the only "arrangement" for planetary systems we were aware of. In the wake of discovering thousands of planets in other systems, one mystery we face is how the arrangement of the Solar System fits into the bigger picture of planetary systems - are we unusual for stars like the Sun? Should planets like Earth have more water, or less? Where should we expect to find systems with a planet like Jupiter together with a planet like Earth? 

We approach this mystery by exploring how planets form, and specifically with a process called "pebble accretion" that has developed along with the discoveries of other planetary systems. In the classic picture, planets grow by smashing large asteroid-sized bodies together until there remain the last planets standing. With this newer idea, small bits of rock the size of M&Ms or golf balls rain out of the early material of a Solar System and fall towards the star in the middle. A planet that is just starting to grow can collect tons of these bits along the way and reach planet-sizes.

The PPOLs Model looks at this process over a wide range of stars (some bigger and smaller than the Sun), and considers where those pebbles are dry or contain water, and tracks how big the planets can grow and how much water they receive. We found that three arrangements of planetary systems were typical: a small collection of roughly Earth-size planets tucked in close to the star and no gas giants; systems of only gas giants located farther away from the star; and middling arrangement that does have both an Earth-size planet and a gas giant planet together. Most interesting of all - the conditions that led to the Earth + gas giant together are a relatively narrow set of conditions, which could explain why our Solar System (as of yet) appears to be uncommon.



Features of the Model
---------------------

Users can easily set and run planet formation models for a variety of disk models. Major features include:

* a disk model with explicit stellar mass and disk mass
   
   * disk mass can be fraction of stellar mass, or have the solid disk mass set exactly

* dust mass that depletes as it converts into pebbles
* temperature profile that is a simple power law or account for viscous/irradiation heating (e.g. Ida et al. 2016)
* a snow line (water-ice) based on temperature, set explicitly, or self-consistent with disk parameters
   
   * the snow line that can evole (or not) based on disk parameters

* any number of planetesimal seed masses, with any initial mass, introduced at any time
* planetesimal seed masses that grow simultaneously via pebble accretion and remove ("filter") pebble mass as it drifts inward
* accounting for dry vs water-rich solid mass that is accreted with time

The model was built to test outcomes of pebble accretion, specifically during first ~10 million years of the protoplanetary disk when gas is still present. It may serve well as inputs for longer n-body simulations.