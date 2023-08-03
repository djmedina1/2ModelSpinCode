
<p> The goal of this project is to simulate a system with two spins under an applied field. We are interested in plotting a hysteresis loop. The mechanics of this system can be better understood through the reading provided in the respective subsection. A simple GUI is implemented to allow users to vary parameters of the simulation. </p>

# Single Spin System (Uniaxial Case)
<p> The process of plotting the hysteresis loop of a single spin system is as follows:</p>
<ul> 
  <li>Define an anisotropic axis angle.</li>
  <li>Define a range of field strengths H to apply.</li>
  <li>Locate the energy where the angle in which the energy is at a minimum for an applied field strength of 0T. This is the uniaxial angle - If there are two angles that corrospond to a minimum take the first angle found. </li>
  <li> Calculate the energy value of the system at an applied non-zero field strength, and track the evolution of this minimum with respect to the uniaxial angle minimum. If an inflection point is encontured before a minimum is found, include this in the tracking.  </li>
  <li>Repeating this process for increasing H fields and decreasing H fields will yield two different curves. </li>
  <li> Plot the cosine of the minimum field angles with respect to the applied field strength to produce a hysteresis loop.</li>
</ul>

 <a href ="public/Simulation_Diagram.png"></a>

=======
## This repo only contains a working simulation of a **one** spin system. 2 spins is currently being developed (as of 6/30/23). 
