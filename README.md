This code takes DICEBOX output and creates GEANT4 input with angular correlations. Needs ROOT (root.cern.ch).

The DICEBOX output has to be in fixed form and contain energies, spins, internal conversion flags, and mixing ratios. This implies the switches in DICEBOX input have to be like the following
* Switches (fixing the regime of the run) (ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS)
  1   0   1   1   0   1   1   0   0  !NOTE this needs to be exactly like this for the angle_generator code
DICEBOX will by default create files named EVENTS.S001.R001, which contain cascades, each spanning 5 lines and holding all necessary information, and can be used as inputs. 

The output file for GEANT4 (.out) then holds one cascade per line, starting with the number of gammas, N, then N energies of these gammas, and then (N-1) angles between these gammas.
The direction of the first gamma is advised to be random.
An arbitrary angle value of -1000 is given for uncorrelated gammas and it is left to the user's choice how to deal with these in their GEANT4 simulations.

A stub of the GEANT4 code, which processes the output file, will be provided.
