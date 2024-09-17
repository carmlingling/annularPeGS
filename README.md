# diskSolve
This module uses a nonlinear least-squares solver to match the photoelastic fringe pattern of each *experimental* particle to a *theoretical* fringe pattern. The best fit theoretical fringe pattern provides the magnitude and direction of each force acting on a particle.

The input of this module is the standardised MATLAB **particle** structure array. Each particle in the image has the following fields:

- `id` : assigned ID number of the particle
- `x`: x coordinate of the particle centre, in *pixels*
- `y` : y coordinate of the particle centre, in *pixels*
- `r` : radius of the particle, in *pixels*
- `rm` : radius of the particle, in *metres*
- `color` : assigned particle color (for plotting)
- `fsigma` : photoelastic stress coefficient of the particle
- `z` : number of contacts on the particle
- `f` : average force on the particle calculated using the G<sup>2</sup> method
- `g2` : G<sup>2</sup> value of the particle image
- `forces` : array of contact force magnitudes 
- `betas` : array of contact force azimuthal angles
- `alphas` : array of contact force 'contact angles'. 0 is a purely normal force, $\pm \pi/2$ is a purely tangential force
- `neighbours` : array containing particle IDs of contacting particles
- `contactG2s` : array of G<sup>2</sup> values, corresponding to area around each contact force
- `forceImage` : cropped image of particle from original experimental image

diskSolve requires all fields of **particle** to be populated *other than* `forces` and `alphas`. diskSolve also adds the following fields:

- `fitError` : error of the least-squares fit
- `synthImg` : fitted theoretical fringe pattern

In addition to this, the following user inputs are requried:

- `directory` : file location of **particle**
- `algorithm` : algorithm used for the least-squares fit. Default set to `levenberg-marquardt`
- `maxIterations` : maximum number of iterations for solver. Default set to 200
- `functionTolerance` : error tolerance for solver. Default set to 0.01
- `scaling` : scale factor to change size of experimental image for solver input. Default set to 0.5
- `maskradius` : percentage of particle radius to use for fit. i.e. 1 - `maskradius` is fraction of radius that is discarded

The only necessary input for the module to run is `directory`. The other inputs can be tuned to improve accuracy and/or run time.

---

## Usage
To run this module, input the `directory` variable on line 6 of `diskSolve.m` and run the script. This will call the following functions:

- `solver_original.m` : calculates initial conditions for least-squares fit and performs fit
- `forceBalance.m` : applies force balance to the particle. **UNMODIFIED FROM PEGSV1 VERSION**
- `fringe_pattern_original.m` : creates fringe pattern for given force magnitudes and angles
- `stress_engine_original.m` : principal stress calculation to provide pixel intensities for fringe pattern 
 
--- 
### PLEASE READ
This version functions the same as PeGS Version 1.0. It includes C Lee's parallelisation modification and B McMillan's fix to the stress tensor calculation. The following features will be added soon:

- Vectorised version of the fringe pattern calculation to improve efficiency
- Line contact solution for forces. Currently all forces are assumed to be point contacts
- The ability to switch between original version, vectorised version, line contact and point contact solutions
- The ability to run diskSolve on multiple video frames
- The ability to turn on/off force balance requirement

The ability to watch the fringe pattern converge on screen has been removed due to the parallelisation requirements. Comparing the experimental fringes to theoretical ones should now be done in post-processing.

**I have not validated or modified the force balance calculations** 
