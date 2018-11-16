# Modified SelfCentering Material for OpenSees

This repository contains a modified implementation of the OpenSees `SelfCentering`
uniaxial material, developed for the simulation of Resilient Slip Friction Joint
(RSFJ) based self-centring damping systems.

The development was carried out as part of a collaborative research project with
researchers from The University of Auckland, New Zealand, focused on the numerical
modelling and seismic performance assessment of low-damage structural systems
equipped with RSFJ technology.

## About the material

The original OpenSees `SelfCentering` material was extended to represent the
multi-stage force–deformation response observed in RSFJ devices. Additional
parameters and constitutive rules were introduced to improve the simulation of
activation mechanisms, stiffness transitions, energy dissipation, and self-centring
behaviour under cyclic loading.

The formulation supports the modelling of resilient structural components and
systems incorporating RSFJ devices, including self-centring dampers and tension-only
bracing systems. The implemented model enables numerical simulation of a stable
flag-shaped hysteretic response with minimal residual deformation.

Validation was performed using experimental data from prototype RSFJ devices and
full-scale structural testing programs. The associated research included the
development and assessment of a self-centring tension-only brace system utilising
RSFJ technology, where experimental investigations demonstrated repeatable
flag-shaped hysteresis, reliable energy dissipation, and self-centring performance
under both quasi-static and dynamic loading protocols.

The implementation is intended for researchers and engineers interested in the
numerical modelling of damage-avoidance and resilience-based seismic design
solutions within the OpenSees framework.

## Usage

```tcl
uniaxialMaterial SelfCentering $tag $k1 $k2 $k3 $SecondFuse $epsA $epsB $ForcB \
    $a3 $a4 $a5 $a6 $a7 $a8 $ActF $beta <$SlipDef $BearDef $rBear>
```