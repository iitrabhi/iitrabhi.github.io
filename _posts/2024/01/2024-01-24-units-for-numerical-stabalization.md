---
layout: post
title: Understanding the Use of Unconventional Units in Computational Ice Dynamics
description: How does using an unconventional unit system comprising megapascals, meters, and years in computational ice dynamics facilitate alignment with the enormous spatial and temporal scales of glaciological phenomena and enhance the numerical stability and precision of simulations?
categories:
  - simulation
tags:
  - glacier
  - firedrake
  - fenics
---
In the realm of computational ice dynamics, particularly in models like Elmer/ICE, a peculiar unit system is often employed: megapascals (MPa), meters (m), and years. This choice deviates from the more conventional metric system (MKS: meters, kilograms, seconds). This blog post delves into the rationale behind this choice, illustrating how it aids in dealing with the peculiarities of glaciological computations and examines whether this approach impacts numerical stability.

## The Problem with Conventional Units in Glaciology

### Scale of Glaciological Phenomena
Glaciological phenomena occur over vast spatial scales (kilometers) and long temporal scales (years to millennia). This already sets a challenge for conventional MKS units, where seconds as a time unit are too granular for practical computation.

### Physical Constants in MKS Units
A key issue arises with the values of physical constants in glaciological modeling. For instance, the Glen flow law rate factor, a crucial constant in ice flow modeling, exhibits extremely large or small values in MKS units. Such extreme values can introduce computational difficulties, particularly in numerical stability and precision.

## Advantages of Megapascals - Meters - Years Unit System

### Sensible Range of Physical Constants
By employing megapascals, meters, and years, physical constants like the Glen flow law rate factor fall within a more manageable numerical range. This choice reduces the risk of computational errors from handling very large or small numbers.

### Alignment with Physical Processes
The choice of years as a time unit aligns better with the temporal scale of ice dynamics. Similarly, using megapascals, a higher pressure unit, matches the scale of pressures typically encountered in glaciological studies.

## Numerical Stability and the Unit System

### Impact on Computations
Numerical stability in computational models is crucial to ensure accurate and reliable results. The choice of unit system can influence this stability. Extreme values (either too large or too small), common in MKS units for glaciological constants, can lead to numerical issues like underflow, overflow, or loss of precision.

### How Megapascals - Meters - Years Help
By normalizing the range of constant values, the megapascals - meters - years system mitigates some of these numerical issues. It helps balance the magnitude of numbers processed during simulations, contributing to improved numerical stability.

## Examples
Consider the Glen flow law, which describes the relationship between stress and strain rate in ice. In MKS units, the rate factor (A) can be extremely small (e.g., $10^{-24}  s^{-1}Pa^{-3}$), posing challenges in computational precision. In contrast, when using the $MPa-m-years$ system, this value is more manageable (e.g., $10^{-3} year^{-1}MPa^{-3}$), reducing the risk of numerical errors.

## Conclusion
Adopting megapascals - meters - years as a unit system in glaciological modeling is a pragmatic response to the unique challenges posed by the scale and nature of ice dynamics. It ensures that physical constants remain within a sensible numerical range and enhances the numerical stability of the simulations. This approach exemplifies how adapting conventional practices to the specific needs of a field can lead to more effective and accurate scientific modeling.
