---
layout: post
title: "Understanding viscoelasticity: Why are some materials stretchy but slow to snap back?"
description: What is viscoelastic material? Why do we study these materials and how do we design products using them?
categories:
  - coding
  - presentation
tags:
---

Imagine you're playing with two different types of materials: a spring (representing an elastic material) and a piece of silly putty (representing a viscoelastic material). When you pull and then release the spring, it snaps back immediately to its original shape. This instant response to pulling and releasing is characteristic of **elastic materials** – they return to their original shape as soon as the force is removed, just like a rubber band.

Now, think about what happens when you stretch or squish the silly putty. It stretches or deforms slowly under your pull and retains its deformed shape for a while before slowly returning to its original shape after you let go. This behavior – the slow, time-dependent return to the original shape – is a hallmark of **viscoelastic materials**. They have a memory of sorts, taking time to go back to how they were, combining both the "instant snap" of elastic materials and the "slow flow" of viscous materials (like honey).

The difference between elastic and viscoelastic materials fundamentally lies in how they respond over time. Elastic materials don't care about time; they react the same way whether you stretch them quickly or slowly. Viscoelastic materials, on the other hand, care a lot about how fast or slow things happen. Stretch them quickly, and they might resist more strongly; stretch them slowly, and they might gradually give way.

This time-dependent behavior of viscoelastic materials is crucial in their design and application. For example:

- **Shock Absorption**: In products like car bumpers or sports helmets, viscoelastic materials are used because they can absorb and dissipate energy over time, providing better protection against impacts.
- **Comfort Products**: Memory foam mattresses and pillows use viscoelastic materials to conform to the shape of your body over time, offering personalized support and comfort.
- **Medical Devices**: Viscoelastic materials are used in prosthetics and implants because they can mimic the time-dependent behavior of human tissues, making the devices more comfortable and functional for users.
- **Vibration Damping**: In buildings and machinery, viscoelastic materials help damp vibrations, not just instantly but over a period, reducing noise and wear and tear.

Designers and engineers leverage the unique "time-aware" properties of viscoelastic materials to create products that perform better under specific conditions, providing benefits that purely elastic materials can't match. They carefully consider how these materials will be stretched, compressed, or twisted over time to ensure the final product behaves exactly as needed, whether that's in ensuring comfort, safety, durability, or any other desired property.

Understanding viscoelastic materials and their behavior involves diving into the world of material science and mechanics, where the properties of these materials are described and predicted through mathematical equations. Let's look at the foundational equations that help explain why viscoelastic materials stretch and slowly return to their original shape.

### 1. Stress-Strain Relationship

In viscoelastic materials, the relationship between stress (force per unit area) and strain (deformation) is time-dependent. Two primary models often used to describe this behavior are the Maxwell model and the Kelvin-Voigt model.

#### Maxwell Model

The Maxwell model represents viscoelastic behavior as a combination of a purely elastic spring and a purely viscous dashpot in series. The constitutive equation is:

$$
\sigma + \frac{\eta}{E} \frac{d\sigma}{dt} = \eta \frac{d\varepsilon}{dt}
$$

where:
- $\sigma$ is the stress,
- $\varepsilon$ is the strain,
- $E$ is the elastic modulus of the spring,
- $\eta$ is the viscosity of the dashpot,
- $t$ is time.

This model captures fluid-like behavior, showing that under a constant stress, the material will continue to deform over time.

#### Kelvin-Voigt Model

The Kelvin-Voigt model describes viscoelastic behavior using a spring and dashpot in parallel. Its constitutive equation is:

$$
\sigma = E\varepsilon + \eta \frac{d\varepsilon}{dt}
$$

This model captures solid-like behavior, where the material deforms immediately under stress and slowly returns to its original shape when the stress is removed, but it does not flow over time as the Maxwell model predicts.

### 2. Creep and Recovery

Creep describes how a material deforms under a constant load over time, while recovery describes how a material returns to its original shape after the load is removed. The creep function $J(t)$ for a viscoelastic material is often represented as:

$$
J(t) = \frac{\varepsilon(t)}{\sigma_0}
$$

where $\sigma_0$ is the constant stress applied. The recovery process, especially in nonlinear viscoelastic materials, might not follow the same function due to the material's history-dependent behavior.

### 3. Complex Modulus

For viscoelastic materials subjected to oscillatory loading, the complex modulus $E^*$ is used to describe the material's response, combining storage modulus $E'$ (elastic response) and loss modulus $E''$ (viscous response):

$$
E^* = E' + iE''
$$

where $i$ is the imaginary unit. The phase angle $\delta$ between stress and strain gives insight into the material's viscoelastic nature:

$$
\tan(\delta) = \frac{E''}{E'}
$$

High values of $\tan(\delta)$ indicate more viscous behavior, while lower values indicate more elastic behavior.

These equations and models provide a mathematical framework to understand and predict how viscoelastic materials behave under various loading conditions. By analyzing these time-dependent stress-strain relationships, scientists and engineers can design materials and products that leverage the unique properties of viscoelasticity for applications requiring both elasticity and viscosity, such as impact absorption, comfort, and durability.

## Variational form

In the context of viscoelastic materials, in addition to the elastic energy stored in the material, there's an energy dissipation component due to the material's viscous behavior. This dissipation energy represents the loss of mechanical energy to heat or internal friction as the material deforms. Unlike purely elastic materials, where the energy is fully recoverable, viscoelastic materials convert a portion of the energy into heat, leading to a hysteresis effect in stress-strain cycles.

### Dissipation Energy in Viscoelasticity

Dissipation energy in viscoelastic materials is a crucial aspect that differentiates them from purely elastic materials. It quantifies the energy lost during cyclic loading due to the material's inherent viscosity. The rate of energy dissipation is related to the material's relaxation processes and is a function of both the applied strain rate and the material's viscoelastic properties.

The dissipated energy $\Delta W_{diss}$ over a loading cycle can be expressed as the area enclosed by the hysteresis loop in a stress-strain diagram, which can be calculated for a cycle from $t_1$ to $t_2$ as:

$$
\Delta W_{diss} = \int_{t_1}^{t_2} \sigma(t) \frac{d\varepsilon(t)}{dt} dt
$$

where $\sigma(t)$ is the stress and $\varepsilon(t)$ is the strain at the time $t$.

### Weak Form for Viscoelasticity

The weak form for viscoelastic materials extends the concept of the weak form in linear elasticity by incorporating the time-dependent behavior. In linear elasticity, the weak form is derived by taking the first variation of the total potential energy. For viscoelasticity, we need to consider both the stored elastic energy and the dissipative energy.

For a viscoelastic material undergoing small deformations, the weak form of the equilibrium equation, incorporating the dissipative effects, can be formulated as follows:

Given a body $\Omega$ with boundary $\partial\Omega$, the weak form states that for all admissible variations $\delta\varepsilon$ of the strain field, the following condition holds:

$$
\int_{\Omega} \left( \sigma_{ij}(\varepsilon) \delta\varepsilon_{ij} - \rho f_i \delta u_i \right) d\Omega + \int_{\Omega} \eta \frac{d\varepsilon_{ij}}{dt} \delta\varepsilon_{ij} d\Omega = 0
$$

where:
- $\sigma_{ij}(\varepsilon)$ is the stress tensor, which is a function of the strain tensor $\varepsilon_{ij}$,
- $\rho$ is the density,
- $f_i$ are body forces,
- $u_i$ are displacements,
- $\eta$ represents the viscosity coefficient, and
- $\frac{d\varepsilon_{ij}}{dt}$ is the strain rate.

The first integral represents the internal and external work done on the system, similar to the elastic case. The second integral accounts for the energy dissipation due to the viscous nature of the material, with $\eta \frac{d\varepsilon_{ij}}{dt}$ representing the viscous stress component.

This formulation captures the essence of viscoelastic behavior by including both the instantaneous elastic response and the time-dependent viscous response, allowing for the analysis of viscoelastic materials under various loading conditions. The weak form is essential for implementing finite element methods for viscoelastic problems, where the time-dependent nature of the material's response must be accurately captured.