porosity
$$
\varepsilon=\frac{V_{\text {void }}}{V_{\text {cell }}}
$$
Notion of saturation for phase $i$: 
$$
S_i = \frac{V_i}{V_{void}}
$$
Total saturation for a two-phase flow in porous media (a non-wetting phase $a$ and wetting phase $b$):
$$
1 = S_a + S_b
$$
Considering an incompressible multiphase flow in a porous medium, the macro-scale mass balance equation for each phase $i$ reads:
$$
\varepsilon \frac{\partial S_{i}}{\partial t}+\nabla \cdot \mathbf{U}_{i}=q_{i}
$$
where $ \mathbf{U}_{i} $ stands for the superficial velocity and $ q_{i}$ is a source term, used for injection or extraction wells.

The generalized Darcy's model [13] with superficial velocity of each phase $ i $:
$$
\mathbf{U}_{i}=-\frac{\mathrm{K}_{i}}{\mu_{i}} \cdot\left(\nabla p_{i}-\rho_{i} \mathbf{g}\right),
$$

- $  \mathrm{K}_{i}$ --  apparent permeability
  $$
  \mathrm{K}_{i}=\mathrm{K} k_{r i}\left(S_{b}\right)
  $$

- $ \mathrm{K}$ -- permeability tensor of the porous medium

- $S_{b}$ -- the local saturation of the wetting phase

- $k_{r i}\left(S_{b}\right)$  is the relative permeability of the phase $ i$ depending on the local saturation of the wetting phase





## `pEqn.H`

To simplify the formulation, we define phase mobility $  M_{i} $ and gravitational contribution $ L_{i} $ as follows:
$$
M_{i}=\frac{K k_{r i}\left(S_{w}\right)}{\mu_{i}}, \quad \text { phase mobility }
$$

$$
L_{i}=\frac{K k_{r i}\left(S_{w}\right)}{\mu_{i}} \rho_{i} . \quad \text { gravitational contribution }
$$

Even if in the generalized Darcy's law the relation $ L_{i}=\rho_{i} M_{i} $ is satisfied, we have found it convenient to separate each contribution. 



```c++
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc)*activateCapillarity
            ==
            // event source terms
            - sourceTerm
        );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc*activateCapillarity;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    phib = Fbf*phi+(((Lbf-Fbf*Lf)&g) & mesh.Sf())+(1-Fbf)*phiPc*activateCapillarity;
    phia = phi - phib;

    Ub = fvc::reconstruct(phib);
    Ua = U-Ub;

    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();  

}
```

