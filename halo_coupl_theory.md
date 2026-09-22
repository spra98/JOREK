# Introduction and notation

In this note we show how to extend the Virtual Casing Principle ideas to
the case where there are shared currents between the inner and outer
domain. It is convenient to set out immediately the notation. Here the
equivalent currents to the plasma are represented via the degrees of
freedom $\mathbf{{I_{eq}}}$, while the external currents are represented
via the degrees of freedom $\mathbf{{I_w}}$. The basis functions for the
CARIDDI equivalent shell and wall include in this case those currents
which cross the Coupling Surface at specified electrodes, *e.g.* for the
equivalent shell: $$\begin{equation}
    X_{eq} = \{ \mathbf{w}_k \in L_2(\Omega_{eq}) : \nabla \cdot \mathbf{w}_k = 0 \} \, .
\end{equation}$$ The currents which are solenoidal within the equivalent
shell are contained in the subspace $$\begin{equation}
    X_{eq,0} = \{ \mathbf{w}_k \in L_2(\Omega_{eq}) : \nabla \cdot \mathbf{w}_k = 0 , \, \mathbf{w}_k \cdot \mathbf{\hat{n}} = 0\} \, \subset X_{eq} .
\end{equation}$$

For each set of currents we may adopt an orthogonal decomposition, which
we illustrate again for the case of the equivalent shell:
$$\begin{equation}
    X_{eq} = X_{eq,0} \oplus X_{eq,D}
\end{equation}$$ It is possible to show that $$\begin{equation}
    X_{eq,0} = \ker \mathbf{{F}}_{eq}^T, \quad X_{eq,D} = \mathrm{im} \mathbf{{F}}_{eq}  \, .
\end{equation}$$ We collect as columns of the matrix $\mathbf{{K_{eq}}}$
basis vectors for $X_{eq,0}$. In particular $$\begin{equation}
    \dim X_{eq,0} = N_{eq,0} = N_{eq} - (N_{e}-1), \quad \dim X_D = N_e -1
\end{equation}$$ where $N_{eq}$ is the overall number of degrees of
freedom for the current in the equivalent shell and $N_e$ is the number
of electrodes between equivalent shell and conducting domain. Since the
overall current (*i.e.* in the whole euclidean space) has to be
solenoidal, we may specify arbtrarily the current fluxes only at $N_e-1$
electrodes. A generic equivalent current can be written explicitly in
terms of a solenoidal component and a component which determines current
fluxes through electrodes: $$\begin{equation}
    \label{eq:decomposition}
    \mathbf{{i}}_{eq} = \mathbf{{K}}_{eq} \, \mathbf{{i}}_{eq,i} + \mathbf{{F}}_{eq} \, \mathbf{{i}}_{eq,e} \, ,
\end{equation}$$ where $\mathbf{{i}}_{eq,i} \in \mathbb{R}^{N_{eq,0}}$
represents the currents purely internal to the equivalent shell and
$\mathbf{{i}}_{eq,e} \in \mathbb{R}^{N_{e}-1}$ represents current
distributions determining electrode current fluxes.

As a side note we highlight that the overall current has to be
solenoidal. In case we can clearly identify perfectly facing electrode
couples between shell and passive structures it is easy to set up the
Kirchoff Current Laws at these corresponding electrodes. If $N_e$ is the
overall number of electrodes on the boundary of each domain, we can
write a set of $N_e-1$ independent Kirchhoff's current laws, imposing
the continuity of current densities from the equivalent shell to the
actual conducting domain: $$\begin{equation}
    \label{eq:LKC}
    \mathbf{{F}}_w^T \mathbf{{i}}_w = \mathbf{{F}}_{eq}^T \mathbf{{i}}_{eq}
\end{equation}$$

# The virtual experiment

The first step of the virtual casing principle experiment consists in
identifying a certain domain $V_{in}$ which contains a current source
(IN) which we want to represent via a set of equivalent currents at the
boundary of the domain itself $\partial V_{in}$. We imagine a
**superconducting shell** (SS) on the domain boundary. Turning on the
current source within the control domain $V_{in}$ a surface current
immediately appears on the superconducting sheet which immediately
cancels any induced **electric** field in the outer domain.

In the standard case, where the current in the internal domain is
therein solenoidal, the electric field is associated solely with
variations of magnetic fluxes, *i.e.* the surface current on the sheet
cancels perfectly any magnetic field or flux variation in the outer
domain.

In the more general case where the current source is actually injecting
and draining current from the surrounding superconducting sheet the
situation is slightly more complicated. Indeed in this case the electric
field in the outer conductors is associated also to the contact
potentials due to the current injection. Hence the condition of zero
electric field on the superconducting shell (SS) shall be considered in
its integrity:

$$\begin{align}
    \mathbf{E} = - %
\frac{\partial^{}\mathbf{A}_{SS}}%
{\partial t^{}} - %
\frac{\partial^{}\mathbf{A}_{IN}}%
{\partial t^{}} - \nabla \phi_{SS} = 0 \label{eq:integral_equation_01}\\
    \nabla_S \cdot \mathbf{k}_{SS} = J_{n,IN} \label{eq:integral_equation_02}\\
    \mathbf{A}_{SS} = \frac{\mu_0}{4 \pi} \int_{\partial V_{in}}{ G(\mathbf{x},\mathbf{x'}) \mathbf{k}_{SS}(\mathbf{x'}) }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{S} \label{eq:integral_equation_03}
\end{align}$$

In particular we may integrate in time equation
([\[eq:integral_equation_01\]](#eq:integral_equation_01){reference-type="ref"
reference="eq:integral_equation_01"}) and plug
([\[eq:integral_equation_03\]](#eq:integral_equation_03){reference-type="ref"
reference="eq:integral_equation_03"}) therein: $$\begin{align}
    \frac{\mu_0}{4 \pi} \int_{\partial V_{in}}{ G(\mathbf{x},\mathbf{x'}) \mathbf{k}_{SS}(\mathbf{x'}) }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{S} + \nabla \int_{0}^{t} {\phi_{SS}(\mathbf{x,\tau})} \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{\tau} = - \Delta \mathbf{A}_{IN} \quad \mathrm{on} \, \partial V_{IN} \label{eq:integral_equation_final_01}\\
    \nabla_S \cdot \mathbf{k}_{SS} = J_{n,IN} \quad \mathrm{in} \, \partial V_{IN}\label{eq:integral_equation_final_02}
\end{align}$$

Here the unkown terms are the sheet surface current density
$\mathbf{k}_{SS}$ and the potentials $\phi_{SS}$, whereas the assigned
quantities are the magnetic vector potential due to internal sources
$\mathbf{A}_{IN}$ and the injected currents in the superconducting shell
$\mathbf{k}_{SS}$. In particular the current $\mathbf{k}_{SS}$ shields
perfectly the magnetic field due to internal currents in the outer
domain. The equivalent current on $\partial V_{IN}$ which reproduces the
internal currents is clearly $\mathbf{k}_{eq}=-\mathbf{k}_{SS}$.

# The equivalent current to the plasma

The weak form for problem
([\[eq:integral_equation_final_01\]](#eq:integral_equation_final_01){reference-type="ref"
reference="eq:integral_equation_final_01"})-([\[eq:integral_equation_final_02\]](#eq:integral_equation_final_02){reference-type="ref"
reference="eq:integral_equation_final_02"}), projecting onto the space
of CARIDDI (or STARWALL) basis functions reads as follows:
$$\begin{align}
    \mathbf{{L}}_{eq} \, \mathbf{{i}}_{eq} + \mathbf{{F}}_{eq} \, \mathbf{{\chi}}_{eq} = \mathbf{{a}}_{pl} \label{eq:equivalent_current_problem_01}\\
    \mathbf{{F}}_{eq}^T \, \mathbf{{i}}_{eq} = \mathbf{{i}}_{el,eq} \label{eq:equivalent_current_problem_02}
\end{align}$$

Here we indicated by $\mathbf{{i}}_{eq}$ the unknown degrees of freedom
associated with the plasma equivalent currents and by
$\mathbf{{\chi}}_{eq}$ the corresponding time integrated potentials on
the superconducting shell. The forcing term are represented by
$\mathbf{{a}}_{pl}$, the magnetic vector potential due to plasma
sources, and the currents injected into superconducting shell electrodes
$\mathbf{i}_{el,eq}$. In our case the plasma contribution to the
magnetic vector potential may be evaluated in terms of the overall
vector potential and external currents as follows:

$$\begin{equation}
    \mathbf{{a}}_{pl} = \mathbf{{H}} \, \mathbf{{a}}_{tot} - \mathbf{{M}}_{eq,w} \, \mathbf{{i}}_w \, .
\end{equation}$$

We stress that the inductance matrix here, as well as the projection
matrix, include all the equivalent shell degrees of freedom, both
representing solenoidal currents in the shell and currents crossing the
electrodes.

## Block inversion

A first possibility to express the equivalent currents
$\mathbf{{i}}_{eq}$ in terms of the forcing terms is given by direct
inversion of the problem
([\[eq:equivalent_current_problem_01\]](#eq:equivalent_current_problem_01){reference-type="ref"
reference="eq:equivalent_current_problem_01"})-([\[eq:equivalent_current_problem_02\]](#eq:equivalent_current_problem_02){reference-type="ref"
reference="eq:equivalent_current_problem_02"}), which in matrix form
looks like

$$\begin{equation}
    \label{eq:equivalent_current_step1}
    \left[ 
    \begin{array}{cc}
        \mathbf{{L}}_{eq} & \mathbf{{F}}_{eq} \\
        \mathbf{{F}}_{eq}^T & \mathbf{{0}}
    \end{array}
    \right] \, \left[ 
    \begin{array}{c}
        \mathbf{{i}}_{eq}  \\
        \mathbf{{\chi}}_{eq}
    \end{array} \right] = 
    \left[
    \begin{array}{c}
        \mathbf{{H}}   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{a}}_{tot} - 
    \left[
    \begin{array}{c}
        \mathbf{{M}}_{eq,w}   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{i}}_{w}   +
    \left[
    \begin{array}{c}
        \mathbf{{0}}   \\
        \mathbf{{i}}_{el,eq}
    \end{array} \right]
\end{equation}$$

We can exploit the block structure of the matrix to invert:

$$\begin{equation}
     \label{eq:block_inversion_main}
    \left[ 
    \begin{array}{cc}
        \mathbf{{L}}_{eq} & \mathbf{{F}}_{eq} \\
        \mathbf{{F}}_{eq}^T & \mathbf{{0}}
    \end{array}
    \right] ^{-1} = \left[ 
    \begin{array}{cc}
        \mathbf{{B}}_{11} & \mathbf{{B}}_{12} \\
        \mathbf{{B}}_{21} & \mathbf{{B}}_{22}
    \end{array}    
    \right]
\end{equation}$$ where: $$\begin{align}
 \mathbf{{B}}_{11} & = \mathbf{{L}}_{eq}^{-1} + \mathbf{{L}}_{eq}^{-1} \, \mathbf{{F}}_{eq} \, \mathbf{{B}}_{22} \, \mathbf{{F}}_{eq}^T \, \mathbf{{L}}_{eq}^{-1}    \\
 \mathbf{{B}}_{12} & = - \mathbf{{L}}_{eq}^{-1} \, \mathbf{{F}}_{eq} \, \mathbf{{B}}_{22}    \\
 \mathbf{{B}}_{21} & = - \mathbf{{B}}_{22} \, \mathbf{{F}}_{eq}^T \, \mathbf{{L}}_{eq}^{-1} \\
 \mathbf{{B}}_{22} & = - \left( \mathbf{{F}}_{eq}^T \, \mathbf{{L}}_{eq}^{-1} \, \mathbf{{F}}_{eq} \right)^{-1}
\end{align}$$ According to above definitions the equivalent currents to
the plasma can be expressed in terms of the *forcing terms*
$\mathbf{{a}}_{tot}$ and $\mathbf{{i}}_{el,eq}$ and of the wall currents
$\mathbf{{i}}_w$ as $$\begin{equation}
\label{eq:equivalent_current_block_inversion}
    \mathbf{{i}}_{eq} = \left( \mathbf{{B}}_{11} \, \mathbf{{H}} \right) \mathbf{{a}}_{tot} - \left( \mathbf{{B}}_{11} \, \mathbf{{M}}_{eq,w} \right) \mathbf{{i}}_{w} + \mathbf{{B}}_{12} \mathbf{{i}}_{el,eq}
\end{equation}$$

## Null scheme

Based on the decomposition into solenoidal currents and
electrodes-crossing currents introduced in equation
([\[eq:decomposition\]](#eq:decomposition){reference-type="ref"
reference="eq:decomposition"}) we may proceed at the inversion problem
per each subspace. The idea is to express the overall current as sum of
the solenoidal and cross-electrode contribution, and to project equation
([\[eq:equivalent_current_problem_01\]](#eq:equivalent_current_problem_01){reference-type="ref"
reference="eq:equivalent_current_problem_01"}) on $X_{eq,0}$, which
corresponds to left multiplication by $\mathbf{{K}}_{eq}^T$:
$$\begin{equation}
    \label{eq:equivalent_current_with_null_01}
    \left[ 
    \begin{array}{cc}
        \mathbf{{L}}_{eq,ii} & \mathbf{{L}}_{eq,ie} \\
        \mathbf{{0}} & \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq}
    \end{array}
    \right] \, \left[ 
    \begin{array}{c}
        \mathbf{{i}}_{eq,i}  \\
        \mathbf{{i}}_{eq,e}
    \end{array} \right] = 
    \left[
    \begin{array}{c}
        \mathbf{{H}}_i   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{a}}_{tot} - 
    \left[
    \begin{array}{c}
        \mathbf{{M}}_{eqi,w}   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{i}}_{w} + \left[
    \begin{array}{c}
        \mathbf{{0}}  \\
        \mathbf{{1}} 
    \end{array}
    \right] \, \mathbf{{i}}_{el,eq}
\end{equation}$$ The intuitive definition of the matrices appearing in
equation
([\[eq:equivalent_current_step1\]](#eq:equivalent_current_step1){reference-type="ref"
reference="eq:equivalent_current_step1"}) is reported in the Appendix,
see equations
([\[eq:definition_1a\]](#eq:definition_1a){reference-type="ref"
reference="eq:definition_1a"})-([\[eq:definition_1h\]](#eq:definition_1h){reference-type="ref"
reference="eq:definition_1h"}). In particular we can solve first for the
electric current density which actually cross electrodes:
$$\begin{equation}
\label{eq:LKC_explicit}
    \left( \mathbf{{F_{eq}}}^T \, \mathbf{{F_{eq}}} \right ) \, \mathbf{{i}}_{eq,e} = \mathbf{{i}}_{el,eq} \, .
\end{equation}$$

Given equation
([\[eq:LKC_explicit\]](#eq:LKC_explicit){reference-type="ref"
reference="eq:LKC_explicit"}), you can proceed at finding the solenoidal
part of the equivalent current density distribution: $$\begin{equation}
    \mathbf{{i}}_{eq,i} = \left( \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_i \right) \mathbf{{a}}_{tot} - \left( \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \right) \mathbf{{i}}_w - \left[ \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eq,ie} \, \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1} \right] \mathbf{{i}}_{el,eq}
\end{equation}$$ Finally combining the solenoidal and cross-electrodes
current via equation
([\[eq:decomposition\]](#eq:decomposition){reference-type="ref"
reference="eq:decomposition"}): $$\begin{align}
\label{eq:equivalent_current_with_null_02}
    \mathbf{{i}}_{eq} = & \left( \mathbf{{K}}_{eq} \,\mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_i \right) \mathbf{{a}}_{tot} - \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \right) \mathbf{{i}}_w \nonumber \\
    & - \left[ \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eq,ie} - \mathbf{{F}}_{eq} \right)\, \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1} \right] \mathbf{{i}}_{el,eq}
\end{align}$$ Clearly the matrices defined in equation
([\[eq:equivalent_current_block_inversion\]](#eq:equivalent_current_block_inversion){reference-type="ref"
reference="eq:equivalent_current_block_inversion"}) and equation
([\[eq:equivalent_current_with_null_02\]](#eq:equivalent_current_with_null_02){reference-type="ref"
reference="eq:equivalent_current_with_null_02"}) shall be identical.

As a side note, we mention that in principle we may project our
equivalent current equation
([\[eq:equivalent_current_problem_01\]](#eq:equivalent_current_problem_01){reference-type="ref"
reference="eq:equivalent_current_problem_01"}) as well on the space of
cross-electrodes currents in order to extrapolate the time integrated
potential $\mathbf{{\chi}}_{eq}$.

$$\begin{equation}
\label{eq:equivalent_current_per_components}
    \left[ \begin{array}{ccc}
        \mathbf{{L}}_{eq,ii} & \mathbf{{L}}_{eq,ie} & 0  \\
        \mathbf{{L}}_{eq,ei} & \mathbf{{L}}_{eq,ee} & \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq}  \\
        0 &   \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & 0 \\
    \end{array} \right] \, \left[ 
    \begin{array}{c}
         \mathbf{{i}}_{eq,i}  \\
         \mathbf{{i}}_{eq,e}  \\
         \mathbf{{\chi}}_{eq} \\
    \end{array}
    \right]  = 
    \left[
    \begin{array}{c}
        \mathbf{{H_i}}   \\
        \mathbf{{H_e}}   \\
        \mathbf{{0}}     \\
    \end{array}
    \right] \, \mathbf{{a}}_{tot} -
    \left[
    \begin{array}{c}
        \mathbf{{M}}_{eqi,w}   \\
        \mathbf{{M}}_{eqe,w}   \\
        \mathbf{{0}}            \\
    \end{array}
    \right] \, \mathbf{{i}}_{w}   +
    \left[ \begin{array}{c}
        \mathbf{{0}}  \\
        \mathbf{{0}}  \\
        \mathbf{{1}}
    \end{array} 
    \right] \, \mathbf{{i}}_{el,eq}
\end{equation}$$ The matrices here are defined again in Appendix, see
equations
([\[eq:definition_2a\]](#eq:definition_2a){reference-type="ref"
reference="eq:definition_2a"})-([\[eq:definition_2c\]](#eq:definition_2c){reference-type="ref"
reference="eq:definition_2c"}). The possibility to eliminate some
constraints above is strictly related to the possibility of determining
some unknown variable in advance. The imposition of electrode currents
via the last set of equations could be in principle avoided if the time
integrated potentials were assigned.

It is worth to notice that a formulation is possible where only the
electrode currents in the shell or the electrode currents in the wall
are specified, rather than both. In this case the wall-electrode to
shell-electrode Kirchoff Current Laws shall be specified. In the simple
case mentioned in equation ([\[eq:LKC\]](#eq:LKC){reference-type="ref"
reference="eq:LKC"}), where there is a one-to-one correspondence between
shell and wall electrodes, this would reduce to $$\begin{equation}
    \label{eq:LKC_equivalent_wall}
    (\mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq}) \mathbf{{i}}_{eq,e} =     (\mathbf{{F}}_{w}^T \, \mathbf{{F}}_{w}) \mathbf{{i}}_{w,e} \, .
\end{equation}$$

We will specify the current injected into shell and wall electrodes
separately, hence we will need to take particular care in guaranteeing
the solenoidality of the overall current density distribution.

# External currents evolution

The general evolution equations for the currents in the external domain
are given by $$\begin{align}
\label{eq:wall currents equation}\mathbf{{L}}_w \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_w + \mathbf{{R}}_w \, \mathbf{{i}}_w + \mathbf{{M}}_{w,eq} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{eq} + \mathbf{{F}}_w \mathbf{{\phi}}_w = \mathbf{{v}} \\
\label{eq:wall source-sink relation}\mathbf{{F}}_w^T \, \mathbf{{i}}_w = \mathbf{{i}}_{el,w}
\end{align}$$ Here we have introduced explicitly the currents injected
into the wall electrodes, assumed as forcing term from JOREK, *i.e.*
$\mathbf{{i}}_{el,w}$. In previous Section we have illustrated how to
relate the equivalent currents $\mathbf{{i}}_{eq}$ directly to primary
variables. Considering equations
([\[eq:equivalent_current_block_inversion\]](#eq:equivalent_current_block_inversion){reference-type="ref"
reference="eq:equivalent_current_block_inversion"}) and
([\[eq:equivalent_current_with_null_02\]](#eq:equivalent_current_with_null_02){reference-type="ref"
reference="eq:equivalent_current_with_null_02"}), We may define in
particular the matrices: $$\begin{align}
    \mathbf{{B}} & = - \mathbf{{B}}_{11} \, \mathbf{{M}}_{eq,w} = - \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \\
    \mathbf{{C}} & = \mathbf{{B}}_{11} \mathbf{{H}} = \mathbf{{K}}_{eq} \,\mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_i \\
    \mathbf{{D}} & = -\mathbf{{B}}_{12} = -\left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eq,ie} - \mathbf{{F}}_{eq} \right)\, \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1}
\end{align}$$ we have synthetically $$\begin{equation}
    \mathbf{{i}}_{eq} = \mathbf{{B}} \, \mathbf{{i}}_{w} + \mathbf{{C}} \, \mathbf{{a}}_{tot} + \mathbf{{D}} \, \mathbf{{i}}_{el,eq} \, .
\end{equation}$$ If we define $$\begin{align}
\mathbf{{B}}^* = -\mathbf{{M}}_{w,eq} \, \mathbf{{B}}, \quad \mathbf{{C}}^* = -\mathbf{{M}}_{w,eq} \, \mathbf{{C}}, \quad \mathbf{{D}}^* = -\mathbf{{M}}_{w,eq} \, \mathbf{{D}}, \\
\mathbf{{L}}_w^* = \mathbf{{L}}_w - \mathbf{{B}}^*,
\end{align}$$ the evolution equation for wall currents in matrix form
can be written: $$\begin{equation}
    \left[ 
    \begin{array}{cc}
        \mathbf{{L}}_w^* & \mathbf{{0}} \\
        \mathbf{{0}} & \mathbf{{0}}
    \end{array}
    \right] \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[
    \begin{array}{c}
         \mathbf{{i}}_w  \\
         \mathbf{{\phi}}_w 
    \end{array}
    \right] + 
    \left[ 
    \begin{array}{cc}
        \mathbf{{R}}_w & \mathbf{{F}}_w \\
        \mathbf{{F}}_w^T & \mathbf{{0}}
    \end{array}
    \right] \, \left[
    \begin{array}{c}
         \mathbf{{i}}_w  \\
         \mathbf{{\phi}}_w 
    \end{array}
    \right] = 
    \left[ 
    \begin{array}{cc}
        \mathbf{{C}}^{*} & \mathbf{{D}}^{*} \\
        \mathbf{{0}}  & \mathbf{{0}}
    \end{array}
    \right] \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[
    \begin{array}{c}
         \mathbf{{a}}_{tot}  \\
         \mathbf{{i}}_{el,eq} 
    \end{array}
    \right] + 
        \left[ 
    \begin{array}{c}
        \mathbf{{0}}  \\
        \mathbf{{1}} 
    \end{array}
    \right] \, 
    \mathbf{{i}}_{el,w}
\end{equation}$$ where $\mathbf{{\phi}}$ is the vector of applied
potentials to the electrodes of interface between equivalent shell and
conducting domain. Projection onto $X_{w,0}$ results clearly in
$$\begin{align}
\mathbf{{L}}^{*}_{w,ii} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = - \mathbf{{R}}_{w,ie} \, \mathbf{{i}}_{w,e} - \mathbf{{L}}^{*}_{w,ie} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,e} + \mathbf{{C}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot}   + \mathbf{{D}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,eq} + \mathbf{{K}}_{w} \, \mathbf{{v}} \label{eq:evolution_current_null_01} \\
\mathbf{{F}}_w^T  \mathbf{{F}}_w \, \mathbf{{i}}_{w,e} = \mathbf{{i}}_{el,w} \label{eq:evolution_current_null_02}
\end{align}$$

Inversion of equation
([\[eq:evolution_current_null_02\]](#eq:evolution_current_null_02){reference-type="ref"
reference="eq:evolution_current_null_02"}) allows to determinate the
last forcing term at the r.h.s. of equation
([\[eq:evolution_current_null_01\]](#eq:evolution_current_null_01){reference-type="ref"
reference="eq:evolution_current_null_01"}): $$\begin{equation}
    \mathbf{{i}}_{w,e} = \left( \mathbf{{F}}_w^T  \mathbf{{F}}_w \right)^{-1} \, \mathbf{{i}}_{el,w}
\end{equation}$$ Defining finally $$\begin{align}
    \mathbf{{E}}_i^* = - \mathbf{{L}}_{w,to tie}^{*} \, \left( \mathbf{{F}}_w^T  \mathbf{{F}}_w \right)^{-1} \\
    \mathbf{{G}}_i^* = - \mathbf{{R}}_{w,ie} \, \left( \mathbf{{F}}_w^T  \mathbf{{F}}_w \right)^{-1}
\end{align}$$ we get the evolution equation $$\begin{equation}
    \mathbf{{L}}^{*}_{w,ii} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = \mathbf{{C}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot} + \mathbf{{D}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,eq} + \mathbf{{E}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,w} + \mathbf{{G}}_i^* \, \mathbf{{i}}_{el,w} + \mathbf{{K}}_{w} \, \mathbf{{v}} \label{eq:evolution_eq_final}
\end{equation}$$ Modal decomposition of the dynamical problem formulated
as in
([\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="ref"
reference="eq:evolution_eq_final"}) is finally possible, via the
classical algorithms for the evaluation of generalized eigenvalues and
eigenvectors.

For the sake of clarity we resume the definitions of matrices appearing
in
([\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="ref"
reference="eq:evolution_eq_final"}) in terms of standard CARIDDI
computable quantities in the following Table
[1](#tab:resume_final_equation){reference-type="ref"
reference="tab:resume_final_equation"}.

::: {#tab:resume_final_equation}
           Matrix                                                                                                           Definition                                                                                                         Dimensions
  ------------------------- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------------------
   $\mathbf{{L}}_{w,ii}^*$                                   $\mathbf{{K_w}}^T \left[ \mathbf{{L}}_{w} - \mathbf{{M}}_{w,eqi} \mathbf{{L}}_{eq,ii}^{-1} \mathbf{{M}}_{eqi,w}  \right] \mathbf{{K}}_w$                                    $N_{w,i}$ x $N_{w,i}$
    $\mathbf{{R}}_{w,ii}$                                                                       $\mathbf{{K_w}}^T \left[ \mathbf{{R}}_{w}  \right] \mathbf{{K}}_w$                                                                       $N_{w,i}$ x $N_{w,i}$
    $\mathbf{{C}}_{i}^*$                                                    $- \mathbf{{K}}_w^T \left[ \mathbf{{M}}_{w,eqi} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_{i}  \right]$                                                    $N_{w,i}$ x $N_{J}$
    $\mathbf{{D}}_{i}^*$     $\mathbf{{K}}_w^T \left[  \left(  \mathbf{{M}}_{w,eqi} \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eqi,e} - \mathbf{{M}}_{w,eqe} \right)\right] \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1}$   $N_{w,i}$ x $N_{el,J}$
    $\mathbf{{E}}_{i}^*$       $- \mathbf{{K_w}}^T \left[ \mathbf{{L}}_{w} - \mathbf{{M}}_{w,eqi} \, \mathbf{{L}}_{eq,ii}^{-1} \mathbf{{M}}_{eqi,w}  \right] \mathbf{{F}}_w \left( \mathbf{{F}}_{w}^T \mathbf{{F}}_{w} \right)^{-1}$     $N_{w,i}$ x $N_{el,w}$
    $\mathbf{{G}}_{i}^*$                                           $- \mathbf{{K_w}}^T \left[ \mathbf{{R}}_{w}  \right] \mathbf{{F}}_w \left( \mathbf{{F}}_{w}^T \mathbf{{F}}_{w} \right)^{-1}$                                          $N_{w,i}$ x $N_{el,w}$

  : Resume of matrices appearing in the evolution equation
  ([\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="ref"
  reference="eq:evolution_eq_final"}). $N_{w,i} = N_{w} - N_{el,w}$ is
  the number of degrees of freedom for solenoidal currents in the
  conducting domain, where $N_{el,w}$ is the number of *independent*
  wall electrodes (*e.g.* for a continuous electrical conductor the
  overall number of electrodes except one). $N_J$ is the number of
  degrees of freedom for the magnetic vector potential in JOREK.
:::

Note that, in the static limit for the wall, the above equation
[\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="eqref"
reference="eq:evolution_eq_final"} is reduced to $$\begin{equation}
    \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = \mathbf{{G}}^*_i \, \mathbf{{i}}_{el,w}
\end{equation}$$

# Eigenvalue problem

The matrices $\mathbf{{L}}_{w,ii}^*$ and $\mathbf{{R}}_{w,ii}$ are used
to solve the reduced-space eigenvalue problem:

$$\begin{equation}
 \label{eq:reduced_eigs}
    \mathbf{{L}}_{w,ii}^* \, \mathbf{{S}}_{\lambda} = \lambda \, \mathbf{{R}}_{w,ii} \, \mathbf{{S}}_{\lambda}
\end{equation}$$

But what if we solve the eigenvalue problem in the full space? Meaning
$$\begin{equation}
\label{eq:full_eigs}
    \mathbf{{L}}_{w}^* \, \mathbf{{S}}_{\lambda} = \lambda \, \mathbf{{R}}_{w} \, \mathbf{{S}}_{\lambda}
\end{equation}$$ In equation
[\[eq:reduced_eigs\]](#eq:reduced_eigs){reference-type="eqref"
reference="eq:reduced_eigs"} we know that we are solving the eigenvalues
of just the circuit in the \"internal\" space, and we assign the rest of
the eigenvalues accounting for the fact that the L-R circuit has a
current impressed by $\mathbf{{i}}_{el,w}$ in the nodes. In case
[\[eq:full_eigs\]](#eq:full_eigs){reference-type="eqref"
reference="eq:full_eigs"} we instead consider the circuit as
voltage-impressed by the unknows $\Phi_w$.\
**ToDo**: *Check if the eigenvalues obtained via reduced system are a
subset of the full ones, or if they are completely different due to the
fact that one system is current-impressed and the other is
voltage-impressed. The latter reason should indeed change the state
variable space of our R-L circuit.*

## Alternative solution method via augmented system {#alternative-solution-method-via-augmented-system .unnumbered}

An \"augmented\" system can be constructed for the solution of the same
problem as
[\[eq:equivalent_current_with_null_02\]](#eq:equivalent_current_with_null_02){reference-type="eqref"
reference="eq:equivalent_current_with_null_02"},
[\[eq:wall currents equation\]](#eq:wall currents equation){reference-type="eqref"
reference="eq:wall currents equation"},
[\[eq:wall source-sink relation\]](#eq:wall source-sink relation){reference-type="eqref"
reference="eq:wall source-sink relation"}.\
This choice mainly hinges on the line of reasoning that the inversion of
$\mathbf{{L}}_{eq,ii}$, which has a $O(n^3)$ complexity and it is
possibly suffering from ill-conditioning, could be avoided. Since the
latter comes from the need to eliminate $\mathbf{{i}}_{eq}$ as a
variable to obtain Equation
[\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="eqref"
reference="eq:evolution_eq_final"}, the augmented method consists in
circumventing this problem by including $\mathbf{{i}}_{eq}$ to the
vector of unknowns, at the obvious cost of dealing with the inversion a
larger system.\
Let construct the augmented vector: $$\begin{equation}
    \mathbf{{i}}_{a} = \left[ 
    \begin{array}{cc}  
\mathbf{{i}}_w \\
\mathbf{{i}}_{eq}
\end{array}
    \right]
\end{equation}$$ then the time evolution equation for this augmented
system becomes $$\begin{align}
\label{eq:augmented time evol 1}\mathbf{{L}}_a %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_a + \mathbf{{R}}_a \mathbf{{i}}_a + \mathbf{{F}}_a \mathbf{{\Phi}}_a = \mathbf{{v}}_a \\
\label{eq:augmented source-sink 1}\mathbf{{F}}^T_a \mathbf{{i}}_a = \mathbf{{i}}_{el,a}
\end{align}$$ The matrix form is: $$\begin{align}
\left[ \begin{array}{cc}  
\mathbf{{L}}_w & \mathbf{{M}}_{w,eq} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{R}}_w & \mathbf{{0}} \\
\mathbf{{M}}_{eq,w}   & \mathbf{{L}}_{eq} 
\end{array}
\right] \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{F}}_w & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{F}}_{eq} 
\end{array} \right] \left[ \begin{array}{cc} \mathbf{{\Phi}}_w \\
\mathbf{{\chi}}_{eq} \end{array} \right] = \left[ \begin{array}{cc} \mathbf{{v}} \\
\mathbf{{H}} \, \mathbf{{a}}_{tot} \end{array}  
\right] \\
\left[ \begin{array}{cc} \mathbf{{F}}^T_w & \mathbf{{0}} \\
\mathbf{{0}} & \mathbf{{F}}^T_{eq} \end{array} \right] \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] =  \left[ \begin{array}{cc} \mathbf{{i}}_{el,w} \\
\mathbf{{i}}_{el,eq} \end{array} 
\right]
\end{align}$$ We can hence construct the augmented matrix whose columns
span the generalized null space $X_{a,0}$: $$\begin{equation}
    \mathbf{{K}}_a = \left[ \begin{array}{cc} \mathbf{{K}}_w & \mathbf{{0}} \\
\mathbf{{0}} & \mathbf{{K}}_{eq} \end{array} 
\right]
\end{equation}$$ and perform left multiplication by $\mathbf{{K}}^T_a$.
The projection to $X_{a,0}$ yields: $$\begin{equation}
\begin{split}
    \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ii} & \mathbf{{M}}_{wi,eqi} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_{w,i} \\
\mathbf{{i}}_{eq,i} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ie} & \mathbf{{M}}_{wi,eqe} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] + \\
+ \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ii} & \mathbf{{0}} \\
\mathbf{{M}}_{eqi,wi}   & \mathbf{{L}}_{eq,ii} 
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{i}}_{w,i} \\
\mathbf{{i}}_{eq,i} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ie} & \mathbf{{0}} \\
\mathbf{{M}}_{eqi,we}   & \mathbf{{L}}_{eq,ie} 
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] = \left[ \begin{array}{cc}  
\mathbf{{K}}^T_{w} & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{K}}^T_{eq} 
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{v}} \\
\mathbf{{H}} \, \mathbf{{a}}_{tot} \end{array} 
\right]
\end{split}
\end{equation}$$ $$\begin{equation}
    \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] = \left[ \begin{array}{cc}  
\left( \mathbf{{F}}^T_{w} \mathbf{{F}}_{w} \right)^{-1} & \mathbf{{0}} \\
\mathbf{{0}}   & \left( \mathbf{{F}}^T_{eq} \mathbf{{F}}_{eq} \right)^{-1}
\end{array} \right] \left[ \begin{array}{cc} \mathbf{{i}}_{el,w} \\
\mathbf{{i}}_{el,eq} \end{array} 
\right]
\end{equation}$$

With a compressed notation, one can write: $$\begin{equation}
    \mathbf{{L}}_{a,ii} %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}\mathbf{{i}}_{a,i}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} +  \mathbf{{L}}_{a,ie} \,\mathbf{{F}}^{-1}_2 %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}\mathbf{{i}}_{el,a}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} + \mathbf{{R}}_{a,ii} \mathbf{{i}}_{a,i} + \mathbf{{R}}_{a,ie} \, \mathbf{{F}}^{-1}_2 \, \mathbf{{i}}_{el,a} = \mathbf{{v}}_a
\end{equation}$$

**[Augmented matrices:]{.underline}**

$$\begin{equation}
     \mathbf{{L}}_{a,ii} = \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ii} & \mathbf{{M}}_{wi,eqi} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]
\end{equation}$$

$$\begin{equation}
     \mathbf{{L}}_{a,ie} = \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ie} & \mathbf{{M}}_{wi,eqe} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]
\end{equation}$$

$$\begin{equation}
     \mathbf{{F}}^{-1}_{2} = \left[ \begin{array}{cc}  
\left( \mathbf{{F}}^T_{w} \mathbf{{F}}_{w} \right)^{-1} & \mathbf{{0}} \\
\mathbf{{0}}   & \left( \mathbf{{F}}^T_{eq} \mathbf{{F}}_{eq} \right)^{-1}
\end{array} \right]
\end{equation}$$

$$\begin{equation}
    \mathbf{{R}}_{a,ii} = \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ii} & \mathbf{{0}} \\
\mathbf{{M}}_{eqi,wi}   & \mathbf{{L}}_{eq,ii} 
\end{array}
\right]
\end{equation}$$ $$\begin{equation}
    \mathbf{{R}}_{a,ie} = \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ie} & \mathbf{{0}} \\
\mathbf{{M}}_{eqi,we}   & \mathbf{{L}}_{eq,ie} 
\end{array}
\right]
\end{equation}$$

$$\begin{equation}
    \mathbf{{v}}_{a} =  \left[ \begin{array}{cc} \mathbf{{K}}^T_w \, \mathbf{{v}} \\
\mathbf{{K}}^T_{eq} \,\mathbf{{H}} \, \mathbf{{a}}_{tot} \end{array} 
\right]
\end{equation}$$

**Alternative form of the augmented system**

$$\begin{align}
\left[ \begin{array}{cc}  
\mathbf{{L}}_w & \mathbf{{M}}_{w,eq} \\
\mathbf{{M}}_{eq,w}   & \mathbf{{L}}_{eq}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{F}}_w & \mathbf{{0}} \\
\mathbf{{0}}_{}   & \mathbf{{F}}_{eq}
\end{array}
\right] \left[ \begin{array}{cc} \mathbf{{\Phi}}_w \\
\mathbf{{\Phi}}_{eq} \end{array} \right] + \left[ \begin{array}{cc}  
\mathbf{{R}}_w & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{0}} 
\end{array}
\right] \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] = \left[ \begin{array}{cc} \mathbf{{v}} \\
\mathbf{{H}} \, \mathbf{{\dot{a}}}_{tot} \end{array}  
\right] \\
\left[ \begin{array}{cc} \mathbf{{F}}^T_w & \mathbf{{0}} \\
\mathbf{{0}} & \mathbf{{F}}^T_{eq} \end{array} \right] \left[ \begin{array}{cc} \mathbf{{i}}_w \\
\mathbf{{i}}_{eq} \end{array} 
\right] =  \left[ \begin{array}{cc} \mathbf{{i}}_{el,w} \\
\mathbf{{i}}_{el,eq} \end{array} 
\right]
\end{align}$$

**[Augmented matrices of alternative form]{.underline}**

$$\begin{equation}
\begin{split}
    \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ii} & \mathbf{{M}}_{wi,eqi} \\
\mathbf{{M}}_{eqi,wi}   & \mathbf{{L}}_{eq,ii}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_{w,i} \\
\mathbf{{i}}_{eq,i} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{L}}_{w,ie} & \mathbf{{M}}_{wi,eqe} \\
\mathbf{{M}}_{eqi,we}   & \mathbf{{L}}_{eq,ie}
\end{array}
\right]  %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] + \\
+ \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ii} & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{i}}_{w,i} \\
\mathbf{{i}}_{eq,i} \end{array} 
\right] + \left[ \begin{array}{cc}  
\mathbf{{R}}_{w,ie} & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{0}}
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] = \left[ \begin{array}{cc}  
\mathbf{{K}}^T_{w} & \mathbf{{0}} \\
\mathbf{{0}}   & \mathbf{{K}}^T_{eq} 
\end{array}
\right]  \left[ \begin{array}{cc} \mathbf{{v}} \\
\mathbf{{H}} \, \mathbf{{\dot{a}}}_{tot} \end{array} 
\right]
\end{split}
\end{equation}$$ $$\begin{equation}
    \left[ \begin{array}{cc} \mathbf{{i}}_{w,e} \\
\mathbf{{i}}_{eq,e} \end{array} 
\right] = \left[ \begin{array}{cc}  
\left( \mathbf{{F}}^T_{w} \mathbf{{F}}_{w} \right)^{-1} & \mathbf{{0}} \\
\mathbf{{0}}   & \left( \mathbf{{F}}^T_{eq} \mathbf{{F}}_{eq} \right)^{-1}
\end{array} \right] \left[ \begin{array}{cc} \mathbf{{i}}_{el,w} \\
\mathbf{{i}}_{el,eq} \end{array} 
\right]
\end{equation}$$

# Tangent magnetic field (explicit coupling)

We here present the expression of the tangent magnetic field in presence
of halo currents.\
Note that, in this more general picture, the $\mathbf{B_{tan}}$
expression is sensitive to:

- The injected currents (eq.
  [\[eq:evolution_current_null_01\]](#eq:evolution_current_null_01){reference-type="eqref"
  reference="eq:evolution_current_null_01"});

- The effect of eddy currents as before;

- The coupling terms between eddy and halo in
  [\[eq:evolution_eq_final\]](#eq:evolution_eq_final){reference-type="eqref"
  reference="eq:evolution_eq_final"}.

The magnetic field in the general picture is given by (we consider the
poloidal component as an example but a similar logic is applied to the
toroidal one): $$\begin{equation}
    \underline{\mathbf{B_{\theta}}} =  \underline{\underline{\mathbf{B^*_{\theta,w} }}} \, \underline{\mathbf{I_w}} + \underline{\underline{\mathbf{C^*_{\theta,\phi}}}} \, \underline{\mathbf{\psi}} + \underline{\underline{\mathbf{C^*_{\theta,\theta} }}} \, \underline{\mathbf{A_{\theta}}}
\end{equation}$$

where $\underline{\underline{\mathbf{B^*_{\theta,w} }}}$ is modified by
the interaction of shell currents with wall currents: $$\begin{equation}
 \label{eq:Bthetastar matrix eddys}
    \underline{\underline{\mathbf{B^*_{\theta,w}}}} = \underline{\underline{\mathbf{B_{\theta,w}}}} - \underline{\underline{\mathbf{B_{\theta,eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq}}}} \, \underline{\underline{\mathbf{M_{eq,w}}}}
\end{equation}$$ and this is coming from the equation for the shell
currents in the eddy-only case $$\begin{equation}
    \underline{\underline{\mathbf{L_{eq}}}} \, \underline{\mathbf{I_{eq}}} = \underline{\underline{\mathbf{H}}} \, \underline{\mathbf{A_{tot}}} - \underline{\underline{\mathbf{M_{eq,w}}}} \, \underline{\mathbf{I_w}}
\end{equation}$$

In the halo+eddy picture, the equation becomes the following:
$$\begin{equation}
    \underline{\underline{\mathbf{L_{eq}}}} \, \underline{\mathbf{I_{eq}}} + \underline{\underline{\mathbf{F_{eq}}}} \, \underline{\mathbf{\chi_{eq}}} = \underline{\underline{\mathbf{H}}} \, \underline{\mathbf{A_{tot}}} - \underline{\underline{\mathbf{M_{eq,w}}}} \, \underline{\mathbf{I_w}}
\end{equation}$$ but the $\underline{\mathbf{\chi_{eq}}}$ potential
doesn't affect the tangent magnetic field since the curl operator is
applied.\
We also know from
[\[eq:equivalent_current_with_null_02\]](#eq:equivalent_current_with_null_02){reference-type="eqref"
reference="eq:equivalent_current_with_null_02"} that $$\begin{align}
    \mathbf{{i}}_{eq} = & \left( \mathbf{{K}}_{eq} \,\mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_i \right) \mathbf{{a}}_{tot} - \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \right) \mathbf{{i}}_w \nonumber \\
    & - \left[ \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eq,ie} - \mathbf{{F}}_{eq} \right)\, \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1} \right] \mathbf{{i}}_{el,eq}
\end{align}$$ and $$\begin{align}
    \mathbf{{i}}_{eq} = & \left( \mathbf{{K}}_{eq} \,\mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{H}}_i \right) \mathbf{{a}}_{tot} - \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \right) \mathbf{{K}}_w \, \mathbf{{i}}_{w,i} - \nonumber \\
    & - \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{M}}_{eqi,w} \right) \mathbf{{F}}_w \, \mathbf{{i}}_{w,e}  - \left[ \left( \mathbf{{K}}_{eq} \, \mathbf{{L}}_{eq,ii}^{-1} \, \mathbf{{L}}_{eq,ie} - \mathbf{{F}}_{eq} \right)\, \left( \mathbf{{F}}_{eq}^T \mathbf{{F}}_{eq} \right)^{-1} \right] \mathbf{{i}}_{el,eq}
\end{align}$$ Now imagine we write: $$\begin{equation}
    \underline{\mathbf{B_{\theta}}} =  \underline{\underline{\mathbf{B_{\theta,w} }}} \, \underline{\mathbf{i_w}} + \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\mathbf{i_{eq}}}
\end{equation}$$

Hence the equation in transformed space becomes $$\begin{multline}
     \underline{\mathbf{B_{\theta}}} = \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, \underline{\underline{\mathbf{K_w}}} \right) \, \underline{\mathbf{i_{w,i}}} + \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, \underline{\underline{\mathbf{F_w}}} \right) \, \underline{\mathbf{i_{w,e}}}  +   \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}} \right) \underline{\mathbf{a_{tot}}} - \\
     \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \underline{\underline{\mathbf{K_w}}}   \right) \underline{\mathbf{i_{w,i}}} - \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \underline{\underline{\mathbf{F_w}}}   \right) \underline{\mathbf{i_{w,e}}} - \\
     \left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}}
\end{multline}$$

which is in turn $$\begin{multline}
     \underline{\mathbf{B_{\theta}}} = \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, \underline{\underline{\mathbf{K_w}}} \, \underline{\underline{\mathbf{S_{w,ii}}}} \right) \, \underline{\mathbf{Y_i}} + \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, \underline{\underline{\mathbf{F_w}}} \right) \, \underline{\mathbf{i_{w,e}}}  +   \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}} \right) \underline{\mathbf{a_{tot}}} - \\
     \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \underline{\underline{\mathbf{K_w}}}  \,\underline{\underline{\mathbf{S_{w,ii}}}} \right) \underline{\mathbf{Y_{i}}} - \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \underline{\underline{\mathbf{F_w}}}   \right) \underline{\mathbf{i_{w,e}}} - \\
     \left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}}
\end{multline}$$ and in more compact form $$\begin{multline}
     \underline{\mathbf{B_{\theta}}} = \left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \right) \underline{\underline{\mathbf{K_w}}} \, \underline{\underline{\mathbf{S_{w,ii}}}} \right] \, \underline{\mathbf{Y_i}} + \\
     \left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}   \right) \underline{\underline{\mathbf{F_w}}} \right] \, \underline{\mathbf{i_{w,e}}}  +   \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}} \right) \underline{\mathbf{a_{tot}}}  - \\
     \left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}}
\end{multline}$$

$$\begin{multline}
     \underline{\mathbf{B_{\theta}}} = \left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \right) \underline{\underline{\mathbf{K_w}}} \, \underline{\underline{\mathbf{S_{w,ii}}}} \right] \, \underline{\mathbf{Y_i}} + \\
    \color{blue} \left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}   \right) \underline{\underline{\mathbf{F_w}}} \right] \, \underline{\mathbf{i_{w,e}}}  \color{black} +   \left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}} \right) \underline{\mathbf{a_{tot}}}  - \\
    \color{teal} \left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}}
\end{multline}$$

The new field matrices are

$$\begin{equation}
     \underline{\underline{\mathbf{B^*_{\theta,wi}}}} =  \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}
\end{equation}$$

$$\begin{equation}
     \underline{\underline{\mathbf{B^*_{\theta,eqi}}}} =  \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}}
\end{equation}$$

# Electrode currents in JOREK time-stepping (implicit coupling)

Let consider the two contributions to $B_{tan}$:
$$\left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}   \right) \underline{\underline{\mathbf{F_w}}} \right] \, \underline{\mathbf{i_{w,e}}}$$\
and

$$\left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} + \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}}$$

We now have an expression for the magnetic field that is:
$$\begin{equation}
    \underline{\mathbf{B}}_{tan} = \underline{\underline{\mathbf{M}}}_{e,yi} \, \underline{\mathbf{Y}}_i + \underline{\underline{\mathbf{M}}}_{e,ye} \, \underline{\mathbf{I}}_e + \underline{\underline{\mathbf{M}}}_{e,e} \, \underline{\mathbf{a}}_{tot} - \underline{\underline{\mathbf{{M}_{e,el}}}} \, \underline{\mathbf{i}}_{el,eq}
\end{equation}$$ $$\begin{equation}
    \mathbf{{B}}_{tan} = \mathbf{{M}}_{e,wi} \, \mathbf{{i}}_{w,i} + \mathbf{{M}}_{e,we} \, \mathbf{{i}}_{w,e} + \mathbf{{M}}_{e,e} \, \mathbf{{a}}_{tot} - \mathbf{{M}}_{e,el} \, \mathbf{{i}}_{el,eq}
\end{equation}$$ To put it in JOREK's reduced MHD perspective, as in
equation 2 of [@hoelzl2012coupling]:

$$\begin{equation}
    B_{tan} = \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \Psi_j + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot Y_{k_i} + \sum_{k_e} \hat{M}^{ey_e}_{i,k_e} \cdot I_{k_e}  - \sum_{el_J} \hat{M}^{e,el}_{i,el_J} \cdot I_{el_J}\right)
\end{equation}$$ Let consider that $$\begin{align}
\underline{\mathbf{I}}_{w,e} = \left( \underline{\underline{\mathbf{F}}}^T_w \, \underline{\underline{\mathbf{F}}}_w  \right)^{-1} \underline{\mathbf{I}}_{el,w} = \underline{\underline{\mathbf{P}}}^w \, \underline{\mathbf{I}}_{el,w}
\end{align}$$ hence yielding $$\begin{equation}
    B_{tan} = \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \Psi_j + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot Y_{k_i} + \sum_{k_e} \hat{M}^{ey_e}_{i,k_e}  \cdot \sum_{i_e} P^w_{k_e,i_e} \, I_{el,w_{ie}}    - \sum_{el_{J_i}} \hat{M}^{e,el}_{i,el_{J_i}} \cdot I_{el_{J_i}} \right)
\end{equation}$$ **How are these electrode currents derived?**\
$$\begin{align}
 \label{eq:JntoIelec}
    I_{el,w} =& \underline{\underline{\mathbf{W}}}^w * J_n \Big\rvert_{bez} =\underline{\underline{\mathbf{W}}}^w \cdot \mathcal{F} \left( j \cdot \psi_s   \right) \\
    I_{el,J} =& \underline{\underline{\mathbf{W}}}^J * J_n \Big\rvert_{bez} =\underline{\underline{\mathbf{W}}}^J \cdot \mathcal{F}\left( j \cdot \psi_s   \right)
\end{align}$$

The $J_n\Big\rvert_{bez}$ loses the separate information on $j$ and
$\psi_s$!\
Let consider that the normal current $J_n$ is initially computed in a
generic point in **real space** as: $$\begin{equation}
    J_{n,\mathcal{R}} (x) = \frac{j(x) \, \psi_s (x) + R^2(x) \, p_s(x) }{\sqrt{R_s^2(x) + Z_s^2(x)} \, F_0 \, R(x) } \, \text{sgn}(x) \approx \frac{j(x) \, \psi_s (x) }{\sqrt{R_s^2(x) + Z_s^2(x)} \, F_0 \, R(x) } \, \text{sgn}(x)
\end{equation}$$ which implies that we must factor in the coefficients
that multiply $j \, \psi_s$.\
In the Bezier space, we have:

$$\begin{equation}
    \delta J_{n,\mathcal{B}} = \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}} \, \delta J_{n,\mathcal{R}}
\end{equation}$$ where
$\underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}$ is the matrix
performing the Fourier+Bezier transform. $$\begin{equation}
    \delta J_{n,\mathcal{B}} = \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}} \, \left( \frac{\partial J_{n,\mathcal{R}}}{\partial j_{\mathcal{B}}} \, \delta j_{\mathcal{B}} \, + \, \frac{\partial J_{n,\mathcal{R}}}{\partial \psi_{\mathcal{B}}} \, \delta \psi_{\mathcal{B}}   \right)
\end{equation}$$

Let write the expression of $J_{n,\mathcal{R}}(x)$ in a generic Gauss
point $x$ and at toroidal plane $p$.

$$\begin{align}
    J_{n,\mathcal{R}}(x,p) &=  \frac{1}{\mu_0} \left( \sum_{k=1}^{N_{vert}} \sum_{j=1}^{N_{dof}} \sum_{l=1}^{N_{tor}} \, \, j_{\mathcal{B}}(k,j,l)  \, k_{size}(k,j) \, H_{1}(k,j, x) \, HZ(l,p) \right) \nonumber \\
    &\cdot  \left( \sum_{k=1}^{N_{vert}} \sum_{j=1}^{N_{dof}} \sum_{l=1}^{N_{tor}} \, \psi_{\mathcal{B}}(k,j,l)  \, k_{size}(k,j) \, H_{1,s}(k,j, x) \, HZ(l,p) \right)  \nonumber \\
    & / \left( \sqrt{R_s^2(x) + Z_s^2(x)} \, F_0 \, R(x) \right) \cdot \text{sgn}(x)
\end{align}$$ We can rewrite the expression as $$\begin{equation}
    J_{n,\mathcal{R}}(x,p) = \frac{1}{\mu_0} {S_{j}(x,p) \cdot S_{\psi}(x,p)} \cdot {D(x)}
\end{equation}$$ where
$D(x) = \frac{\text{sgn}(x)}{\sqrt{R_s^2(x) + Z_s^2(x)} \, F_0 \, R(x)}$.\
Note that $$\begin{equation}
    \frac{\partial S_j(x,p)}{\partial j_{B}(k,j,l)} =  k_{\text{size}}(k,j) \, H_{1}(k,j,x) \, HZ(l,p)
\end{equation}$$

$$\begin{equation}
    \frac{\partial S_{\psi}(x,p)}{\partial \psi_{B}(k,j,l)} = k_{\text{size}}(k,j) \, H_{1,s}(k,j,x) \, HZ(l,p)
\end{equation}$$

And if we define the vectors $\mathbf{{S}}_j$ and $\mathbf{{S}}_\psi$
containing sums in all the Gauss points (first ordered poloidally and
then toroidally): $$\begin{equation}
    \frac{\partial \mathbf{{S}}_j}{\partial j_{B}(k,j,l)} =  k_{\text{size}}(k,j) \, \mathbf{{H}}_{1}(k,j,:) \cdot \mathbf{{HZ}}(l,:)
\end{equation}$$

$$\begin{equation}
    \frac{\partial \mathbf{{S}}_{\psi}}{\partial \psi_{B}(k,j,l)} = k_{\text{size}}(k,j) \, \mathbf{{H}}_{1,s}(k,j,:) \, \mathbf{{HZ}}(l,:)
\end{equation}$$ where both $\mathbf{{H}}_1(k,j)$ and $\mathbf{{HZ}}(l)$
are arrays with $N_{gauss}$ elements and express the basis functions'
values in each of the points. The total number of Bezier dofs in JOREK
is $b = k \cdot j \cdot l$.\
Let first write the expression for a single element of the vector
$\delta J_{n,\mathcal{B}}$.

$$\begin{align}
 \label{eq:deltaJnB}
    \delta J_{n,\mathcal{B}}(k,j,l) &= \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(b,:) \left[ \frac{ \left(k_{\text{size}}(k,j) \, \mathbf{{H}}_{1}(k,j,:) \cdot \mathbf{{HZ}}(l,:) \right) \cdot \mathbf{{S}}_{\psi}(:) } {\mathbf{{D}}} \, \delta j_{\mathcal{B}}(k,j,l) \right] \\ \nonumber
    &+ \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(b,:) \left[\frac{ \left(k_{\text{size}}(k,j) \, \mathbf{{H}}_{1s}(k,j,:) \cdot \mathbf{{HZ}}(l,:) \right) \cdot \mathbf{{S}}_{j}(:) } {\mathbf{{D}}}   \, \delta \psi_{\mathcal{B}}(k,j,l)  \right]
\end{align}$$ The products $\cdot$ in
[\[eq:deltaJnB\]](#eq:deltaJnB){reference-type="eqref"
reference="eq:deltaJnB"} are all element-wise. At last, the matrix
product with conversion matrix
$\underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}$ transfers from real
space to Bezier.\
Let us now group the whole multiplicative terms: $$\begin{equation}
    \delta J_{n,\mathcal{B}}(b) =  \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(b,:) \left\{ 
    \left[  \mathbf{{k}}_{size}(b) \cdot \mathbf{{R}}(b,:) \; \delta j_{\mathcal{B}}(b)    \right]  +  \left[  \mathbf{{k}}_{size}(b) \cdot \mathbf{{R}}_s(b,:) \; \delta \psi_{\mathcal{B}}(b)    \right] \right\}
\end{equation}$$ How would the field equation look if we plug these two
expressions inside the equation? $$\begin{align}
    B_{tan} &= \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \Psi_j + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot Y_{k_i} \right) + \\ \notag
     & \left( \sum_{k_e} \hat{M}^{ey_e}_{i,k_e}  \cdot \sum_{i_e} P^w_{k_e,i_e} \, \sum_j W^w_{ie,j} \, J_{n,{j}}    - \sum_{el_{J_i}} \hat{M}^{e,el}_{i,el_{J_i}} \cdot \sum_j W^J_{i,el_{J_i},j} \, J_{n,{j}}  \right)
\end{align}$$

Let define two new matrices: $$\begin{equation}
 \label{eq:Wmatrix}
    \underline{\underline{\mathbf{\hat{W}}}} = \underline{\underline{\mathbf{M}}}^{eye} \, \left( \underline{\underline{\mathbf{F}}}^T_w \, \underline{\underline{\mathbf{F}}}_w \right)^{-1} \, \underline{\underline{\mathbf{W}}}^w
\end{equation}$$ $$\begin{equation}
\label{eq:Umatrix}
    \underline{\underline{\mathbf{\hat{U}}}} = - \underline{\underline{\mathbf{M}}}^{e,el_J} \, \underline{\underline{\mathbf{W}}}^J
\end{equation}$$ and rewrite the equation as $$\begin{align}
    B^{(n)}_{tan} &= \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \Psi_j + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot Y_{k_i} \right) + \\ \notag
     & \sum_i b_i \left( \sum_{j} \hat{W}_{i,j}  \, J_{n,{j}}    + \sum_{j} \hat{U}_{i,j}  \, J_{n,{j}} \right)
\end{align}$$

$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \Psi^{(n)}_j + \delta \Psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot \left( Y_{k_i}^{(n)} + \delta Y_{k_i}^{(n)} \right) \right) + \\ \notag
     + & \sum_i b_i \left( \sum_{j} \hat{W}_{i,j}  \, \left[ J_{n,{j}}^{(n)} + \delta J_{n,{j}}^{(n)} \right]    + \sum_{j} \hat{U}_{i,j}  \, \left[ J_{n,{j}}^{(n)} + \delta J_{n,{j}}^{(n)} \right] \right)
\end{align}$$

We have expressions for the $J^{(n)}_{n,j}$ and $\delta J^{(n)}_{n,j}$:
$$\begin{equation}
    % J_{n,{j}} =  \Amatu{M}_{\mathcal{B,F}}(j,:) \left\{ 
    % \left[  \Avec{k}_{size}(j) \cdot \Amat{R}(j,:) \; j_{\mathcal{B}}(j)    \right]  +  \left[  \Avec{k}_{size}(j)  \cdot \Amat{R}_s(j,:) \; \psi_{\mathcal{B}}(j)    \right] \right\}
    J_{n,{j}} = \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(j,:) \, J_{n,\mathcal{R}}
\end{equation}$$

$$\begin{equation}
    \delta J_{n,{j}} =  \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(j,:) \left\{ 
    \left[  \mathbf{{k}}_{size}(j) \cdot \mathbf{{R}}(j,:) \; \delta j_{\mathcal{B}}(j)    \right]  +  \left[  \mathbf{{k}}_{size}(j) \cdot \mathbf{{R}}_s(j,:) \; \delta \psi_{\mathcal{B}}(j)    \right] \right\}
\end{equation}$$

For every degree of freedom in JOREK, we can define $$\begin{align}
    c_{\text{zj},j} = \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(j,:) \, \cdot k_{size,j} \cdot \,  \cdot \underline{\mathbf{R}}(:,j) \\
      c_{\psi,j} = \underline{\underline{\mathbf{M}}}_{\mathcal{B,F}}(j,:) \, \cdot k_{size,j} \cdot \,  \cdot \underline{\mathbf{R}}_s(:,j)
\end{align}$$ $$\begin{align}
    J^{(n)}_{n,j} = c_{\text{zj},j} \cdot \text{zj}^{(n)}_j + c_{\psi,j} \cdot \psi_j^{(n)} \\
   \delta J^{(n)}_{n,j} = c_{\text{zj},j} \cdot \delta \text{zj}^{(n)}_j + c_{\psi,j} \cdot \delta \psi^{(n)}_j
\end{align}$$ Leveraging an analogous wall current equation as for (8)
in [@hoelzl2012coupling] applied to $Y_{k_i}$, and observing that we can
sum up
$\underline{\underline{\mathbf{\hat{W}}}} + \underline{\underline{\mathbf{\hat{U}}}} = \underline{\underline{\mathbf{\hat{T}}}}$,
we get to the form: $$\begin{align}
 \label{eq:final Btan equation time-stepping}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{E}_{i,j} \, \delta \Psi^{(n)}_j + \sum_{k_i} \hat{F}^i_{i,k_i} \, Y^{(n)}_{k_i} + \sum_{k_i} \hat{G}_{i,k_i} \, \delta Y^{(n-1)}_{k_{i}} + \sum_j \hat{H}_{i,j} \, \Psi^{(n)}_j + \sum_j \hat{J}_{i,j} \delta \Psi^{(n-1)}_j
    \right] \\ \notag
    & + \sum_i b_i \left[  \sum_j \hat{T}_{i,j} \left[  c_{\text{zj}} \cdot \left( \text{zj}_j^{(n)} + \delta \text{zj}_j^{(n)}  \right) + c_{\psi} \cdot \left( \Psi_j^{(n)} + \delta \Psi_j^{(n)}  \right)    \right]    \right]
\end{align}$$

Equation
[\[eq:final Btan equation time-stepping\]](#eq:final Btan equation time-stepping){reference-type="eqref"
reference="eq:final Btan equation time-stepping"} contains the
expression we need to provide to the BCs.\
For reduced MHD, let us now implement them in the JOREK current
definition equation $$\begin{equation}
    \int dV \, \frac{j^*_l}{R^2} \, j + \int dV \, \frac{1}{R^2} \, \nabla j^*_l \cdot \nabla \Psi - \oint dA \, \frac{j^*_l}{R} \, B_{tan} = 0
\end{equation}$$

Let's now plug in expression
[\[eq:final Btan equation time-stepping\]](#eq:final Btan equation time-stepping){reference-type="eqref"
reference="eq:final Btan equation time-stepping"}:

$$\begin{align}
 \label{eq:curr def eq w halos}
    & \sum_{i_{elem}} \int \frac{dV}{R^2} \, \left( j^*_l \, \delta j^{(n)} + \nabla j^*_l \cdot \delta \Psi^{(n)} \right) - \sum_{i_{bnd}} \oint dA \, \frac{j^*_l}{R} \sum_i b_i \left( \sum_j \hat{E}_{i,j} \, \delta\Psi^{(n)}_j + \color{blue} \sum_j \hat{T}_{i,j} \, c_{\Psi} \, \delta\Psi^{(n)}_j \color{black} \right) \\ \notag 
    & \color{blue} -\sum_{i_{bnd}} \oint dA \, \frac{j^*_l}{R} \sum_i b_i  \sum_j \hat{T}_{i,j} \, c_{\text{zj}} \, \delta\text{zj}^{(n)}_j \color{black} = \\ \notag
    & -\sum_{i_{elem}} \int \frac{dV}{R^2} \left( j^*_l \, j^{(n)} + \nabla j^{(n)}_l \cdot \nabla \Psi^{(n)} \right) + \\ \notag
    & + \sum_{i_{bnd}} \oint dA \frac{j^*_l}{R} \sum_i b_i \left[ \sum_{k_i} \left( \hat{F}^i_{i,k_i} \, Y^{(n)}_{k_i} + \hat{G}_{i,k_i} \, \delta Y^{(n-1)}_{k_{i}}  \right)  +  \sum_j \left( \hat{H}_{i,j} \, \Psi^{(n)}_j + \hat{J}_{i,j} \delta \Psi^{(n-1)}_j + \color{red} \hat{T}_{i,j} J^{(n)}_{n,j} \right) \color{black}  \right]
\end{align}$$

From equation
[\[eq:curr def eq w halos\]](#eq:curr def eq w halos){reference-type="eqref"
reference="eq:curr def eq w halos"} we notice two new terms on the LHS
and one new term on the RHS of the equation.\
The geometric factors in
$$J_n (x) = \frac{j(x) \, \psi_s (x) }{\sqrt{R_s^2(x) + Z_s^2(x)} \, F_0 \, R(x) } \, \text{sgn}(x)$$
including the derivative $\frac{\partial \psi}{\partial s}$ are to be
included in the boundary integral as `common_prefactor` for each of
these terms.

## Wall current time stepping

We derive here the time stepping form for the extended wall current
equation $$\begin{equation}
    \mathbf{{L}}^{*}_{w,ii} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = \mathbf{{C}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot} + \mathbf{{D}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,eq} + \mathbf{{E}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,w} + \mathbf{{G}}_i^* \, \mathbf{{i}}_{el,w} + \mathbf{{K}}_{w} \, \mathbf{{v}}
\end{equation}$$ To keep consistency with [@hoelzl2012coupling], we
change the sign of the $\mathbf{{C}}_i^*$ term and pass to a form of
type $$\begin{equation}
 \label{eq:time_dep_Yki_eq}
    \dot{Y_{k_i}} + \sum_j \hat{M}^{y_i,e}_{k_i,j} \dot{\psi_j} - \sum_{l_{eq}} \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \dot{i}_{el,l_{eq}} - \sum_{l_w} \hat{E}^{y_i, w_e}_{k_i,l_w} \dot{i}_{el,l_w} = - \frac{\eta_w}{d_w} \hat{M}^{y_i,y_i}_{k_i,k_i} Y_{k_i}  + \sum_{l_w} \hat{G}^{y_i, w_e}_{k_i,l_w} i_{el,l_w}
\end{equation}$$

Let now apply the temporal discretization to this form
$\partial \mathbf{A(u)} / \partial t = \mathbf{{B}}(\mathbf{u},t)$,
which yields the linearized form $$\begin{equation}
\left[ (1+\xi) \left( \frac{\partial \mathbf{A}}{\partial \mathbf{u}}\right)^{(n)} - \Delta t \theta \left( \frac{\partial \mathbf{B}}{\partial \mathbf{u}}\right)^{(n)}   \right] \delta \mathbf{u}^{(n)} = \Delta t \mathbf{B}^{(n)} + \xi \left( \frac{\partial \mathbf{A}}{\partial \mathbf{u}}\right)^{(n)} \delta \mathbf{u}^{(n-1)}
\end{equation}$$ Applying it to our equation
[\[eq:time_dep_Yki_eq\]](#eq:time_dep_Yki_eq){reference-type="eqref"
reference="eq:time_dep_Yki_eq"} $$\begin{multline}
     (1+\xi) \left[ \delta Y^{(n)}_{k_i} + \sum_j \hat{M}^{y_i,e}_{k_i,j} \delta\psi^{(n)}_j- \sum_{l_{eq}}  \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \delta i^{(n)}_{el,l_{eq}} - \sum_{l_w} \hat{E}^{y_i,w_e}_{k_i,l_w} \delta i^{(n)}_{el,l_{w}} \right] \\
     - \Delta t \theta \left[ - \hat{M}^{y_i,y_i}_{k_i,k_i} \, \delta Y_{k_i}^{(n)} + \sum_{l_w} \hat{G}^{y_i,w_e}_{k_i,l_w} \, \delta i^{(n)}_{el,l_w}  \right] = \\
     = \Delta t \left[ - \hat{M}^{y_i,y_i}_{k_i,k_i} \, Y_{k_i}^{(n)} + \sum_{l_w} \hat{G}^{y_i,w_e}_{k_i,l_w} \, i^{(n)}_{el,l_w}  \right] \\
     + \xi \left( \delta Y^{(n-1)}_{k_i}  + \sum_j \hat{M}^{y_i,e}_{k_i,j} \delta\psi^{(n-1)}_j- \sum_{l_{eq}}  \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \delta i^{(n-1)}_{el,l_{eq}} - \sum_{l_w} \hat{E}^{y_i,w_e}_{k_i,l_w} \delta i^{(n-1)}_{el,l_{w}}\right)
\end{multline}$$ $$\begin{multline}
    \left( 1+\xi+\Delta t \theta \hat{M}^{y_i,y_i}_{k_i,k_i} \right) \, \delta Y^{(n)}_{k_i}  =  (1+\xi) \left[ - \sum_j \hat{M}^{y_i,e}_{k_i,j}  \delta\psi^{(n)}_j  + \sum_{l_{eq}}  \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \delta i^{(n)}_{el,l_{eq}} + \sum_{l_w} \hat{E}^{y_i,w_e}_{k_i,l_w} \delta i^{(n)}_{el,l_{w}}\right] +\\
    + \Delta t \theta \sum_{l_w} \hat{G}^{y_i,w_e}_{k_i,l_w} \, \delta i^{(n)}_{el,l_w} - \Delta t \, \hat{M}^{y_i,y_i}_{k_i,k_i} \, Y_{k_i}^{(n)} + \Delta t \,\sum_{l_w} \hat{G}^{y_i,w_e}_{k_i,l_w} \, i^{(n)}_{el,l_w} + \\
    + \xi \left( \delta Y^{(n-1)}_{k_i}  + \sum_j \hat{M}^{y_i,e}_{k_i,j} \delta\psi^{(n-1)}_j- \sum_{l_{eq}}  \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \delta i^{(n-1)}_{el,l_{eq}} - \sum_{l_w} \hat{E}^{y_i,w_e}_{k_i,l_w} \delta i^{(n-1)}_{el,l_{w}}\right)
\end{multline}$$ The equation for $\delta Y^{(n)}_{k_i}$ may be written
as: $$\begin{multline}
    \delta Y^{(n)}_{k_i} = \sum_j \hat{A}_{k_i,j} \, \delta \psi^{(n)}_j + \sum_{l_{eq}} \widehat{JE}_{k_i,l_{eq}} \, \delta i^{(n)}_{el,l_{eq}} +  \sum_{l_{w}} \widehat{WE}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \sum_{l_w} \widehat{WI}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \hat{B}_{k_i,k_i} \, Y^{(n)}_{k_i} + \\
    + \sum_{l_w} \widehat{WH}_{k_i,l_{w}} \,  i^{(n)}_{el,w} 
    + \hat{C}_{k_i,k_i} \, \delta Y^{(n-1)}_{k_i} + \sum_j \hat{D}_{k_i,j} \, \delta \psi^{(n-1)}_j + \sum_{l_{eq}} \widehat{JP}_{k_i,l_{eq}} \, \delta i^{(n-1)}_{el,l_{eq}} + \sum_{l_w} \widehat{WP}_{k_i,l_w} \delta i^{(n-1)}_{el,l_w}
\end{multline}$$ where $$\begin{align}
    \hat{S}_{k_i,k_i} &= 1+\xi+\Delta t \theta \, \frac{\eta_w}{d_w} \hat{M}^{y_i,y_i}_{k_i,k_i} \\
    \hat{A}_{k_i,j} &= -(1+\xi) \, \hat{M}^{y_i,e}_{k_i,j} / \hat{S}_{k_i,k_i} \\
    \widehat{JE}_{k_i,l_{eq}} &= (1+\xi) \, \hat{D}^{y_i, eq_e} _{k_i,l_{eq}} / \hat{S}_{k_i,k_i} \\
    \widehat{WE}_{k_i,l_{w}} &= (1+\xi) \, \hat{E}^{y_i, w_e}_{k_i,l_{w}} / \hat{S}_{k_i,k_i} \\
    \widehat{WI}_{k_i,l_{w}} &= \Delta t \theta \, \hat{G}^{y_i,w_e}_{k_i,l_w} / \hat{S}_{k_i,k_i} \\
    \hat{B}_{k_i,k_i} &= - \Delta t \frac{\eta_w}{d_w} \hat{M}^{y_i,y_i}_{k_i,k_i} / \hat{S}_{k_i,k_i} \\
    \widehat{WH}_{k_i,l_w} &= \Delta t \,  \hat{G}^{y_i,w_e}_{k_i,l_w} / \hat{S}_{k_i,k_i}  \\
    \hat{C}_{k_i,k_i} &= \xi / \hat{S}_{k_i,k_i} \\
    \hat{D}_{k_i,j} &= \xi \hat{M}^{y_i,e}_{k_i,j} / \hat{S}_{k_i,k_i} \\
    \widehat{JP}_{k_i,l_{eq}} &= - \xi \, \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} / \hat{S}_{k_i,k_i} \\
    \widehat{WP}_{k_i,l_{w}} &= - \xi \, \hat{E}^{y_i, w_e}_{k_i,l_{w}} / \hat{S}_{k_i,k_i}
\end{align}$$

Let us now replace it in the expression for the magnetic field time
stepping

$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left( \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \Psi^{(n)}_j + \delta \Psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot \left( Y_{k_i}^{(n)} + \delta Y_{k_i}^{(n)} \right) \right) + \\ \notag
     + & \sum_i b_i \left( \sum_{j} \hat{W}_{i,j}  \, \left[ J_{n,{j}}^{(n)} + \delta J_{n,{j}}^{(n)} \right]    + \sum_{j} \hat{U}_{i,j}  \, \left[ J_{n,{j}}^{(n)} + \delta J_{n,{j}}^{(n)} \right] \right)
\end{align}$$ which we can make more compact here: $$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \psi^{(n)}_j + \delta \psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot \left( Y_{k_i}^{(n)} + \delta Y_{k_i}^{(n)} \right) \right] + \\ \notag
     + & \sum_i b_i \left[ \sum_{j} \hat{T}_{i,j}  \left( J^{(n)}_{n,j}  + c_{\text{zj}} \delta \text{zj}^{(n)}_j + c_{\psi} \delta \psi^{n}_j \right) \right]
\end{align}$$ Let replace the expression for $\delta Y^{(n)}_{k_i}$:
$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \psi^{(n)}_j + \delta \psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot  Y_{k_i}^{(n)} \right] + \\ \notag
     &  +\sum_i b_i \sum_{k_i} \hat{M}^{ey_i}_{i,k_i}  \left[ \sum_j \hat{A}_{k_i,j} \, \delta \psi^{(n)}_j + \sum_{l_{eq}} \widehat{JE}_{k_i,l_{eq}} \, \delta i^{(n)}_{el,l_{eq}} +  \sum_{l_{w}} \widehat{WE}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \sum_{l_w} \widehat{WI}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \right. \\ \notag 
     + &\hat{B}_{k_i,k_i} \, Y^{(n)}_{k_i} + \sum_{l_w} \widehat{WH}_{k_i,l_{w}} \,  i^{(n)}_{el,w} 
    + \hat{C}_{k_i,k_i} \, \delta Y^{(n-1)}_{k_i} + \sum_j \hat{D}_{k_i,j} \, \delta \psi^{(n-1)}_j + \sum_{l_{eq}} \widehat{JP}_{k_i,l_{eq}} \, \delta i^{(n-1)}_{el,l_{eq}} + \\ \notag 
    +& \left. \sum_{l_w} \widehat{WP}_{k_i,l_w} \delta i^{(n-1)}_{el,l_w}  \right] + \\ \notag 
   +&  \sum_i b_i \left[ \sum_{j} \hat{T}_{i,j}  \left( J^{(n)}_{n,j}  + c_{\text{zj},j} \, \delta \text{zj}^{(n)}_j + c_{\psi,j} \, \delta \psi^{n}_j \right) \right]
\end{align}$$

In order to eliminate $J_n$ as a variable, we write it as:
$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \psi^{(n)}_j + \delta \psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot  Y_{k_i}^{(n)} \right] + \\ \notag
     &  +\sum_i b_i \sum_{k_i} \hat{M}^{ey_i}_{i,k_i}  \left[ \sum_j \hat{A}_{k_i,j} \, \delta \psi^{(n)}_j + \sum_{l_{eq}} \widehat{JE}_{k_i,l_{eq}} \, \delta i^{(n)}_{el,l_{eq}} +  \sum_{l_{w}} \widehat{WE}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \sum_{l_w} \widehat{WI}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \right. \\ \notag 
     + &\hat{B}_{k_i,k_i} \, Y^{(n)}_{k_i} + \sum_{l_w} \widehat{WH}_{k_i,l_{w}} \,  i^{(n)}_{el,w} 
    + \hat{C}_{k_i,k_i} \, \delta Y^{(n-1)}_{k_i} + \sum_j \hat{D}_{k_i,j} \, \delta \psi^{(n-1)}_j + \sum_{l_{eq}} \widehat{JP}_{k_i,l_{eq}} \, \delta i^{(n-1)}_{el,l_{eq}} + \\ \notag 
    +& \left. \sum_{l_w} \widehat{WP}_{k_i,l_w} \delta i^{(n-1)}_{el,l_w}  \right] + \\ \notag 
   +&  \sum_i b_i \left[ \sum_{j} \hat{T}_{i,j}  \left(  c_{\text{zj},j} \, \delta \text{zj}^{(n)}_j + c_{\psi,j} \, \delta \psi^{n}_j \right) \right] + \\ \notag 
   +&  \sum_i b_i \left[ \sum_{l_w} \hat{M}^{e,el_w}_{i,l_w} \, i^{(n)}_{el,l_w} - \sum_{l_{eq}} \hat{M}^{e,el_J}_{i,l_{eq}} \, i^{(n)}_{el,l_{eq}} \right]
\end{align}$$ Let us now express the electrode current variations as
follows: $$\begin{equation}
    \delta i^{(n)}_{el,l_{eq}} = \sum_j W^J_{l_{eq},j} \left( c_{\text{zj},j} \, \delta \text{zj}^{(n)}_j + c_{\psi,j} \, \delta \psi^{(n)}_j \right)
\end{equation}$$ $$\begin{equation}
    \delta i^{(n)}_{el,l_{w}} = -\sum_j W^w_{l_{w},j} \left( c_{\text{zj},j} \, \delta \text{zj}^{(n)}_j + c_{\psi,j} \, \delta \psi^{(n)}_j \right)
\end{equation}$$ Let us now group the terms that multiply each variable:
$$\begin{align}
    & \delta \psi_j^{(n)}: \; \color{red} \hat{E_H}_{i,j} = \\ \notag   & \color{red} = \hat{M}^{e,e}_{i,j}+\sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \left[  \hat{A}_{k_i,j} +  \left(  \sum_{l_{eq}} \widehat{JE}_{k_i,l_{eq}} W^J_{l_{eq},j}    - \sum_{l_w} \left(  \widehat{WE}_{k_i,l_w} + \widehat{WI}_{k_i,l_w} \right) W^w_{l_w,j}    \right) c_{\psi,j}  \right]  + \hat{T}_{i,j} \, c_{\psi,j} \\ 
    & \delta \text{zj}^{(n)} : \; \color{red}  \hat{Z}_{i,j} = \\ \notag & \color{red} = \sum_{k_i} \hat{M}^{e,y_i}_{i,k_i}  \left[ \left(  \sum_{l_{eq}} \widehat{JE}_{k_i,l_{eq}} W^J_{l_{eq},j}  - \sum_{l_w} \left(  \widehat{WE}_{k_i,l_w} + \widehat{WI}_{k_i,l_w} \right) W^w_{l_w,j}    \right) c_{\text{zj},j}   \right]  + \hat{T}_{i,j} \, c_{\text{zj},j} \\
    & Y^{(n)}_{k_i} : \; \color{blue} \hat{F}_{i,k_i} = \hat{M}^{e,y_i}_{i,k_i} \left( 1+\hat{B}_{k_i,k_i}\right) \\
    & \delta Y^{(n-1)}_{k_i} : \; \color{blue} \hat{G}_{i,k_i} = \hat{M}^{e,y_i}_{i,k_i} \, \hat{C}_{k_i,k_i} \\
    & \psi^{(n)}_j : \; \color{green} \hat{H}_{i,j} = \hat{M}^{e,e}_{i,j} \\
    & \delta \psi^{(n-1)}_j : \; \color{green} \hat{J}_{i,j} = \sum_{k_i} \hat{M}^{e,y_i}_{i,ki} \, \hat{D}_{k_i,j} \\
    & i^{(n)}_{el,w} : \;  \color{orange} \hat{P_H}_{i,l_w} = \sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \widehat{WH}_{k_i,l_w} + \hat{M}^{e,el_w}_{i,l_w}  \\
    & \delta i^{(n-1)}_{el,w} : \;  \color{orange} \hat{Q_H}_{i,l_w} = \sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \widehat{WP}_{k_i,l_w} \\ 
    & \delta i^{(n)}_{el,eq} : \; \color{cyan} \hat{X_H}_{i,l_{eq}} = - \hat{M}^{e,el_J}_{i,l_{eq}} \\
    & \delta i^{(n-1)}_{el,eq} : \; \color{cyan} \hat{V_H}_{i,l_{eq}} = \sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \widehat{JP}_{k_i,l_{eq}} \\
\end{align}$$ We write the final expression for the magnetic field
$B^{(n+1)}_{tan}$: $$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[  \sum_j \hat{E_H}_{i,j} \, \delta \psi^{(n)}_j   +  \sum_j \hat{Z}_{i,j} \, \delta \text{zj}^{(n)}_j    + \sum_{k_i} \hat{F}_{i,k_i} \, Y^{(n)}_{k_i}  +  \sum_{k_i} \hat{G}_{i,k_i} \, Y^{(n-1)}_{k_i}  \right.   \\ \notag
   & \left. \sum_j \hat{H}_{i,j} \, \psi^{(n)}_{j}  + \sum_j \hat{J}_{i,j} \, \delta \psi^{(n-1)}_j  + \sum_{l_w} \hat{P_H}_{i,l_w} \, i^{(n)}_{el,w} +  \sum_{l_w} \hat{Q_H}_{i,l_w} \, i^{(n-1)}_{el,w} + \right.  \\ \notag
   & \left. \sum_{l_{eq}} \hat{X_H}_{i,l_{eq}} \, i^{(n)}_{el,eq} +  \sum_{l_{eq}} \hat{V_H}_{i,l_{eq}} \, i^{(n-1)}_{el,eq} \right]
\end{align}$$ Boundary condition for current definition equation:
$$\begin{equation}
    \int dV \, \frac{j^*_l}{R^2} \, j + \int dV \, \frac{1}{R^2} \, \nabla j^*_l \cdot \nabla \Psi - \oint dA \, \frac{j^*_l}{R} \, B_{tan} = 0
\end{equation}$$ $$\begin{align}
 \label{eq:curr def eq w halos tot}
    & \sum_{i_{elem}} \int \frac{dV}{R^2} \, \left( j^*_l \, \delta j^{(n)} + \nabla j^*_l \cdot \delta \Psi^{(n)} \right) - \sum_{i_{bnd}} \oint dA \, \frac{j^*_l}{R} \sum_i b_i \left( \sum_j \hat{E_H}_{i,j} \, \delta\Psi^{(n)}_j    \right) \\ \notag 
    &  - \sum_{i_{bnd}} \oint dA \, \frac{j^*_l}{R} \sum_i b_i \left( \sum_j \hat{Z}_{i,j} \, \delta\text{zj}^{(n)}_j    \right) = \\ \notag
    & -\sum_{i_{elem}} \int \frac{dV}{R^2} \left( j^*_l \, j^{(n)} + \nabla j^{(n)}_l \cdot \nabla \Psi^{(n)} \right) + \\ \notag
    & + \sum_{i_{bnd}} \oint dA \frac{j^*_l}{R} \sum_i b_i \left[ \sum_{k_i} \left( \hat{F}_{i,k_i} \, Y^{(n)}_{k_i} + \hat{G}_{i,k_i} \, \delta Y^{(n-1)}_{k_{i}}  \right)  +  \sum_j \left( \hat{H}_{i,j} \, \Psi^{(n)}_j + \hat{J}_{i,j} \delta \Psi^{(n-1)}_j \right) \right] \\ \notag 
    & + \sum_{i_{bnd}} \oint dA \frac{j^*_l}{R} \sum_i b_i \left[  \sum_{l_w} \left( \hat{P_H}_{i,l_w} \, i^{(n)}_{el,w} + \hat{Q_H}_{i,l_w} \, \delta i^{(n-1)}_{el,w}    \right) + \sum_{l_{eq}} \left(  \hat{X_H}_{i,l_{eq}} \, i^{(n)}_{el,eq} + \hat{V_H}_{i,l_{eq}} \, i^{(n-1)}_{el,eq} \right) \right]
\end{align}$$

## Implicit stepping with only static correction {#implicit-stepping-with-only-static-correction .unnumbered}

Given the equation
[\[eq:time_dep_Yki_eq\]](#eq:time_dep_Yki_eq){reference-type="eqref"
reference="eq:time_dep_Yki_eq"} $$\begin{equation}
    \dot{Y_{k_i}} + \sum_j \hat{M}^{y_i,e}_{k_i,j} \dot{\psi_j} - \sum_{l_{eq}} \hat{D}^{y_i, eq_e}_{k_i,l_{eq}} \dot{i}_{el,l_{eq}} - \sum_{l_w} \hat{E}^{y_i, w_e}_{k_i,l_w} \dot{i}_{el,l_w} = - \hat{M}^{y_i,y_i}_{k_i,k_i} Y_{k_i}  + \sum_{l_w} \hat{G}^{y_i, w_e}_{k_i,l_w} i_{el,l_w}
\end{equation}$$ we consider only the form with the resistive coupling:
$$\begin{equation}
 \label{eq:static_Yki_eq}
    \dot{Y_{k_i}} + \sum_j \hat{M}^{y_i,e}_{k_i,j} \dot{\psi_j}  = - \hat{M}^{y_i,y_i}_{k_i,k_i} Y_{k_i}  + \sum_{l_w} \hat{G}^{y_i, w_e}_{k_i,l_w} i_{el,l_w}
\end{equation}$$ Which, by expanding with the implicit time stepping,
yields: $$\begin{multline}
    \delta Y^{(n)}_{k_i} = \sum_j \hat{A}_{k_i,j} \, \delta \psi^{(n)}_j + \sum_{l_w} \widehat{WI}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \hat{B}_{k_i,k_i} \, Y^{(n)}_{k_i} + \\
    + \sum_{l_w} \widehat{WH}_{k_i,l_{w}} \,  i^{(n)}_{el,w} 
    + \hat{C}_{k_i,k_i} \, \delta Y^{(n-1)}_{k_i} + \sum_j \hat{D}_{k_i,j} \, \delta \psi^{(n-1)}_j
\end{multline}$$ which we replace in the magnetic field expression
$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \psi^{(n)}_j + \delta \psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot \left( Y_{k_i}^{(n)} + \delta Y_{k_i}^{(n)} \right) \right] + \\ \notag
     + & \sum_i b_i \left[ \sum_{j} \hat{T}_{i,j}  \left( J^{(n)}_{n,j}  + c_{\text{zj}} \delta \text{zj}^{(n)}_j + c_{\psi} \delta \psi^{n}_j \right) \right]
\end{align}$$

$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[ \sum_j \hat{M}^{e,e}_{i,j} \cdot \left( \psi^{(n)}_j + \delta \psi^{(n)}_j \right) + \sum_{k_i} \hat{M}^{ey_i}_{i,k_i} \cdot  Y_{k_i}^{(n)} \right] + \\ \notag
     &  +\sum_i b_i \sum_{k_i} \hat{M}^{ey_i}_{i,k_i}  \left[ \sum_j \hat{A}_{k_i,j} \, \delta \psi^{(n)}_j +  \sum_{l_w} \widehat{WI}_{k_i,l_{w}} \, \delta i^{(n)}_{el,l_w} + \right. \\ \notag 
     + &\hat{B}_{k_i,k_i} \, Y^{(n)}_{k_i} + \sum_{l_w} \widehat{WH}_{k_i,l_{w}} \,  i^{(n)}_{el,w} 
    + \hat{C}_{k_i,k_i} \, \delta Y^{(n-1)}_{k_i} + \sum_j \hat{D}_{k_i,j} \, \delta \psi^{(n-1)}_j  + \\ \notag 
   +&  \sum_i b_i \left[ \sum_{j} \hat{T}_{i,j}  \left(  c_{\text{zj},j} \, \delta \text{zj}^{(n)}_j + c_{\psi,j} \, \delta \psi^{n}_j \right) \right] + \\ \notag 
   +&  \sum_i b_i \left[ \sum_{l_w} \hat{M}^{e,el_w}_{i,l_w} \, i^{(n)}_{el,l_w} - \sum_{l_{eq}} \hat{M}^{e,el_J}_{i,l_{eq}} \, i^{(n)}_{el,l_{eq}} \right]
\end{align}$$ The new matrices become $$\begin{align}
    & \delta \psi_j^{(n)}: \; \color{red} \hat{E_H}_{i,j} = \\ \notag   & \color{red} = \hat{M}^{e,e}_{i,j}+\sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \left[  \hat{A}_{k_i,j} +  \left(   - \sum_{l_w} \left(   \widehat{WI}_{k_i,l_w} \right) W^w_{l_w,j}    \right) c_{\psi,j}  \right]  + \hat{T}_{i,j} \, c_{\psi,j} \\ 
    & \delta \text{zj}^{(n)} : \; \color{red} \hat{Z}_{i,j} = \\ \notag & \color{red} = \hat{M}^{e,y_i}_{i,k_i}  \left[ \left(  - \sum_{l_w} \left( \widehat{WI}_{k_i,l_w} \right) W^w_{l_w,j}    \right) c_{\text{zj},j}   \right]  + \hat{T}_{i,j} \, c_{\text{zj},j} \\
    & Y^{(n)}_{k_i} : \; \color{blue} \hat{F}_{i,k_i} = \hat{M}^{e,y_i}_{i,k_i} \left( 1+\hat{B}_{k_i,k_i}\right) \\
    & \delta Y^{(n-1)}_{k_i} : \; \color{blue} \hat{G}_{i,k_i} = \hat{M}^{e,y_i}_{i,k_i} \, \hat{C}_{k_i,k_i} \\
    & \psi^{(n)}_j : \; \color{green} \hat{H}_{i,j} = \hat{M}^{e,e}_{i,j} \\
    & \delta \psi^{(n-1)}_j : \; \color{green} \hat{J}_{i,j} = \sum_{k_i} \hat{M}^{e,y_i}_{i,ki} \, \hat{D}_{k_i,j} \\
    & i^{(n)}_{el,w} : \;  \color{orange} \hat{P_H}_{i,l_w} = \sum_{k_i} \hat{M}^{e,y_i}_{i,k_i} \, \widehat{WH}_{k_i,l_w} + \hat{M}^{e,el_w}_{i,l_w}  \\
    & i^{(n)}_{el,eq} : \; \color{cyan} \hat{X_H}_{i,l_{eq}} = - \hat{M}^{e,el_J}_{i,l_{eq}} \\
\end{align}$$

$$\begin{align}
    B^{(n+1)}_{tan} &= \sum_i b_i \left[  \sum_j \hat{E_H}_{i,j} \, \delta \psi^{(n)}_j   +  \sum_j \hat{Z}_{i,j} \, \delta \text{zj}^{(n)}_j    + \sum_{k_i} \hat{F}_{i,k_i} \, Y^{(n)}_{k_i}  +  \sum_{k_i} \hat{G}_{i,k_i} \, Y^{(n-1)}_{k_i}  \right.   \\ \notag
   & \left. \sum_j \hat{H}_{i,j} \, \psi^{(n)}_{j}  + \sum_j \hat{J}_{i,j} \, \delta \psi^{(n-1)}_j  + \sum_{l_w} \hat{P_H}_{i,l_w} \, i^{(n)}_{el,w}  + \right.  \\ \notag
   & \left. \sum_{l_{eq}} \hat{X_H}_{i,l_{eq}} \, i^{(n)}_{el,eq} \right]
\end{align}$$

# Halo coupling matrices

Here we list the full set of matrices that are being calculated for the
halo coupling. See table
[1](#tab:resume_final_equation){reference-type="ref"
reference="tab:resume_final_equation"}, equation
[\[eq:Solenoidal currents evolution eq.\]](#eq:Solenoidal currents evolution eq.){reference-type="eqref"
reference="eq:Solenoidal currents evolution eq."} and equation
[\[eq:Btheta final form\]](#eq:Btheta final form){reference-type="eqref"
reference="eq:Btheta final form"} for reference.

:::: center
::: tabular
\| P3em \| P3cm\| P4cm \| **Matrix** & **Usage** & **Definition**\
$\underline{\underline{\mathbf{D_{ii}}}}$ & Eigenvalue matrix for wall
currents in the solenoidal sub-space &\
$\underline{\underline{\mathbf{M_{yi,e}}}}$ & $\mathbf{A_{tan}}$
boundary dofs $\rightarrow$ $Y_i$ solenoidal currents &
$\mu_0 \, \underline{\underline{\mathbf{S^T_{ii}}}} \, \underline{\underline{\mathbf{C^*_i}}} \, / \underline{\underline{\mathbf{D_{ii}}}}$
(table [1](#tab:resume_final_equation){reference-type="ref"
reference="tab:resume_final_equation"})\
$\underline{\underline{\mathbf{M_{e,yi}}}}$ & $Y_i$ solenoidal currents
$\rightarrow$ $\mathbf{B_{tan}}$ boundary dofs & Addend 1 in
[\[eq:Btheta final form\]](#eq:Btheta final form){reference-type="eqref"
reference="eq:Btheta final form"}\
$\underline{\underline{\mathbf{M_{e,ye}}}}$ & $I_{we}$ source/sink
currents $\rightarrow$ $\mathbf{B_{tan}}$ boundary dofs & Addend 2 in
[\[eq:Btheta final form\]](#eq:Btheta final form){reference-type="eqref"
reference="eq:Btheta final form"}\
$\underline{\underline{\mathbf{M_{e,e}}}}$ & $\mathbf{A_{tan}}$ boundary
dofs $\rightarrow$ $\mathbf{B_{tan}}$ boundary dofs & Addend 3 in
[\[eq:Btheta final form\]](#eq:Btheta final form){reference-type="eqref"
reference="eq:Btheta final form"}\
$\underline{\underline{\mathbf{M_{eq,el}}}}$ & ${I_{eq,el}}$ shell
electrode currents $\rightarrow$ $\mathbf{B_{tan}}$ boundary dofs &
Addend 4 in
[\[eq:Btheta final form\]](#eq:Btheta final form){reference-type="eqref"
reference="eq:Btheta final form"}\
$\underline{\underline{\mathbf{W^w}}}$ &
$J_{n,bez} \rightarrow I^w_{elec}$ conversion & Equation
[\[eq:JntoIelec\]](#eq:JntoIelec){reference-type="eqref"
reference="eq:JntoIelec"}\
$\underline{\underline{\mathbf{W^J}}}$ &
$J_{n,bez} \rightarrow I^J_{elec}$ conversion & Equation
[\[eq:JntoIelec\]](#eq:JntoIelec){reference-type="eqref"
reference="eq:JntoIelec"}\
$\underline{\underline{\mathbf{\hat{W}}}}$ &
$J_{n,bez} \rightarrow \mathbf{B_{tan}}$ conversion for wall electrodes
&
$\underline{\underline{\mathbf{M_{e,ye}}}} \, \left( \underline{\underline{\mathbf{F}}}^T_w \, \underline{\underline{\mathbf{F}}}_w \right)^{-1} \, \underline{\underline{\mathbf{W}}}^w$
(equation [\[eq:Wmatrix\]](#eq:Wmatrix){reference-type="eqref"
reference="eq:Wmatrix"}\
$\underline{\underline{\mathbf{\hat{U}}}}$ &
$J_{n,bez} \rightarrow \mathbf{B_{tan}}$ conversion for shell electrodes
&
$- \underline{\underline{\mathbf{M_{eq,el}}}} \, \underline{\underline{\mathbf{W}}}^J$
(equation [\[eq:Umatrix\]](#eq:Umatrix){reference-type="eqref"
reference="eq:Umatrix"}\
$\underline{\underline{\mathbf{\hat{T}}}}$ &
$J_{n,bez} \rightarrow \mathbf{B_{tan}}$ conversion &
$\underline{\underline{\mathbf{\hat{W} + \hat{U}}}}$ (equation
[\[eq:curr def eq w halos\]](#eq:curr def eq w halos){reference-type="eqref"
reference="eq:curr def eq w halos"})\
:::
::::

$$\begin{multline}
\label{eq:Btheta final form}
     \underline{\mathbf{B_{\theta}}} = \color{red} \underbrace{  \cancel{\left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}} \, \right) \underline{\underline{\mathbf{K_w}}} \, \underline{\underline{\mathbf{S_{w,ii}}}} \right]}} _{\underline{\underline{\mathbf{M_{e,yi}}}}}  \, \underline{\mathbf{Y_i}} + \\
    \color{blue} \underbrace{\left[ \left( \underline{\underline{\mathbf{B_{\theta,w} }}} \, - \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}   \right) \underline{\underline{\mathbf{F_w}}} \right]}_{\underline{\underline{\mathbf{M_{e,ye}}}}} \, \underline{\mathbf{i_{w,e}}}  \color{black} +   \underbrace{\left( \underline{\underline{\mathbf{B_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{H_i}}} \right)}_{\underline{\underline{\mathbf{M_{ee}}}}} \underline{\mathbf{a_{tot}}}  - \\
    - \color{teal} \underbrace{\left[ \,  \underline{\underline{\mathbf{B_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] }_{\underline{\underline{\mathbf{M_{eq,el}}}}}\underline{\mathbf{i_{el,eq}}}
\end{multline}$$

$$\begin{equation}
 \label{eq:Solenoidal currents evolution eq.}
    \, \underline{\mathbf{\dot{Y_i}}} = -  \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\mathbf{{Y_i}}} + \mu_0 \, \underline{\underline{\mathbf{M_{yi,e}}}} \, {\underline{\mathbf{\dot{a}_{tot}}}}
\end{equation}$$

# Testing the new electrode matrices

The objective of this section is to design a test to validate the new
field matrices (due to the halo currents) against a reference.\
**N.B. Reduced MHD is not fully consistent due to the approximation
$A_{tot} \approx \psi$. If we aim at providing a reference magnetic
field via Biot-Savart, the latter will have all the information for
$\mathbf{A}$. RMHD JOREK will intrinsically carry an error because of
this.**

## Test concept 1 {#test-concept-1 .unnumbered}

Since our true objective is to test the halo field matrices, we design a
test in CARIDDI where:

- $\mathbf{A}_{tot} = \mathbf{A}_w + \mathbf{A}_{pl}$ can be calculated
  via the CARIDDI inductance matrices as: $$\begin{equation}
          \mathbf{A}_{tot} = \underline{\underline{\mathbf{M_{eq,w}}}} \; \mathbf{I}_w + \underline{\underline{\mathbf{L_{eq}}}} \; \mathbf{I}_{eq}
  \end{equation}$$

- Since we are decomposing
  $\mathbf{{I}} = \mathbf{{K}} \, \mathbf{{I}}_i + \mathbf{{F}} \, \mathbf{{I}}_e$
  on both shell and wall, we impose that the internal dofs are zero. If
  that is the case, we obtain $$\begin{equation}
          \mathbf{{A}}_{tot} = \underline{\underline{\mathbf{M_{eq,w}}}} \; \underline{\underline{\mathbf{F_w}}} \; \mathbf{{I}}_{w,e} +  \underline{\underline{\mathbf{L_{eq}}}} \; \underline{\underline{\mathbf{F_{eq}}}} \; \mathbf{{I}}_{eq,e}
  \end{equation}$$ $$\begin{equation}
            \underline{\underline{\mathbf{K^T_{eq}}}} \, \mathbf{{A}}_{tot} = \underline{\underline{\mathbf{M_{eqi,we}}}} \; \mathbf{{I}}_{w,e} +  \underline{\underline{\mathbf{L_{eq,ie}}}} \; \mathbf{{I}}_{eq,e}
  \end{equation}$$

- *Fields to compare:*\
  We compare the Biot-Savart field with the new expression for the
  generalized halo coupling: $$\begin{equation}
   \label{eq:BiotSavartreffield}
          \mathbf{{B}}_{\theta} \Big\rvert_{ref} = \underline{\underline{\mathbf{Q_w}}} \, \mathbf{{I}}_w + \underline{\underline{\mathbf{Q_{eq}}}} \,  \mathbf{{I}}_{eq}
  \end{equation}$$ $$\begin{equation}
   \label{eq:halo_elec_field}
          \mathbf{{B}}_{\theta} \Big\rvert_{coupl} = \underline{\underline{\mathbf{Q^*_{we}}}} \, \mathbf{{I}}_{w,e} + \underline{\underline{\mathbf{Q^*_{eqe}}}} \,  \mathbf{{I}}_{eq,e} + \underline{\underline{\mathbf{Q_{vac}}}} \, \mathbf{{A}}_{tot}
  \end{equation}$$ Exploding these terms and equating
  [\[eq:BiotSavartreffield\]](#eq:BiotSavartreffield){reference-type="eqref"
  reference="eq:BiotSavartreffield"} to
  [\[eq:halo_elec_field\]](#eq:halo_elec_field){reference-type="eqref"
  reference="eq:halo_elec_field"}, we get: $$\begin{align}
   \label{eq:fieldsequate}
      & \Bigg[ \underline{\underline{\mathbf{Q_{\theta,w}}}} \, \mathbf{{I}}_w + \underline{\underline{\mathbf{Q_{\theta,eq}}}} \,  \mathbf{{I}}_{eq}\Big]_{ref}  =  \\ \notag
      & \left[ \left( \underline{\underline{\mathbf{Q_{\theta,w} }}} \, - \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,w}}}}   \right) \underline{\underline{\mathbf{F_w}}} \right] \, \underline{\mathbf{i_{w,e}}} - \\ \notag
      & - \left[ \,  \underline{\underline{\mathbf{Q_{\theta,eq} }}} \,  \left(  \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} - \underline{\underline{\mathbf{F_{eq}}}}     \right)  \left( \underline{\underline{\mathbf{F_{eq}^T}}} \underline{\underline{\mathbf{F_{eq}}}} \right)^{-1}   \right] \underline{\mathbf{i_{el,eq}}} + \\ \notag
      & + \left( \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \, \underline{\underline{\mathbf{K^T_{eq}}}} \, \underline{\underline{\mathbf{H}}} \right) \underline{\mathbf{A}}_{tot}
  \end{align}$$ In the framework of the test: $$\begin{equation}
      \underline{\underline{\mathbf{F_w}}} \; \underline{\mathbf{I_{we}}} = \underline{\mathbf{I_w}}
  \end{equation}$$ $$\begin{equation}
      \underline{\underline{\mathbf{F_{eq}}}} \; \underline{\mathbf{I_{eq,e}}} = \underline{\mathbf{I_{eq}}}
  \end{equation}$$ $$\begin{equation}
   \label{eq:new expr. for A}
            \underline{\underline{\mathbf{K^T_{eq}}}} \, \mathbf{{A}}_{tot} = \underline{\underline{\mathbf{M_{eqi,we}}}} \; \mathbf{{I}}_{w,e} +  \underline{\underline{\mathbf{L_{eq,ie}}}} \; \mathbf{{I}}_{eq,e}
  \end{equation}$$

  And since we are not projecting to the Bezier space, the
  $\underline{\underline{\mathbf{H}}}$ matrix in
  [\[eq:fieldsequate\]](#eq:fieldsequate){reference-type="eqref"
  reference="eq:fieldsequate"} is an identity and we can use expression
  [\[eq:new expr. for A\]](#eq:new expr. for A){reference-type="eqref"
  reference="eq:new expr. for A"}.\
  Hence the identity reads:

  $$\begin{align}
      &\cancel{\underline{\underline{\mathbf{Q_{\theta,w}}}} \, \mathbf{{I}}_w} + \cancel{\underline{\underline{\mathbf{Q_{\theta,eq}}}} \,  \mathbf{{I}}_{eq}} = \\ \notag
      &\cancel{\underline{\underline{\mathbf{Q_{\theta,w}}}} \, \mathbf{{I}}_w} - \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{{K}_{eq}}}} \, \underline{\underline{\mathbf{{L}_{eq,ii}^{-1}}}} \, \underline{\underline{\mathbf{{M}_{eqi,we}}}} \, \mathbf{{I}}_{w,e} - \\ \notag
      & - \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L_{eq,ii}^{-1}}}} \underline{\underline{\mathbf{L_{eq,ie}}}} \, \mathbf{{I}}_{eq,e} + \cancel{\underline{\underline{\mathbf{Q_{\theta,eq}}}} \, \, \mathbf{{I}}_{eq}} + \\ \notag 
      & + \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \; \underline{\underline{\mathbf{M_{eqi,we}}}} \; \mathbf{{I}}_{w,e} +  \underline{\underline{\mathbf{Q_{\theta,eq} }}} \, \underline{\underline{\mathbf{K_{eq}}}} \, \underline{\underline{\mathbf{L^{-1}_{eq,ii}}}} \; \underline{\underline{\mathbf{L_{eq,ie}}}} \; \mathbf{{I}}_{eq,e}
  \end{align}$$ **If the halo matrices are correctly computed, the two
  fields should match exactly.**

  ### [Test results]{.underline} {#test-results .unnumbered}

  We assign a current in the electrodes that follows a sinusoidal trend
  along $\phi$ and a normal distribution along the $\theta$ angle,
  centered on a single electrode. This would mimic the current crossing
  the electrodes stemming from a kinked plasma. The test is anyway
  independent of the particular current distribution.

  <figure id="fig:test_elec" data-latex-placement="h">
  <img src="./Figures/elec_test_curr.jpg" style="width:70.0%" />
  <figcaption>Input electrode currents for testing.</figcaption>
  </figure>

  <figure id="fig:test_fields" data-latex-placement="h">
  <img src="./Figures/test_elec_field.jpg" style="width:60.0%" />
  <figcaption>Reference and halo fields.</figcaption>
  </figure>

  The expressions match very well, with a relative error of $10^{-14}$.

# Electric potential $\Phi_w$ computation

Due to the inductive nature of the circuit, the link between currents
and potentials is not static. Therefore, the potentials should be
computed at each time step by inverting another system, projected along
the crossing dofs and not the solenoidal ones this time.\
The system will be $$\begin{align}
\mathbf{{L}}^{*}_{w,ee} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,e} + \mathbf{{R}}_{w,ee} \, \mathbf{{i}}_{w,e} + \color{blue} \mathbf{{F}}_w^T \mathbf{{F}}_w \Phi_w \color{black} = - \mathbf{{R}}_{w,ei} \, \mathbf{{i}}_{w,i} - \mathbf{{L}}^{*}_{w,ei} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{C}}_e^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot}   + \mathbf{{D}}_e^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,eq} + \mathbf{{F}}^T_{w} \, \mathbf{{v}} \label{eq:evolution_current_cross_01} \\
\mathbf{{F}}_w^T  \mathbf{{F}}_w \, \mathbf{{i}}_{w,e} = \mathbf{{i}}_{el,w} \label{eq:evolution_current_cross_02}
\end{align}$$

# Register

Here we note all the gradual steps followed for the coupling
implementation and testing.

## First coupling test

For the first test we plan to solve a reduced form of the wall current
equation

$$\begin{equation}
    \mathbf{{L}}^{*}_{w,ii} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = \mathbf{{C}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot} + \mathbf{{D}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,eq} + \mathbf{{E}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{el,w} + \mathbf{{G}}_i^* \, \mathbf{{i}}_{el,w} + \mathbf{{K}}_{w} \, \mathbf{{v}}
\end{equation}$$

$$\begin{equation}
    \mathbf{{i}}_{w,e} = \left( \mathbf{{F}}_w^T  \mathbf{{F}}_w \right)^{-1} \, \mathbf{{i}}_{el,w}
\end{equation}$$ Said reduced form is obtained by neglecting the mutual
interaction between subspaces. Note that we don't expect physically
relevant results from this, we only aim to prove that we can relax to
$I_p$-asymmetry constraint by means of injected current.\
This proof-of-concept system is the following: $$\begin{align}
 \label{eq:reduced halo 1}
    \mathbf{{L}}^{*}_{w,ii} \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{i}}_{w,i} + \mathbf{{R}}_{w,ii} \, \mathbf{{i}}_{w,i} = \mathbf{{C}}_i^* \, %
\frac{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}^{}}{\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}t^{}} \mathbf{{a}}_{tot} \\
     \mathbf{{i}}_{w,e} = \left( \mathbf{{F}}_w^T  \mathbf{{F}}_w \right)^{-1} \, \mathbf{{i}}_{el,w}
\end{align}$$ And once we solve separately for $\mathbf{{i}}_{w,i}$ and
$\mathbf{{i}}_{w,e}$, we obtain the extended current $$\begin{equation}
    \mathbf{i}_w = \underline{\underline{\mathbf{K^T_w}}} \, \mathbf{i}_{w.i} + \underline{\underline{\mathbf{F^T_w}}} \, \mathbf{i}_{w.e}
\end{equation}$$

Note that this must be implemented in the transformed space where
$\mathbf{{L}}^{*}_{w,ii}$ is transformed to diagonal form
$\underline{\underline{\mathbf{D_{ii}}}}$ So what we have is:
$$\begin{equation}
    \, \underline{\mathbf{\dot{Y_i}}} = -  \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\mathbf{{Y_i}}} + \mu_0 \, \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\underline{\mathbf{S^T_{ii}}}} \, \underline{\underline{\mathbf{C^*_i}}} \,\dot{\underline{\mathbf{a_{tot}}}}
\end{equation}$$ $$\begin{equation}
    \, \underline{\mathbf{\dot{Y_i}}} = -  \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\mathbf{{Y_i}}} + \mu_0 \, \underline{\underline{\mathbf{M_{yi,e}}}} \, {\underline{\mathbf{\dot{a}_{tot}}}}
\end{equation}$$ When introducing also the new
$\underline{\underline{\mathbf{Gi^*}}}$, we obtain $$\begin{equation}
    \, \underline{\mathbf{\dot{Y_i}}} = -  \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\mathbf{{Y_i}}} + \mu_0 \, \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\underline{\mathbf{S^T_{ii}}}} \, \underline{\underline{\mathbf{C^*_i}}} \,\dot{\underline{\mathbf{a_{tot}}}} + \mu_0 \, \underline{\underline{\mathbf{D^{-1}_{ii}}}} \, \underline{\underline{\mathbf{S^T_{ii}}}} \, \underline{\underline{\mathbf{G^*_i}}} \, \underline{\mathbf{i_{el,w}}}
\end{equation}$$

## Differences between matching quantities

**Problem: Mismatch between `I3D_ideal_wall_response` test for `M_ye`**

<figure id="fig:I3D_resp" data-latex-placement="h">
<img src="./Figures/I3D_resp.png" style="width:70.0%" />
<figcaption>Response wall current from a kick in vector
potential</figcaption>
</figure>

**N.B.** `M_ye` is not involved in free-boundary equilibrium computation
(now working).

<figure id="fig:flux_comp" data-latex-placement="h">
<img src="./Figures/flux_comp.PNG" style="width:70.0%" />
<figcaption>Free-boundary flux surfaces comparison</figcaption>
</figure>

- The `I_3D_ideal_wall_response` test is calculated as

  Eddy: $$\begin{equation}
      \mathbf{{I}}_{resp} = \mathbf{{S}} \cdot \mathbf{{S}}^T \, \mathbf{{M}}_{w,eq} \, \mathbf{{L}}^{-1}_{eq} \, \mathbf{{H}} \,  / \, \mathbf{{D}}
  \end{equation}$$

  Halo: $$\begin{equation}
      \mathbf{{I}}_{resp,i} = \mathbf{{S}}_{ii} \cdot \mathbf{{S}}_{ii}^T \, \mathbf{{K}}^T_w \, \mathbf{{M}}_{w,eq} \, \mathbf{{K}}_{eq} \, \mathbf{{L}}^{-1}_{eq,ii} \, \mathbf{{K}}^T_{eq} \, \mathbf{{H}}
  \end{equation}$$

  **But like this, the response for the halo is just noise!**

- $\mathbf{{M}}_{e,yi} \, \mathbf{{Y}}_i$ and
  $\mathbf{{M}}_{e,y} \, \mathbf{{Y}}$ yield the same equilibrium
  tangent field, but the individual factors are different!\
  $$\begin{equation}
          \mathbf{{M}}_{e,yi} \neq \mathbf{{M}}_{e,y}
  \end{equation}$$ $$\begin{equation}
          \mathbf{{Y}}_i \neq \mathbf{{Y}}
  \end{equation}$$ but $$\begin{equation}
          \mathbf{{B}}_{tan} = \mathbf{{M}}_{e,yi} \mathbf{{Y}}_i \rvert_{eddy+halo}= \mathbf{{M}}_{e,y} \mathbf{{Y}} \rvert_{eddy}
  \end{equation}$$ which implies $$\begin{equation}
          \left\{ \left[ \mathbf{{S}}^{T} \, \mathbf{{M}}_{w,eq} \, \mathbf{{L}}^{-1}_{eq} \, \mathbf{{H}} \right] / \mathbf{{D}} \right\} \, \cdot \mathbf{{S}}^{-1} \, \mathbf{{I}}_{eddy}  =    \left\{ \left[ \mathbf{{S}}_{ii}^{T} \, \mathbf{{K}}^T_w \, \mathbf{{M}}_{w,eq} \, \mathbf{{K}}_{eq} \, \mathbf{{L}}^{-1}_{eq,ii} \, \mathbf{{H}}_i \right] / \mathbf{{D}}_{ii} \right\} \cdot \mathbf{{S}}_{ii}^{-1} \, \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right) \, \mathbf{{K}}^T_w \; \mathbf{{I}}_{halo}
  \end{equation}$$ And the relation between the transformed currents
  $\mathbf{{Y}}$ and $\mathbf{{Y}}_i$ is: $$\begin{equation}
          \mathbf{{Y}}_i = \textit{pinv}(\mathbf{{M}}_{e,yi}) \,  \mathbf{{M}}_{e,y} \; \mathbf{{Y}}
  \end{equation}$$ Let observe that $$\begin{equation}
          \mathbf{{Y}}_i = \mathbf{{S}}_{ii}^{-1} \, \left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w \right] \, \mathbf{{I}}_{{eddy+halo}}
  \end{equation}$$ $$\begin{equation}
          \mathbf{{Y}} = \mathbf{{S}}^{-1} \, \mathbf{{I}}_{{eddy}}
  \end{equation}$$

- **Matrix `M_ye`**

  In the eddy coupling, we use the following matrix:

  $$\begin{equation}
      \mathbf{{M}}_{y,e} = 
      \left[ \mathbf{{S}}^{T} \, \mathbf{{M}}_{w,eq} \, \mathbf{{L}}^{-1}_{eq} \, \mathbf{{H}} \right] / \mathbf{{D}}
  \end{equation}$$ where $\mathbf{{D}}$ is the diagonal matrix with the
  eigenvalues of the system. The **eigenvalues match between the eddy
  coupling scheme and the internal dofs of the halo coupling scheme**.

  In the halo coupling, we use the following matrix: $$\begin{equation}
      \mathbf{{M}}_{yi,e} = 
      \left[ \mathbf{{S}}_{ii}^{T} \cdot \mathbf{{K}}^T_w \, \mathbf{{M}}_{w,eq} \, \mathbf{{K}}_{eq} \, \mathbf{{L}}^{-1}_{eq,ii} \, \mathbf{{H}}_i \right] / \mathbf{{D}}_{ii}
  \end{equation}$$

  Is this pre-multiplication by $\mathbf{{S}}^T_{ii}$ correct? We have
  that $$\begin{align}
      \mathbf{{Y}} &= \mathbf{{S}}^{-1} \, \mathbf{{I}}_{eddy} \\
      \mathbf{{Y}}_i &= \mathbf{{S}}_{ii}^{-1} \, \mathbf{{I}}_{i}
  \end{align}$$ with
  $\mathbf{{I}}_i = \left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w  \right] \, \mathbf{{I}}_{eddy+halo}$

The complete transformation reads as $$\begin{equation}
    \mathbf{Y}_i = \underbrace{\mathbf{S}^{-1}_{ii}  \left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w  \right] } \, \mathbf{{I}}
\end{equation}$$ And to go back: $$\begin{equation}
    \mathbf{{S}}_{ii} \, \mathbf{{Y}}_i = \left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w  \right] \, \mathbf{{I}}
\end{equation}$$ Since
$\mathbf{{K}}^{+}=\left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w  \right]$
is the left pseudoinverse, then $$\begin{equation}
    \mathbf{{K}}^{+} \, \mathbf{{K}} = identity
\end{equation}$$ Let us express $$\begin{equation}
    \mathbf{{I}} = \mathbf{{K}} \, \mathbf{{I}}_i + \mathbf{{F}} \, \mathbf{{I}}_e
\end{equation}$$ then this yields $$\begin{equation}
    \mathbf{{S}}_{ii} \, \mathbf{{Y}}_i = \left[ \left( \mathbf{{K}}^T_w \, \mathbf{{K}}_w  \right)^{-1} \, \mathbf{{K}}^T_w  \right] \, \mathbf{{K}} \, \mathbf{{I}}_i
\end{equation}$$ hence $$\begin{equation}
    \mathbf{{S}}_{ii} \, \mathbf{{Y}}_i = \mathbf{{I}}_i
\end{equation}$$ therefore $$\begin{equation}
    \mathbf{{I}} = \mathbf{{K}} \, \mathbf{{S}}_{ii} \, \mathbf{{Y}}_i
\end{equation}$$

# JOREK boundary conditions for the velocity

We investigate first boundary conditions for the velocity in the full
MHD framework. In first paragraphs of this Section we recall the
representation of JOREK variables and the weak formulation of the
conservation of momentum equation, taking results from Reference
[@Hoelzl2021]. The vector field for the velocity is represented by:

$$\begin{equation}
    \label{eq:representation_velocity}
    \mathbf{v} = v_R \mathbf{e}_R + v_Z \mathbf{e}_Z + v_{||} \mathbf{B}
\end{equation}$$ The Weyl gauge is adopted, so the magnetic vector
potential has three degrees of freedom and is represented by:
$$\begin{equation}
    \label{eq:representation_A}
    \mathbf{A} = A_R \mathbf{e}_R + A_Z \mathbf{e}_Z + \frac{\psi}{R} \mathbf{e}_{\varphi}
\end{equation}$$

Each of the scalar potentials $(v_R,v_Z,v_{||})$ and $(A_R,A_Z,\psi)$
introduced in equations
([\[eq:representation_velocity\]](#eq:representation_velocity){reference-type="ref"
reference="eq:representation_velocity"})-([\[eq:representation_A\]](#eq:representation_A){reference-type="ref"
reference="eq:representation_A"}) will be approximated via the standard
Bézier-Fourier basis functions typical for scalar fields in the JOREK
finite element space. The strong form for the force balance reads:
$$\begin{equation}
    \label{eq:strong_form}
    \rho %
\frac{\partial^{}\mathbf{v}}%
{\partial t^{}} = -\rho \mathbf{v}\cdot \nabla \mathbf{v} - \nabla p + \mathbf{j} \times \mathbf{B} + \nabla \cdot \mathbf{\tau} + \mathbf{S}_{v}
\end{equation}$$ Clearly the weak form resulting from projecting
([\[eq:strong_form\]](#eq:strong_form){reference-type="ref"
reference="eq:strong_form"}) on a test field which is given by a
standard scalar test function and any of the directions
$\mathbf{e}_X \in \{ \mathbf{e}_R,\mathbf{e}_z,\mathbf{B} \}$ is valid:
$$\begin{equation}
    \int_{V}{v^* \mathbf{e}_X \cdot \rho %
\frac{\partial^{}\mathbf{v}}%
{\partial t^{}} }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{V} = \int_{V}{v^* \mathbf{e}_X \cdot [-\rho \mathbf{v}\cdot \nabla \mathbf{v} - \nabla p + \mathbf{j} \times \mathbf{B} + \nabla \cdot \mathbf{\tau} + \mathbf{S}_{v}] }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{V}
\end{equation}$$ In particular, for what concerns the
$\mathbf{j}\times \mathbf{B}$ product, using Ampére's law results in
[@Braams1991]: $$\begin{align}
    \int_{V}{v^* \mathbf{e}_X \cdot \mathbf{j} \times \mathbf{B} }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{V} = & \int_{V}{ \frac{B^2}{2 \mu_0} \nabla \cdot \left[ v^* \mathbf{e}_X \right] - \mathbf{B} \cdot \left( \mathbf{B} \cdot \nabla \right) \left( v^* \mathbf{e}_X \right) }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{V} \nonumber \\
    - & \int_{+ \partial V}{ \frac{B^2}{2 \mu_0}  \left[ v^* \mathbf{e}_X \cdot \hat{\mathbf{n}} \right]  - \left( \mathbf{B} \cdot v^* \mathbf{e}_X \right)\left( \mathbf{B} \cdot \hat{\mathbf{n}} \right) }\@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{S}
\end{align}$$ Hence a boundary term appears, besides the usual
volumetric term. Using Ampére's law is recommended to not deal with
$\mathbf{j}$ directly (hence with second order derivatives and
pontential discontinuous functions.

# Divergence of J

$$\begin{equation}
    \mathbf{j} = J_R(R,\varphi,z) \nabla R + j_Z(R,\varphi,z) \nabla z + R j_\varphi(R,\varphi,z) \nabla \varphi
\end{equation}$$

Rappresentazione di Fourier lungo l'angolo toroidale $$\begin{equation}
    \mathbf{j} = \sum_{k=0}^{\infty} j_{R,k}(R,z) \cos{\left(k \varphi + \phi_{r,k}(R,z)\right)} \nabla R + j_{z,k}(R,z) \cos{\left(k \varphi + \phi_{z,k}(R,z)\right)} \nabla z + R j_{\varphi,k}(R,z) \cos{\left(k \varphi + \phi_{\varphi,k}(R,z)\right)} \nabla \varphi
\end{equation}$$

In notazione esponenziale: $$\begin{equation}
    \mathbf{j} = \sum_{k=-\infty}^{\infty} \left[ \mathbf{J}_{R,k}(R,z) \nabla R + \mathbf{J}_{z,k}(R,z) \nabla z + R \mathbf{J}_{\varphi,k}(R,z) \nabla \varphi\right] \exp{\left(-i k \varphi\right)}
\end{equation}$$

Divergenza di J: $$\begin{equation}
    \nabla \cdot \mathbf{j}_k = \left[ \frac{1}{R} %
\frac{\partial^{}}%
{\partial R^{}} \left(R \mathbf{J}_{R,k}\right)  + %
\frac{\partial^{}}%
{\partial Z^{}} \left( \mathbf{J}_{Z,k}\right) - \frac{i k}{R} \mathbf{J}_{\varphi,k}(R,Z)   \right] \exp{\left( - {i k} \varphi \right)}
\end{equation}$$

The imposition of a divergence-free current density hence reads:
$$\begin{equation}
    \nabla \cdot \left[  \mathbf{J}_{R,k}(R,z) \nabla R + \mathbf{J}_{z,k}(R,z) \nabla z \right] - \frac{i k}{R} \mathbf{J}_{\varphi,k}(R,z) = 0
\end{equation}$$ Let us integrate over the volume: $$\begin{equation}
    \int_{0}^{2 \pi} \int_{\Gamma} \mathbf{J}_{n,k}(R,z) R d \ell d \varphi - {(i k)} \int_{0}^{2\pi} \int_{\Omega}{ \mathbf{J}_{\varphi,k}(R,z) }  \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{S} \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{\varphi} = 0
\end{equation}$$ Hence: $$\begin{equation}
    \int_{\Gamma}{\mathbf{J}_{n,k}(R,z) R \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{\ell}} - {(i k )} \int_{\Omega}{ \mathbf{J}_{\varphi,k}(R,z) }  \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{S} = 0
\end{equation}$$ and this implies that, for harmonics coming out of beat
frequencies, for which the $J_{\phi}$ current is not simulated, we
should simply have $$\begin{equation*}
    \int_{\Gamma}{\mathbf{J}_{n,k}(R,z) R \@ifnextchar^{\DIfF}{%
\mathop{\mathrm{\mathstrut d}}%
\nolimits^{}%
\futurelet\diffarg%
\let\DiffSpace\!%
\ifx\diffarg(%
\let\DiffSpace\relax
\else
\ifx\diffarg[%
\let\DiffSpace\relax
\else
\ifx\diffarg\{%
\let\DiffSpace\relax
\fi\fi\fi\DiffSpace}{\ell}} = 0
\end{equation*}$$

# Definitions {#definitions .unnumbered}

In equation
([\[eq:equivalent_current_step1\]](#eq:equivalent_current_step1){reference-type="ref"
reference="eq:equivalent_current_step1"}) we introduced:

$$\begin{align}
    \mathbf{{L}}_{eq,ii} & = \mathbf{{K}}_{eq}^T \, \mathbf{{L}}_{eq} \, \mathbf{{K}}_{eq} \label{eq:definition_1a}\\
    \mathbf{{L}}_{eq,ie} & = \mathbf{{K_{eq}}}^T \, \mathbf{{L}}_{eq} \, \mathbf{{F}}_{eq} \\
    \mathbf{{L}}_{eq,ei} & = \mathbf{{F}}_{eq}^T \, \mathbf{{L}}_{eq} \, \mathbf{{K}}_{eq} \\
    \mathbf{{L}}_{eq,ee} & = \mathbf{{F}}_{eq}^T \, \mathbf{{L}}_{eq} \, \mathbf{{F}}_{eq} \\
    \mathbf{{H}}_{i} & = \mathbf{{K}}_{eq}^T \, \mathbf{{H}} \\
    \mathbf{{H}}_{e} & = \mathbf{{F}}_{eq}^T \, \mathbf{{H}} \\
    \mathbf{{M}}_{eqi,w} & = \mathbf{{K}}_{eq}^T \, \mathbf{{M}}_{eq,w} \\
    \mathbf{{M}}_{eqe,w} & = \mathbf{{F}}_{eq}^T \, \mathbf{{M}}_{eq,w} \label{eq:definition_1h} \\
\end{align}$$

In equation
([\[eq:evolution_current_null_01\]](#eq:evolution_current_null_01){reference-type="ref"
reference="eq:evolution_current_null_01"}) we have introduced the
matrices:

$$\begin{align}
    \mathbf{{L}}^{*}_{w,ii} = \mathbf{{K}}_w^T \, \mathbf{{L}}^{*}_{w} \, \mathbf{{K}}_w = \mathbf{{K}}_w^T \, \left( \mathbf{{L}}_w - \mathbf{{M}}_{w,eqi} \, \mathbf{{L}}_{eq,ii} \, \mathbf{{M}}_{eqi,w}  \right) \, \mathbf{{K}}_w    \\
    \mathbf{{L}}^{*}_{w,ie} = \mathbf{{K}}_w^T \, \mathbf{{L}}^{*}_{w} \, \mathbf{{F}}_w \\
    %\Amat[]{M_{wi,eqe}} = \Amat[]{M_{eqe,wi}}^T \\
    \mathbf{{R}}_{w,ii} = \mathbf{{K}}_w^T \, \mathbf{{R}}_{w} \, \mathbf{{K}}_w \\
    \mathbf{{R}}_{w,ie} = \mathbf{{K}}_w^T \, \mathbf{{R}}_{w} \, \mathbf{{F}}_w
\end{align}$$

# Various cuts to delete {#various-cuts-to-delete .unnumbered}

$$\begin{equation}
\label{eq:equivalent_current_per_components}
    \left[ \begin{array}{cccc}
        \mathbf{{L}}_{eq,ii} & \mathbf{{L}}_{eq,ie} & 0 & \mathbf{{M}}_{eqi,we} \\
        \mathbf{{L}}_{eq,ei} & \mathbf{{L}}_{eq,ee} & \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & \mathbf{{M}}_{eqe,we} \\
        0 &   \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & 0 & -\mathbf{{F}}_{w}^T \, \mathbf{{F}}_{w} \\
    \end{array} \right] \, \left[ 
    \begin{array}{c}
         \mathbf{{i}}_{eq,i}  \\
         \mathbf{{i}}_{eq,e}  \\
         \mathbf{{\chi}}_{eq} \\
         \mathbf{{i}}_{w,e}
    \end{array}
    \right]  = 
    \left[
    \begin{array}{c}
        \mathbf{{H_i}}   \\
        \mathbf{{H_e}}   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{a}}_{tot} -
    \left[
    \begin{array}{c}
        \mathbf{{M}}_{eqi,wi}   \\
        \mathbf{{M}}_{eqe,wi}   \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{i}}_{wi}
\end{equation}$$ The matrices here are defined again in Appendix, see
equations
([\[eq:definition_2a\]](#eq:definition_2a){reference-type="ref"
reference="eq:definition_2a"})-([\[eq:definition_2c\]](#eq:definition_2c){reference-type="ref"
reference="eq:definition_2c"}). Clearly the above system has much more
unknowns than constraints, as we did not explicitly impose the value of
injected electrode currents or of time integrated potentials.
$$\begin{equation}
\label{eq:equivalent_current_per_components_to_delete}
    \left[ \begin{array}{cccc}
        \mathbf{{L}}_{eq,ii} & \mathbf{{L}}_{eq,ie} & 0 & \mathbf{{M}}_{eqi,we} \\
        \mathbf{{L}}_{eq,ei} & \mathbf{{L}}_{eq,ee} & \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & \mathbf{{M}}_{eqe,we} \\
        0 &   \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & 0 & -\mathbf{{F}}_{w}^T \, \mathbf{{F}}_{w} \\
        0 &   \mathbf{{F}}_{eq}^T \, \mathbf{{F}}_{eq} & 0 & 0 \\
    \end{array} \right] \, \left[ 
    \begin{array}{c}
         \mathbf{{i}}_{eq,i}  \\
         \mathbf{{i}}_{eq,e}  \\
         \mathbf{{\chi}}_{eq} \\
         \mathbf{{i}}_{w,e}
    \end{array}
    \right]  = 
    \left[
    \begin{array}{c}
        \mathbf{{H_i}}   \\
        \mathbf{{H_e}}   \\
        \mathbf{{0}}     \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{a}}_{tot} -
    \left[
    \begin{array}{c}
        \mathbf{{M}}_{eqi,wi}   \\
        \mathbf{{M}}_{eqe,wi}   \\
        \mathbf{{0}}            \\
        \mathbf{{0}}
    \end{array}
    \right] \, \mathbf{{i}}_{wi}   +
    \left[ \begin{array}{c}
        \mathbf{{0}}  \\
        \mathbf{{0}}  \\
        \mathbf{{0}}  \\
        \mathbf{{i}}_{el,eq}
    \end{array}
    \right]
\end{equation}$$ In equation
([\[eq:equivalent_current_per_components\]](#eq:equivalent_current_per_components){reference-type="ref"
reference="eq:equivalent_current_per_components"}) we introduced:

$$\begin{align}
    \mathbf{{M_{eqi,wi}}} = \mathbf{{M_{eqi,w}}} \, \mathbf{{K_w}} \label{eq:definition_2a}\\
    \mathbf{{M_{eqi,we}}} = \mathbf{{M_{eqi,w}}} \, \mathbf{{F_w}} \\
    \mathbf{{M_{eqe,wi}}} = \mathbf{{M_{eqe,w}}} \, \mathbf{{K_w}} \\
    \mathbf{{M_{eqe,we}}} = \mathbf{{M_{eqe,w}}} \, \mathbf{{F_w}} \label{eq:definition_2c}
\end{align}$$
