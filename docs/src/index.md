# SOCP Trajectory Optimization Documentation

The key functionality is the primal SOCP solver for Trajectory Optimization.

```@meta
CurrentModule = TrajOptSOCPs
```

The SOCP trajectory optimization problem is formulated as follows.

```math
\begin{aligned}
\underset{x_{1:N}, u_{1:N-1}}{\text{minimize}} \quad& \frac{1}{2} (x_N - x_{N, ref})^\top Q (x_N - x_{N, ref}) + \\
& \frac{1}{2} \Sigma_{n = 1}^{N-1} (x_n - x_{n, ref})^\top Q (x_n - x_{n, ref})
+ (u_n - u_{n, ref})^\top R (u_n - u_{n, ref}))  \\
\text{subject to} \quad& x_{k + 1} = A x_k + B u_k + C \\
& H [x; u] ≤ h \\
& ||u_k|| ≤ u_{max}
\end{aligned}
```

Table of Contents:
```@contents
Pages   = ["UI/ui.md"]
```


```@docs
TrajOptSOCPs
```
