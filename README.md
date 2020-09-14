![PyPI - License](https://img.shields.io/pypi/l/JumpDiff)
![PyPI](https://img.shields.io/pypi/v/JumpDiff)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/JumpDiff)
[![Build Status](https://travis-ci.org/LRydin/JumpDiff.svg?branch=master)](https://travis-ci.org/LRydin/JumpDiff)
[![codecov](https://codecov.io/gh/LRydin/JumpDiff/branch/master/graph/badge.svg)](https://codecov.io/gh/LRydin/JumpDiff)
[![Documentation Status](https://readthedocs.org/projects/jumpdiff/badge/?version=latest)](https://jumpdiff.readthedocs.io/en/latest/?badge=latest)

# JumpDiff
`JumpDiff` is a `python` library with non-parametric Nadaraya─Watson estimators to extract the parameters of jump-diffusion processes.
With `JumpDiff` one can extract the parameters of a jump-diffusion process from one-dimensional timeseries, employing both a kernel-density estimation method combined with a set on second-order corrections for a precise retrieval of the parameters for short timeseries.

## Installation
To install `JumpDiff`, run

```
   pip install JumpDiff
```

Then on your favourite editor just use

```python
   import JumpDiff as jd
```

## Dependencies
The library parameter estimation depends on `numpy` and `scipy` solely. The mathematical formulae depend on `sympy`. It stems from [`kramersmoyal`](https://github.com/LRydin/KramersMoyal) project, but functions independently from it<sup>3</sup>.

## Documentation
You can find the documentation [here](https://jumpdiff.readthedocs.io/).

# Jump-diffusion processes
## The theory
Jump-diffusion processes<sup>1</sup>, as the name suggest, are a mixed type of stochastic processes with a diffusive and a jump term.
One form of these processes which is mathematically traceable is given by the [Stochastic Differential Equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation)

<img src="/Others/SDE_1.png" title="A jump diffusion process" height="25"/>

which has 4 main elements: a drift term <img src="/Others/a_xt.png" title="drift term" height="18"/>, a diffusion term <img src="/Others/b_xt.png" title="diffusion term" height="18"/>, and jump amplitude term <img src="/Others/xi.png" title="jump amplitude term" height="18"/>, which is given by a Gaussian distribution, and finally a jump rate <img src="/Others/lambda.png" title="jump rate term" height="14"/>.
You can find a good review on this topic in Ref. 2.

## Integrating a jump-diffusion process
Let us use the functions in `JumpDiff` to generate a jump-difussion process, and subsequently retrieve the parameters. This is a good way to understand the usage of the integrator and the non-parametric retrieval of the parameters.

First we need to load our library. We will call it `jd`
```python
import JumpDiff as jd
```
Let us thus define a jump-diffusion process and use `jdprocess` to integrate it. Do notice here that we need the drift <img src="/Others/a_xt.png" title="drift term" height="18"/> and diffusion <img src="/Others/b_xt.png" title="diffusion term" height="18"/> as functions.

```python
# integration time and time sampling
t_final = 10000
delta_t = 0.001

# A drift function
def a(x):
    return -0.5*x

# and a (constant) diffusion term
def b(x):
    return 0.75

# Now define a jump amplitude and rate
xi = 2.5
lamb = 1.75

# and simply call the integration function
X = jd.jdprocess(t_final, delta_t, a=a, b=b, xi=xi, lamb=lamb)
```

This will generate a jump diffusion process `X` of length `int(10000/0.001)` with the given parameters.

<img src="/Others/X_trajectory.png" title="A jump-difussion process" height="200"/>

## Using `JumpDiff` to retrieve the parameters
### Moments and Kramers─Moyal coefficients
Take the timeseries `X` and use the function `moments` to retrieve the conditional moments of the process.
For now let us focus on the shortest time lag, so we can best approximate the Kramers─Moyal coefficients.
For this case we can simply employ

```python
edges, moments = jd.moments(timeseries = X)
```
In the array `edges` are the limits of our space, and in our array `moments` are recorded all 6 powers/order of our conditional moments.
Let us take a look at these before we proceed, to get acquainted with them.

We can plot the first moment with any conventional plotter, so lets use here `plotly` from `matplotlib`

```python
import matplotlib.plotly as plt

# we want the first power, so we need 'moments[1,...]'
plt.plot(edges, moments[1,...])
```
The first moment here (i.e., the first Kramers─Moyal coefficient) is given solely by the drift term that we have selected `-0.5*x`

<img src="/Others/1_moment.png" title="The 1st Kramers─Moyal coefficient" height="200"/>

And the second moment (i.e., the second Kramers─Moyal coefficient) is a mixture of both the contributions of the diffusive term <img src="/Others/b_xt.png" title="diffusion term" height="18"/> and the jump terms <img src="/Others/xi.png" title="jump amplitude term" height="18"/> and <img src="/Others/lambda.png" title="jump rate term" height="14"/>.

<img src="/Others/2_moment.png" title="The 2nd Kramers─Moyal coefficient" height="200"/>

You have this stored in `moments[2,...]`.

### Retrieving the jump-related terms
Naturally one of the most pertinent questions when addressing jump-diffusion processes is the possibility of recovering these same parameters from data. For the given jump-diffusion process we can use the `jump_amplitude` and `jump_rate` functions to non-parametrically estimate the jump amplitude <img src="/Others/xi.png" title="jump amplitude term" height="18"/> and jump rate <img src="/Others/lambda.png" title="jump rate term" height="18"/> terms.

After having the `moments` in hand, all we need is

```python
# first estimate the jump amplitude
xi_est = jd.jump_amplitude(moments = moments)

# and now estimated the jump rate
lamb_est = jd.jump_rate(moments = moments)
```
which resulted in our case in `(xi_est) ξ = 2.43 ± 0.17` and `(lamb_est) λ = 1.744 * delta_t` (don't forget to divide `lamb_est` by `delta_t`)!

### Other functions and options
Include in this package is also the [Milstein scheme](https://en.wikipedia.org/wiki/Milstein_method) of integration, particularly important when the diffusion term has some spacial `x` dependence. `moments` can actually calculate the conditional moments for different lags, using the parameter `lag`.

In `formulae` the set of formulas needed to calculate the second order corrections are given (in `sympy`).

# Contributions
We welcome reviews and ideas from everyone. If you want to share your ideas, upgrades, doubts, or simply report a bug, open an [issue](https://github.com/LRydin/JumpDiff/issues) here on GitHub, or contact us directly.
If you need help with the code, the theory, or the implementation, drop us an email.
We abide to a [Conduct of Fairness](contributions.md).

# Changelog
- Version 1.0 - releasing the first full version.
- Version 0.4 - Designing a set of self-consistency checks, the documentation, examples, and a trial code.
- Version 0.3 - Designing a straightforward procedure to retrieve the jump amplitude and jump rate functions, alongside with a easy `sympy` displaying the correction.
- Version 0.2 - Introducing the second-order corrections to the moments
- Version 0.1 - Design an implementation of the `moments` functions, generalising `kramersmoyal` `km`.

# Literature and Support

### History
This project was started in 2017 at the [neurophysik](https://www.researchgate.net/lab/Klaus-Lehnertz-Lab-2) by Leonardo Rydin Gorjão, Jan Heysel, Klaus Lehnertz, and M. Reza Rahimi Tabar, and separately by Pedro G. Lind, at the Department of Computer Science, Oslo Metropolitan University. Pedro, Leonardo, and Dirk developed in 2019 and 2020 a set of corrections and an implementation for python, presented here.

### Funding
Helmholtz Association Initiative _Energy System 2050 - A Contribution of the Research Field Energy_ and the grant No. VH-NG-1025 and *STORM - Stochastics for Time-Space Risk Models* project of the Research Council of Norway (RCN) No. 274410.

---
##### Bibliography

<sup>1</sup> Tabar, M. R. R. *Analysis and Data-Based Reconstruction of Complex Nonlinear Dynamical Systems.* Springer, International Publishing (2019), Chapter [*Stochastic Processes with Jumps and Non-vanishing Higher-Order Kramers–Moyal Coefficients*](https://doi.org/10.1007/978-3-030-18472-8_11).

<sup>2</sup> Friedrich, R., Peinke, J., Sahimi, M., Tabar, M. R. R. *Approaching complexity by stochastic methods: From biological systems to turbulence,* [Physics Reports 506, 87–162 (2011)](https://doi.org/10.1016/j.physrep.2011.05.003).

<sup>3</sup> Rydin Gorjão, L., Meirinhos, F. *kramersmoyal: Kramers–Moyal coefficients for stochastic processes.* [Journal of Open Source Software, **4**(44) (2019)](https://doi.org/10.21105/joss.01693).

##### Extended Literature
You can find further reading on SDE, non-parametric estimatons, and the general principles of the Fokker–Planck equation, Kramers–Moyal expansion, and related topics in the classic (physics) books

- Risken, H. *The Fokker–Planck equation.* Springer, Berlin, Heidelberg (1989).
- Gardiner, C.W. *Handbook of Stochastic Methods.* Springer, Berlin (1985).

And an extensive review on the subject [here](http://sharif.edu/~rahimitabar/pdfs/80.pdf)
