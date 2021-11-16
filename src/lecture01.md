---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: SageMath 9.3
  language: sage
  name: sagemath
---

# Lecture 1 : Basic ideas

The concept of differentiation and integration to noninteger order goes as far back as the concept we are used to work with. Leibniz mentiones it in a letter to L'Hospital in 1695:

> "John Bernoulli seems to have told you of my having mentioned to him a
> marvelous analogy which makes it possible to
>				say in a way that successive differentials are
>				in geometric progression. **One can ask what
>				would be a differential having as its exponent
>				a fraction**. You see that the result can be expressed
>				by an infinite series. Although this
>				seems removed from Geometry, which does
>				not yet know of such fractional exponents, **it
>				appears that one day these paradoxes will
>				yield useful consequences, since there is hardly
>				a paradox without utility**. Thoughts that
>				mattered little in themselves may give occasion
>				to more beautiful ones.''

On the other hand, we have to wait until the beginning of the 19th century to have a somewhat systematic study of the matters, specifically the works by Liouville and Riemann, even if earlier cotributions by Euler, Lagrange and others had already been made. For an historical survey one could refere to Chapter 1 in {cite}`MR0361633`.

Since we have a limited amount of time, we will build the instrument we need
in an ahistorical way, that is, preferring a utilitarian approach.

We start by introducing a **special function** that will be ubiquitous in our
computations: the $\Gamma$ function.

````{prf:definition}
:label: gammafunction

The **$\Gamma$ function** $\Gamma(z)$ is defined for complex numbers with a
positive real part via the convergent improper integral:
:::{math}
\Gamma(z) = \int_{0}^{+\infty} x^{z-1} e^{-x}\,{\rm d}x, \qquad \Re(z) > 0,
:::
and then extended by *analyitic continuation* to a meromorphic function that is holomorphic in the whole complex plane except zero and the negative integers, where the function has simple poles.

````

```{code-cell} sagemath
:tags: [gammafunction]

plot(gamma(x), (x,-5,5), ymin=-5,ymax=5,exclude = [-5..5])
```

This is the most "popular" and useful extension of the factorial function over the complex, indeed other do exist, but the $\Gamma$ is the one we need for our tasks.

The next ingredient we need is a formula for performing repeated integration of
a function.

````{prf:lemma}
:label: swappingintegrals
If $G(x,t)$ is jointly continuous on $[c,b]\times[c,b]$:
```{math}
\int_{c}^{x} {\rm d}x_1 \int_{c}^{x_1} G(x_1,t) {\rm d}t = \int_{c}^{x} {\rm d}t \int_{t}^{x}G(x_1,t) {\rm d}x_1.
```
````

That we can use to **prove by induction** a formula for the $n$th-fold integral
of a function.

````{prf:proposition}
Given a **continuous** function $f : [c,b] \rightarrow \mathbb{R}$ we have for $x < b$
```{math}
I^n_c(x) = \int_{c}^{x}{\rm d}x_1 \int_{c}^{x_1} {\rm d}x_2 \int_{c}^{x_2} {\rm d}x_3 \cdots \int_{c}^{x_{n-1}}f(t)dt = \frac{1}{(n-1)!} \int_{c}^{x} (x-t)^{n-1} f(t){\rm d}t
```
````

````{prf:proof}
We want to use {prf:ref}`swappingintegrals`
````


## Riemann-Liouville Derivative



### Finite differences approximation

## Caputo Derivative

### Finite differences approximation

## Riesz Derivative

### Finite differences approximation

+++

# Bibliography

```{bibliography}
:filter: docname in docnames
```
