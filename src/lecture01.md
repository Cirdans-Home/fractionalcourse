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

## Fractional Integrals

We start by introducing a **special function** that will be ubiquitous in our
computations: the $\Gamma$ function.

````{prf:definition}
:label: Gamma-Function

The **$\Gamma$ function** $\Gamma(z)$ is defined for complex numbers with a
positive real part via the convergent improper integral:
:::{math}
\Gamma(z) = \int_{0}^{+\infty} x^{z-1} e^{-x}\,{\rm d}x, \qquad \Re(z) > 0,
:::
and then extended by *analyitic continuation* to a meromorphic function that is holomorphic in the whole complex plane except zero and the negative integers, where the function has simple poles.

````

```{code-cell} sagemath
:tags: [Gamma-Function]

plot(gamma(x), (x,-5,5), ymin=-5,ymax=5,exclude = [-5..0])
```

This is the most "popular" and useful extension of the factorial function over the complex, due to the fact that
by Wielandt's Theorem we know that the Gamma function $\Gamma$ is the only holomorphic function in the right half plane such that $\Gamma(z + 1) = z \Gamma(z)$ and that is bounded in the strip $S = \{z \in \mathbb{C} \,:\, \Re z \in [1,2)\}$. Indeed other do exist, but the $\Gamma$ is the one we need for our tasks.

The next ingredient we need is a formula for performing repeated integration of
a function.

````{prf:lemma}
:label: Lemma-Swap-Integrals
If $G(x,t)$ is jointly continuous on $[c,b]\times[c,b]$:
```{math}
\int_{c}^{x} {\rm d}x_1 \int_{c}^{x_1} G(x_1,t) {\rm d}t = \int_{c}^{x} {\rm d}t \int_{t}^{x}G(x_1,t) {\rm d}x_1.
```
````

That we can use to **prove by induction** a formula for the $n$th-fold integral
of a function.

````{prf:proposition}
:label: Proposition-Formula-for-Repeated-Integrals

Given a **continuous** function $f : [c,b] \rightarrow \mathbb{R}$ we have for $x < b$
```{math}
I^n_c(x) = \int_{c}^{x}{\rm d}x_1 \int_{c}^{x_1} {\rm d}x_2 \int_{c}^{x_2} {\rm d}x_3 \cdots \int_{c}^{x_{n-1}}f(t)dt = \frac{1}{(n-1)!} \int_{c}^{x} (x-t)^{n-1} f(t){\rm d}t
```
````

````{prf:proof}
We want to use {prf:ref}`Lemma-Swap-Integrals`, for this we start from the case $n = 2$ with $G(x_1,t) \equiv f(t)$:
```{math}
      I_c^2f(x) = & \int_{c}^{x}d x_1 \int_{c}^{x_1} f(t)dt = \int_{c}^{x} f(t)\,dt \int_{t}^{x} dx_1 = \\
      = & \int_{c}^{x}(x-t) f(t)dt,
```
Then for $n = 3$
```{math}
I_c^3 f(x) = & \int_{c}^{x}d x_1 \int_{c}^{x_1} \,dx_2 \int_{c}^{x_2} f(t)dt \\ = & \int_{c}^{x}d x_1 \left[ \int_{c}^{x_1} \,dx_2 \int_{c}^{x_2} f(t)dt \right] \\
= & \int_{c}^{x}d x_1  \int_{c}^{x_1}(x_1-t) f(t)dt \\
= & \int_{c}^{x} f(t)\,dt \int_{t}^{x} (x_1 - t) dx_1 = \int_{c}^{x} f(t) \frac{(x-t)^2}{2} dt
```
Then for a generic $n$ by induction one gets:
```{math}
  I_c^n f(x) = \frac{1}{(n-1)!} \int_{c}^{x} (x-t)^{n-1} f(t)dt,
```
we have reduced the $n$--fold integral of $f$ to a single integral.
````
Nevertheless, we are still dealing with an _integer number_ $n$ of repetitions. We can now use the **Gamma function** to move from the integer case to the **real one**.

````{prf:definition}
:label: Definition-Riemann-Liouville-FractionalIntegral

Let $\Re\alpha > 0$, and let $f \in \mathbb{L}^1([a,b])$. Then for $t \in [a,b]$ we call
```{math}
  I_a^\alpha(t) = {}_{a}{D_t^{-\alpha}}f(t) = \frac{1}{\Gamma({\alpha})} \int_{0}^{t} (t-\xi)^{{\alpha} - 1} f(\xi) \,{\rm d}\xi.
```
the **Riemann--Liouville** fractional integral of $f$ of order $\alpha$, we set it to be the identity operator whenever $\alpha = 0$.
````
Let us now look at some properties before attempting the evaluation for a simple function. First of all, we have to prove that the object we have just
defined is meaningful.

::::{prf:theorem}
:label: Theorem-Riemann-Liouville-Existence
Lef $f \in \mathbb{L}^{1}[a,b]$, and $\alpha > 0$. Then, the integral $I_a^{\alpha}f(t)$ exists for almost every $t \in [a,b]$. Moreover, the function $I_a^{\alpha}f$ itself is also an element of $\mathbb{L}^{1}[a,b]$.
::::

:::{margin} Integrability of convolutions
This is a consequence of Fubini's Theorem relating integration over a set in a product space to integration in the factor spaces; see your favourite book on
Lebesgue integrals.
:::
::::{prf:proof}
It is sufficient to recognize that we can write the integral in question as
a convolution on $\mathbb{R}$, indeed:
```{math}
\int_{0}^{t} (t-\xi)^{{\alpha} - 1} f(\xi) {\rm d}\xi = \int_{-\infty}^{+\infty} \Phi_1(t-\xi) \Phi_2(\xi)\,{\rm d}\xi,
```
where
```{math}
\Phi_1(u) = \begin{cases}
u^{\alpha - 1}, & \text{ for } 0 < t \leq b-a,\\
0, \text{ otherwise},
\end{cases} \text{ and } \Phi_2(u) = \begin{cases}
f(u), & \text{ for } u \in [a,b],\\
0, & \text{ otherwise}.
\end{cases}
```
By construction both the $\Phi_j$, $j=1,2$, are in $\mathbb{L}^1(\mathbb{R})$, and thus the integral exists and is a member of $\mathbb{L}^1$ as a convolution of $\mathbb{L}^1$ functions.
::::



````{prf:proposition}
:label: Proposition-Riemann-Liouville-Semigroup
The fractional integral in {prf:ref}`Definition-Riemann-Liouville-FractionalIntegral`
operators $\{ I_c^{\alpha} : \mathbb{L}^1[a,b] \rightarrow \mathbb{L}^1[a,b], \; \alpha \geq 0 \}$ form a commutative semigroup with respect to the concatenation operation, that is
```{math}
I_c^{\alpha}(I_c^{\beta} f(t)) =  I_c^{\alpha+\beta} f(t)),\\
I_c^{\beta}(I_c^{\alpha} f(t)) =  I_c^{\alpha+\beta} f(t)).
```
The neutral element of this semigroup is the $I^0_a$ operator.
````

````{prf:proof}
We need to compute
```{math}
I_a^\alpha I_a^\beta f(x) = \frac{1}{\Gamma(\alpha)\Gamma(\beta)} \int_{a}^{x}(x-t)^{\alpha - 1} \int_{a}^{t}(t-\tau)^{\beta - 1}f(\tau)\,{\rm d}\tau\,{\rm d}t,
```
by using {prf:ref}`Theorem-Riemann-Liouville-Existence` we know that the integral exists, then
by using Fubini's theorem we can interchange the order of integration
```{math}
I_a^\alpha I_a^\beta f(x) = & \frac{1}{\Gamma(\alpha)\Gamma(\beta)} \int_{a}^{x}\int_{\tau}^{x} (x-t)^{\alpha - 1} (t-\tau)^{\beta-1} f(\tau)\,{\rm d}t
\,{\rm d}\tau  \\
= & \frac{1}{\Gamma(\alpha)\Gamma(\beta)} \int_{a}^{x} f(\tau) \int_{\tau}^{x} (x-t)^{\alpha - 1} (t-\tau)^{\beta-1} \,{\rm d}t
\,{\rm d}\tau.
```
We now use the substitution $t = \tau + s(x-\tau)$, ${\rm d}t = (x-\tau){\rm d}s$, to obtain
```{math}
I_a^\alpha I_a^\beta f(x) = & \frac{1}{\Gamma(\alpha)\Gamma(\beta)} \int_{a}^{x} f(\tau) \int_{0}^{1}[(x-\tau)(1-s)]^{\alpha - 1} [s(x-\tau)]^{\beta-1}(x-\tau) \,{\rm d}s\,{\rm d}\tau \\
= & \frac{1}{\Gamma(\alpha)\Gamma(\beta)} \int_{a}^{x} f(\tau) (x-\tau)^{m+n-1} \int_{0}^{1} (1-s)^{m-1}s^{n-1} \,{\rm d}s\,{\rm d}\tau.
```
To proceed, we now need to recognize that the second integral is indeed the **Euler $\beta$ Function**:
```{math}
  \beta(x,y) \triangleq \int_{0}^{1} u^{x-1} (1-u)^{y-1}\,du = \frac{\Gamma(x) \Gamma(y)}{\Gamma(x+y)} \quad \Re x > 0, \Re y > 0,
```
and thus
```{math}
I_a^\alpha I_a^\beta f(x) = \frac{1}{\Gamma(\alpha + \beta)}  \int_{a}^{x} (x-\tau)^{\alpha+\beta-1} f(\tau) \,{\rm d}\tau = I_a^{\alpha + \beta} f(x),
```
almost everywhere on $[a,b]$. The same proof works also if we exchange the role
of $\alpha$ and $\beta$, while we have the $0$th order operator being the neutral
element by definition.
````

````{admonition} A note on regularity.
:class: tip, dropdown

Observe that in the proof of {prf:ref}`Proposition-Riemann-Liouville-Semigroup` we could say something
more on the regularity of the resulting functions. Indeed if $f$ is a continuous
function on $[a,b]$, then also $I_a^\alpha f$ is continuous. Therefore we have
that also the concatetion $I^\alpha_a I^\beta_a$ and $I_a^{\alpha+\beta}$ are
continuous. Then what we have proved is that we have two continuos function that
are *almost everywhere* equal, and therefore they most coincide everywhere.
Furthermore, if $f \in \mathbb{L}^1[a,b]$ and $\alpha+\beta \geq 1$ we can use
{prf:ref}`Proposition-Riemann-Liouville-Semigroup` to write
```{math}
I_a^\alpha I_a^\beta f = I_a^{\alpha + \beta} f  = I_a^{\alpha + \beta - 1} I_a^{1} f, \textit{ a.e.}
```
Now, since $I_a^1 f$ is continuos, we also get that the other two way of writing
it are continuous, and thus we can conclude the equality everywhere by the same
argument as before.
````

We can now compute our first Riemann-Liouville fractiona integral.
::::{admonition} An Example of Riemann--Liouville Fractional Integral
```{math}
{}_{0}{D_t^{-{\alpha}}}t^\mu = \frac{1}{\Gamma({\alpha})} \int_{0}^{t} (t-\xi)^{{\alpha} - 1} \xi^\mu \,d\xi,
```
that should be the simplest possible example, and indeed it is as simple as using again the **Euler $\beta$ Function**:
```{math}
  \beta(x,y) \triangleq \int_{0}^{1} u^{x-1} (1-u)^{y-1}\,du = \frac{\Gamma(x) \Gamma(y)}{\Gamma(x+y)} \quad \Re x > 0, \Re y > 0.
```
To obtain it, we do the substitution for $u = \frac{\xi}{t}$, then
```{math}
{}_{0}{D_t^{-{\alpha}}}t^\mu = & \frac{t^{\alpha+\mu}}{\Gamma(\alpha)} \int_{0}^{1} u^\mu (1-u)^{\alpha-1} \,du \\ = & \frac{t^{\alpha+\mu}}{\Gamma(\alpha)} \frac{\Gamma(\mu+1)\Gamma(\alpha)}{\Gamma(\alpha+\mu+1)} \\ = &
\frac{\Gamma(\mu+1)}{\Gamma(\alpha+\mu+1)}t^{\alpha+\mu}.
```
::::

Observe that we could have computed it also with few lines of SageMath code
```{code-cell} sagemath
t,xi,alpha,mu = var('t,xi,alpha,mu')
assume(t>0)
assume(alpha>0)
assume(mu>0)
1/gamma(alpha)*integrate((t-xi)^(alpha-1)*xi^mu,xi,0,t)
```
that we can use also to look at a plot of the obtained function
```{code-cell} sagemath
def RLInt(t,alpha,mu):
     xi = var('x')
     return 1/gamma(alpha)*integrate((t-xi)^(alpha-1)*xi^mu,xi,0,t)
plot((t,RLInt(t,1.5,1.0),1/3* t^3),(t,0,5),legend_label=['t','Fractional Integral','1/3 t^3'])
```
that matches the _intuition_ of having something in-between the original function $t$ and the integral $\frac{1}{3}t^3$.

### Quadratures based on polynomial interpolation

As you have seen, calculating the integral of a simple function such as a
polynomial already requires resorting to the use of a special function such as
Euler's $\beta$. In general, thinking of being able to compute them for
arbitrary functions is a more difficult undertaking than in the case of
ordinary integrals.

However, as in the ordinary case, we can resort to quadrature formulas to obtain numerical approximations.

Let us assume that $f(t)$ in Definition {prf:ref}`Definition-Riemann-Liouville-FractionalIntegral` is suitably smooth on an interval $(a,b)$. Then we let $h = \frac{b-a}{N}$ with $N \in \mathbb{N}$, and denote with
```{math}
t_k = a + kh, \quad \text{ with } k=0,1,2,\ldots,N,
```
then we can approximate for $t = t_n$ the fractional integral as
```{math}
\left. {}_{a}D_{b}^{-\alpha} f(t)\right\rvert_{t = t_n} = & \frac{1}{\Gamma(\alpha)} \int_{a}^{t_n} (t_n - \xi)^{\alpha - 1} f(\xi)\,{\rm d}\xi \\
= & \frac{1}{\Gamma(\alpha)} \sum_{k=0}^{N-1} \int_{t_k}^{t_{k+1}} (t_k - \xi)^{\alpha-1} f(\xi)\,{\rm d}\xi.
```
We can now approximate $f(x)$ with a specific polynomial $p(x)$ on each interval, such that we can compute exactly the integrals $\int_{t_k}^{t_{k+1}} (t_k - \xi)^{\alpha-1} f(\xi)\,{\rm d}\xi$, $0 \leq k \leq N-1$. This yields quadratures by the usual look
```{math}
\left. {}_{a}D_{b}^{-\alpha} f(t)\right\rvert_{t = t_n}  \approx \sum_{k=0}^{N-1} \omega_k f(t_k),
```
that are the straightforward extension of those for integer order integrals.

The **simplest** case is obtained when we select a piecewise constant
approximation for $f(t)$ on the intervals $[t_k,t_k+1)$, $k=0,\ldots,N-1$. As
we do for the rectangular formula in the ordinary case, we select
```{math}
f(t) \approx p(t) \equiv p(t_k), \quad t \in [t_k,t_k+1), \; k=0,1,\ldots,N-1,
```
from which we get the formula
```{math}
\left. {}_{a}D_{b}^{-\alpha} f(t)\right\rvert_{t = t_n} \approx & \frac{1}{\Gamma(\alpha)} \sum_{k=0}^{N-1} f(t_k) \int_{t_k}^{t_k+1}(t_n - \xi)^{\alpha - 1}\,{\rm d}\xi   \\
= & \frac{1}{\Gamma(\alpha)} \sum_{k=0}^{N-1} f(t_k) \left[ -\frac{1}{\alpha}(t_n - \xi)^{\alpha} \right]_{t_k}^{t_{k+1}} \\
= & \sum_{k=0}^{N-1} f(t_k) \frac{1}{\alpha \Gamma(\alpha)} \left[ (t_n - t_k)^{\alpha} - (t_n - t_{k+1})^{\alpha} \right] \\
= & \sum_{k=0}^{N-1} f(t_k) \frac{1}{\alpha \Gamma(\alpha)} \left[ ( a + hn - a -kh)^{\alpha} - (a + hn -a -(k+1)h)^{\alpha} \right] \\
= & \sum_{k=0}^{N-1} f(t_k) \frac{h^\alpha}{\Gamma(\alpha+1)} \left[ (n-k)^\alpha - (n-k-1)^\alpha \right] \\
= & \sum_{k=0}^{N-1} b_{n-k-1}f(t_k),
```
where we have defined
```{math}
:label: eq:bkrectangular
b_k = \frac{h^\alpha}{\Gamma(\alpha + 1)} [ (k+1)^\alpha - k^\alpha], \qquad 0 \leq k \leq N-1.
```
With analogous computations, we can get the case in which we select the *right*
approximation
```{math}
f(t) \approx p(t) \equiv p(t_{k+1}), \quad t \in [t_k,t_k+1), \; k=0,1,\ldots,N-1,
```
and, more generally, for the *weighted* formula in which we select
```{math}
f(t) \approx p(t) \equiv \lambda p(t_{k}) + (1-\lambda)p(t_{k+1}), \quad t \in [t_k,t_k+1), \; k=0,1,\ldots,N-1, \; \lambda \in [0,1],
```
getting, for the same $b_k$ coefficients in {eq}`eq:bkrectangular`,
```{math}
\left. {}_{a}D_{b}^{-\alpha} f(t)\right\rvert_{t = t_n} \approx \sum_{k=0}^{N-1} b_{n-k-1} \left[ \lambda p(t_{k}) + (1-\lambda)p(t_{k+1}) \right], \quad \lambda \in [0,1].
```

:::{admonition} Fractional trapezoidal formula
:class: seealso

With a little patience, you can obtain for example the fractional trapezoidal
formula by choosing:
```{math}
p(t) = \frac{t_{k+1} - t}{t_{k+1}-t_{k}} f(t_k) + \frac{x - t_k}{t_{k+1}-t_k} f(t_{k+1}), \quad t \in [t_k,t_{k+1}],
```
from which one obtains
```{math}
\left. {}_{a}D_{b}^{-\alpha} f(t)\right\rvert_{t = t_n} \approx \sum_{k=0}^{N} a_{k,N} f(t_k),
```
where
```{math}
a_{k,N} = \frac{h^\alpha}{\Gamma(\alpha + 2)} \begin{cases}
(N-1)^{\alpha + 1} - (N-1-\alpha)N^\alpha, & k = 0,\\
(N-k+1)^{\alpha+1} + (N-1-k)^{\alpha+1} - 2(N-k)^{\alpha + 1}, & 1 \leq k \leq N-1,\\
1, & k = N.
\end{cases}
```
:::

## Riemann-Liouville Derivative

## Bibliography

```{bibliography}
:filter: docname in docnames
```
