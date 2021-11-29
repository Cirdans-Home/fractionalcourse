# An introduction to fractional calculus: fundamental ideas and numerics

Fractional calculus, i.e., the definition of integrodifferential operators of
non-integer order, has a long history. It runs parallel to that of modern
differential calculus as it started with the first correspondence between
Leibniz, John Bernoulli, and De l'Hôpital in which it appears as a pure
formula-generating curiosity far removed from the geometric interpretation of
ordinary derivatives. Then, it evolved towards more modern works by Riemann,
Hardy, and Littlewood, finally reaching the contemporary efforts to extend and
interpret its objects {cite}`MillerRoss,MR0361633`. Where the possibility of
doing rigorous mathematics and proving results within the new theory may be
reason enough to invest time in its study, the push of the last few years in
this endeavor is due, as often happens, to the birth of different areas of
application outside of mathematics. Consider, for example,
the modeling of the transport of fluids within highly non-homogeneous media, the
modeling of non-local phenomena, or the search for evolutionary models of a
non-Markovian nature. These needs have led to the definition of new operators or
to the reinterpretation of classical instruments and, at the same time, to the
search for numerical methods to simulate these models efficiently and accurately.

The purpose of this course is, on the one hand, to introduce the theoretical and
analytical tools necessary to formulate the constitutive laws of these new models
(definitions of operators, existence and uniqueness theorems, regularity of the
solution, etc.), from other, to show how it is possible to construct numerical
methods for their simulation. In particular, we will focus on the definition of
fractional integral operators, with the associated fractional derivative
operations according to Riemann-Liouville, Caputo, and Riesz. With these
instruments we will first discuss the formulation of ordinary fractional
equations and their solution {cite}`DiethelmBook`, then their use in the
definition of fractional diffusion problems as a fractional partial differential
equation {cite}`PodlubnyBook`. In the final part, we will briefly touch on some
questions concerning the definition of the fractional Laplacian and some of its
applications within modeling for complex networks {cite}`ArrigoDurastante,cnaa017,RiascosBook`.


## How to obtain these notes and run the code

This material has been collected in order to build a *Ph.D. course* on
numerical methods related to Fractional Calculus. The material is
already divided into **10 lectures** and combines together theoretical
and implementative aspects that one has to face in solving problems
related to fractional derivatives and their applications.


:::{note}
The numerical example in here are a mix of **SageMath** and **Python**, so to
run them you'll need both. There is the possibility of using the *Binder*
service online to run them.
:::

You can obtain everything to compile this notes from [my GitHub account](https://github.com/Cirdans-Home):
```bash
git clone https://github.com/Cirdans-Home/fractionalcourse
```
To **actually compile** it you'll need a working installation of the Jupyter-book
toolchain together with SageMath, and, if you want to produce a LaTeX version
the TeXLive toolchain to produce and compile everything as `.pdf`.

.. only:: latex

## Bibliography

There exist many complete and good books covering all the material discussed
here, from classical text like {cite}`MR0361633,MillerRoss,PodlubnyBook`, to
more recent works {cite}`DiethelmBook,JinBook`. A classic problem in this sector
is the proliferation of notations and attributions of objects to sequences by
different mathematicians. The material contained here draws on several of these
sources but makes an effort to translate them into a single notation. Notation
which obviously is far from being a non-existent accepted standard, but which in
the eyes of the writer is sufficiently simple.

```{bibliography}
:filter: docname in docnames
```
