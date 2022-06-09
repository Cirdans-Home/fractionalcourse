# An introduction to fractional calculus: fundamental ideas and numerics

Fractional calculus, i.e., the definition of integrodifferential operators of non-integer order, has a long
history. It runs parallel to that of modern differential calculus as it started with the first correspondence
between Leibniz, John Bernoulli, and De l’Hôpital in which it appears as a pure formula-generating
curiosity far removed from the geometric interpretation of ordinary derivatives. Then, it evolved towards
more modern works by Riemann, Hardy, and Littlewood, finally reaching the contemporary efforts
to extend and interpret its objects [5, 6]. Where the possibility of doing rigorous mathematics and
proving results within the new theory may be reason enough to invest time in its study, the push of the
last few years in this endeavor is due, as often happens, to the birth of different areas of application
outside of mathematics. Consider, for example, the modeling of the transport of fluids within highly
non-homogeneous media, the modeling of non-local phenomena, or the search for evolutionary models of
a non-Markovian nature. These needs have led to the definition of new operators or to the reinterpretation
of classical instruments and, at the same time, to the search for numerical methods to simulate these
models efficiently and accurately.
The purpose of this course is, on the one hand, to introduce the theoretical and analytical tools
necessary to formulate the constitutive laws of these new models (definitions of operators, existence and
uniqueness theorems, regularity of the solution, etc.), from other, to show how it is possible to construct
numerical methods for their simulation. In particular, we will focus on the definition of fractional integral
operators, with the associated fractional derivative operations according to Riemann-Liouville, Caputo,
and Riesz. With these instruments we will first discuss the formulation of ordinary fractional equations
and their solution [3], then their use in the definition of fractional diffusion problems as a fractional
partial differential equation [7]. Regarding numerical methods, we will focus on finite
difference and multi-step methods for simulating initial value problems with time-fractional
derivatives. Regarding the solution of boundary value problems with space-fractional derivatives,
we will focus on the use of finite difference methods and on the connection between
the properties of the resulting discretization matrices and the analytic
properties of the operators.

In the final part, we will briefly touch on some questions concerning
the definition of the fractional Laplacian and some of its applications within modeling for complex
networks [1, 2, 4].

## Bibliography

1. F. Arrigo and F. Durastante. “Mittag–Leffler Functions and their Applications in Network Science”. In: SIAM J.
Matrix Anal. Appl. 42.4 (2021), pp. 1581–1601. issn: 0895-4798. doi: 10.1137/21M1407276.
2. M. Benzi, D. Bertaccini, F. Durastante, and I. Simunec. “Non-local network dynamics via fractional graph
Laplacians”. In: Journal of Complex Networks 8.3 (Aug. 2020). cnaa017. issn: 2051-1329.
3. K. Diethelm. The analysis of fractional differential equations. Vol. 2004. Lecture Notes in Mathematics. An
application-oriented exposition using differential operators of Caputo type. Springer-Verlag, Berlin, 2010,
pp. viii+247. isbn: 978-3-642-14573-5. doi: 10.1007/978-3-642-14574-2.
4. T. Michelitsch, A. P. Riascos, B. Collet, A. Nowakowski, and F. Nicolleau. Fractional Dynamics on Networks and
Lattices. John Wiley & Sons, Ltd, 2019, pp. xvi+305. isbn: 9781119608165. doi: 10.1002/9781119608165.
5. K. S. Miller and B. Ross. An introduction to the fractional calculus and fractional differential equations. A
Wiley-Interscience Publication. John Wiley & Sons, Inc., New York, 1993, pp. xvi+366. isbn: 0-471-58884-9.
6. K. B. Oldham and J. Spanier. The fractional calculus. Mathematics in Science and Engineering, Vol. 111. Theory and
applications of differentiation and integration to arbitrary order, With an annotated chronological bibliography
by Bertram Ross. Academic Press [Harcourt Brace Jovanovich, Publishers], New York-London, 1974, pp. xiii+234.
7. I. Podlubny. Fractional differential equations. Vol. 198. Mathematics in Science and Engineering. An introduction
to fractional derivatives, fractional differential equations, to methods of their solution and some of their
applications. Academic Press, Inc., San Diego, CA, 1999, pp. xxiv+340. isbn: 0-12-558840-2.
