#set text(
 font: "Noto Sans SignWriting Regular",
)
#set math.equation(numbering: "1", supplement: [Eq.])
= Formalism
- An operator $Q$ is hermitian operator if
$ angle.l f bar.v hat(Q) g angle.r = angle.l hat(Q) f bar.v g angle.r space "for all" space f space "and" space g $
- *Observables are represented by hermitian operators*

- hermitian conjugate of $hat(Q)$ satisfies
$ angle.l f | hat(Q) g angle.r = angle.l hat(Q)^dagger f | g angle.r $
- hermitian operator is equal to its hermitian conjugate

- $(hat(Q)hat(R))^dagger = hat(R)^dagger hat(Q)^dagger$
$ angle.l f | hat(Q) hat(R) g angle.r = angle.l hat(Q) f | hat(R) g angle.r = angle.l hat(R)hat(Q)f | g angle.r $
- the product of two hermitian operators is hermitian if and only if they are commutative
$ angle.l f bar.v hat(Q)hat(R) g angle.r = angle.l (hat(Q)hat(R))^dagger f bar.v g angle.r\
(hat(Q)hat(R))^dagger = hat(R)^dagger hat(Q)^dagger\
hat(Q)hat(R) =  hat(R)hat(Q)
  $
- $x^dagger = x$

- $(d/(d x))^dagger  = - d/(d x)$
$ angle.l f | d/(d x) g angle.r &= integral f d/(d x) g d x\ &= integral d f g - integral g d/(d x) f d x\ &= angle.l (-d/(d x)) f | g angle.r $
cause we are in hilbert space, $integral_(-oo)^oo d f g$ vanishes

== Eigenfunctions of hermitian operators
=== Discrete Spectra
- eigenvalues are real
$ hat(Q) f = q f\
  angle.l f | hat(Q)f angle.r = angle.l hat(Q)f | f angle.r\
  q angle.l f | f angle.r  = q^* angle.l f | f angle.r\
  q = q^* space qed $
- eigenfunctions belonging to distinct eigenvalues are orthogonal
$ hat(Q) f = q f, space hat(Q)g = q' g\
  angle.l f | hat(Q)g angle.r  = angle.l hat(Q)f | g angle.r\
  q' angle.l f | g angle.r = q^* angle.l f | g angle.r\
  "since" space q' eq.not q, space angle.l f | g angle.r = 0 space qed $
while this only applies to indegenerate states, we can always use *Gram-Schmidt orthogonalization prodedure* to construct orthogonal eigenfunctions within each degenerate subspace
- eigenfunctions of an observable operator are _*complete*_
=== Continuous Spectra
- eigenfunctions are non-normalizable

- eigenfunctions with real eigenvalues are Dirac orthonormalizable and complete
more specifically, for momentum operator
$ -i planck.reduce d/(d x) f = p f\
  f = A e^(i p x slash planck.reduce)\
  integral_(-oo)^oo f_(p')^* f_p d x = |A|^2 integral_(-oo)^oo e^(i(p - p')x slash planck.reduce) d x = |A|^2 2 pi planck.reduce delta (p - p')\
  A = 1/sqrt(2 pi planck.reduce)\
  angle.l f_(p') | f_p angle.r = delta (p-p') $
position operator
$ hat(x) g_y (x) = y g_y (x)\
  g_y(x) = delta(x - y)\
  angle.l g_y | g_(y') angle.r = delta(y - y') $

== Uncertainty Principle
$ sigma^2 = angle.l (Q - angle.l Q angle.r)^2 angle.r = angle.l Psi | (hat(Q) - q)^2 Psi angle.r $
for observables A and B,
$ sigma_A^2 = angle.l (hat(A) - angle.l A angle.r) Psi | (hat(A) - angle.l A angle.r) Psi angle.r = angle.l f|f angle.r\
  sigma_B^2 = angle.l (hat(B) - angle.l B angle.r) Psi | (hat(B) - angle.l B angle.r) Psi angle.r = angle.l g|g angle.r $
where $f := (hat(A) - angle.l A angle.r)Psi, space g := (hat(B) - angle.l B angle.r)Psi$
$ sigma_A^2 sigma_B^2 = angle.l f|f angle.r angle.l g|g angle.r gt.eq |angle.l f|g angle.r|^2 space "(Schwarz innequality)"\
  |z|^2 = Re^2(Z) + Im^2(z) gt.eq Im^2(z) = [1/(2 i)](z - z^*)^2\
  sigma_A^2 sigma_B^2 gt.eq (1/(2 i)[angle.l f|g angle.r - angle.l g|f angle.r])^2 $
$ angle.l f|g angle.r &= angle.l (hat(A) - angle.l A angle.r)Psi|(hat(B) - angle.l B angle.r)Psi angle.r = angle.l Psi|(hat(A) - angle.l A angle.r)(hat(B) - angle.l B angle.r)Psi angle.r\
  &= angle.l Psi|(hat(A)hat(B) - hat(A) angle.l B angle.r - hat(B)angle.l A angle.r + angle.l A angle.r angle.l B angle.r)Psi angle.r\
  &= angle.l Psi|hat(A)hat(B)Psi angle.r - angle.l B angle.r angle.l Psi|hat(A)Psi angle.r - angle.l A angle.r angle.l Psi|hat(B)Psi angle.r + angle.l A angle.r angle.l B angle.r angle.l Psi|Psi angle.r\
  &= angle.l hat(A)hat(B) angle.r - angle.l A angle.r angle.l B angle.r $
$ angle.l g|f angle.r = angle.l hat(B)hat(A) angle.r - angle.l A angle.r angle.l B angle.r $
$ angle.l f|g angle.r - angle.l g|f angle.r = angle.l hat(A)hat(B) angle.r - angle.l hat(B)hat(A) angle.r = angle.l [hat(A), hat(B)] angle.r $
- $sigma_A^2 sigma_B^2 gt.eq (1/(2 i)angle.l [hat(A), hat(B)]angle.r)^2 $

- $sigma_x^2 sigma_p^2 gt.eq (planck.reduce/2)^2$
== Comutator

- $[hat(A) + hat(B), hat(C)] = [hat(A), hat(C)] + [hat(B), hat(C)]$

- $[hat(A)hat(B), hat(C)] = hat(A)[hat(B), hat(C)] + [hat(A), hat(C)]hat(B)$
$ (hat(A)hat(B)hat(C) - hat(C)hat(A)hat(B))psi\
  (hat(A)(hat(B)hat(C)-hat(C)hat(B)) + (hat(A)hat(C) - hat(C)hat(A)) hat(B))psi $
- $[x^n, hat(p)] = i planck.reduce n x^(n-1)$
$ (-x^n i planck.reduce d/(d x)  + i planck.reduce d/(d x) x^n) psi &= -x^n i planck.reduce d/(d x) psi + i planck.reduce n x^(n-1) x psi + i planck.reduce x^n d/(d x) psi\ &= i planck.reduce n x^(n-1) psi $
- $[f(x), hat(p)] = i planck.reduce (d f)/(d x)$

- two commuting operators can share a same set of complete functions
if $hat(P)$ and $hat(Q)$ have a complete set of common eigenfunctions, then $[hat(P), hat(Q)]f=0$ holds for any function in hilbert space
