#set text(
 font: "Noto Sans SignWriting Regular",
)
#set math.equation(numbering: "1", supplement: [Eq.])

//#show ref: it => {
//  let eq = math.equation
//  let el = it.element
//  if el != none and el.func() == eq {
//    // Override equation references.
//    numbering(
//      el.numbering,
//      ..counter(eq).at(el.location())
//    )
//  } else {
//    // Other references as usual.
//    it
//  }
//}
//
= One dimensional Square Well

== Even potential 

- Let $V(-x) = V(x)$, if $phi(x)$ is a solution, then $phi(-x)$ is also a solution

- For the eigenvalue $E$, there always exists a set of solutions that are *complete*, and each eigenfunction holds a specific parity.

- Solutions to the eigenvalue $E$ are indegenerate

== General properties in 1 dimension

- for two degenerate solutions $psi_1$ and $psi_2$ belonging to the same eigenvalue, $psi_1 psi_2^' - psi_1^' psi_2 = C$
$ [-planck.reduce^2 / (2 m) d^2/(d x^2) + V] psi_1 = V psi_1 $ <1>
$ [-planck.reduce^2 / (2 m) d^2/(d x^2) + V] psi_2 = V psi_2 $ <2>
$phi_1 times$ @2 $-$ $phi_2 times$ @1:
$ psi_1 psi_2^' - psi_1^' psi_2 = (psi_1 psi_2^' - psi_1^' psi_2)^' = 0 space qed $

- for a poential without singularity, bound states are indegenerate
$ psi_1 psi_2^' - psi_1^' psi_2 = C $
for the bound state, $limits(psi)_(x -> oo) = 0$, $C = 0$, $(psi_1^')/psi_1 = (psi_2^')/psi_2$, so they are the same state

== Infinite Square Well
$ V(x) = cases(
  0"," space 0 lt.eq x lt.eq a,
  oo"," space "otherwise"
  ) $
boundary conditions
$ psi(0) = 0\
  psi(a) = 0 $
Indise well, 
$ psi(x) = A e^(i/planck.reduce p x) + B e^(-i/planck.reduce p x)\
  cases(
    A + B = 0,
    A e^(i/planck.reduce p a) + B e^(-i/planck.reduce p a) = 0
    )
$
$ B = -A\
  sin p/planck.reduce a = 0
$
- $psi_n(x) = A_n sin (n pi)/a x, space E_n = (n^2 planck.reduce^2 pi^2)/(2 m a^2)$
$ integral_0^a abs(A_n)^2 sin^2 (n pi)/a x d x &= abs(A_n)^2 / 2 integral_0^a 1 - cos (2 n pi)/a x d x\
  &= (abs(A_n)^2 a)/2 = 1 $
- $E_n prop n^2$
$
  Delta E_n prop n\
  (Delta E_n)/(E_n) prop 1/n
$
- $psi_n prop sin (n pi)/a x$

- ground state $E_1, psi_1$

- $psi(x) = sum a_n psi_n (x)$
- density of states
$ rho(E) = (delta N)/(delta E) = 1/(d E slash d n)\
  E = A n^2\
  (d E)/(d n) = 2 A n\
  rho(E) = 1/(2 sqrt(E A)) = (a sqrt(2 m))/(2 pi planck.reduce sqrt(E))
  $
- $psi(x, t_0) = sqrt(2/a) sin pi/a x$


//Inside the well, time-independent Schrodinger equation
//$ -planck.reduce^2 /(2m)(d^2 psi)/(d x^2) = E psi $
//let $k eq.triple sqrt(2m E)/planck.reduce$, we have
//$ (d^2 psi)/(d x^2) = -k^2 psi $
//where the general solution is
//$ psi(x) = A sin k x + B cos k x $ <sol1>
//_A_ and _B_ are, of course, *complex constants*. As for the continuity,
//$ psi(0) = psi(a) = 0 $
//plugging into @sol1,
//$ A sin 0 + B cos 0 = B = 0 $
//so $psi(x) = A sin k x $. And again
//$ A sin k a = 0 $
//which means either $A=0$ (discarded), or $sin k a = 0$, i.e.
//$ k a = 0, plus.minus pi, plus.minus 2pi, plus.minus 3pi, dots.h $
//So here comes the quantum number, let's say, $n$, we have
//$ k_n = (n pi)/a, space "with" n in II $
//and the energy
//$ E_n = (planck.reduce^2 k_n^2)/(2m) = (n^2 pi^2 planck.reduce^2)/(2m a^2) $
//at last, the normalization
//$ integral_0^a abs(A)^2 sin^2 (k x) d  x = abs(A)^2 a/2 =  1, space space "so" abs(A)^2 = 2/a $
//and the time-independent solutions are
//$ psi_n (x) = sqrt(2/a) sin((n pi)/a x) $

== Harmonic Oscillator
$ V(x) = 1/2 m omega^2 x^2 $

=== Algebraic method
$ hat(H) = 1/(2m)[hat(p)^2 + (m omega x)^2] $
$ hat(a)_(plus.minus) eq.triple 1/sqrt(2 planck.reduce m omega)(minus.plus i hat(p) + m omega x)  $
$ hat(a)_minus hat(a)_plus &= 1/(2 planck.reduce m omega)(i hat(p) + m omega x)(-i hat(p)+m omega x)\ &= 1/(2planck.reduce m omega)[hat(p)^2 + (m omega x)^2 - i m omega(x hat(p) - hat(p) x)]\ &= 1/(2planck.reduce m omega)[hat(p)^2 + (m omega x)^2]-i/(2planck.reduce)[x, hat(p)] $
$ [x, hat(p)] f(x) &= [x(-i planck.reduce)d/(d x)f - (-i planck.reduce)d/(d x)(x f)] = -i planck.reduce(x (d f)/(d x) - x(d f)/(d x) - f)\ &= i planck.reduce f(x) $
$ [x, hat(p)] = i planck.reduce $
$ hat(a)_minus hat(a)_plus = 1/(planck.reduce omega)hat(H) + 1/2 $
$ hat(H) = planck.reduce omega (hat(a)_- hat(a)_+ - 1/2) $
$ hat(a)_plus hat(a)_minus = 1/(planck.reduce omega)hat(H) - 1/2 $
$ hat(H) = planck.reduce omega (hat(a)_+ hat(a)_- + 1/2) $
$ hat(H) = planck.reduce omega (hat(a)_plus.minus hat(a)_minus.plus plus.minus 1/2) $ <ham>
$ hat(H)(hat(a)_+ psi) &= planck.reduce omega(hat(a)_plus hat(a)_minus + 1/2)(hat(a)_plus psi) = planck.reduce omega(hat(a)_+ hat(a)_- hat(a)_+ + 1/2 hat(a)_+)psi\ &= planck.reduce omega hat(a)_+ (hat(a)_- hat(a)_+ + 1/2)psi = hat(a)_+ [planck.reduce omega (hat(a)_+ hat(a)_- + 1 + 1/2)psi]\ &= hat(a)_+ (hat(H) + planck.reduce omega)psi = hat(a)_+ (E + planck.reduce omega)psi = (E + planck.reduce omega)(hat(a)_+ psi) $
$ hat(H)(hat(a)_- psi) &= planck.reduce omega(hat(a)_- hat(a)_+ - 1/2)(hat(a)_- psi) = planck.reduce omega hat(a)_- (hat(a)_+ hat(a)_- - 1/2)psi\ &= hat(a)_- [planck.reduce omega (hat(a)_- hat(a)_+ -1 - 1/2)psi]\ &= hat(a)_- (hat(H) - planck.reduce omega)psi = hat(a)_- (E - planck.reduce omega)psi\ &= (E - planck.reduce omega)(hat(a)_- psi) $
find the ground state with lowest energy such that
$ hat(a)_- psi_0 = 0 $
i.e.
$ 1/sqrt(2 planck.reduce m omega) (planck.reduce d/(d x) + m omega x) psi_0 = 0\ (d psi_0)/(d x) = -(m omega)/(planck.reduce)x psi_0\ psi_0 (x) = A e^(-(m omega)/(2 planck.reduce)x^2) $
to normalize it
$ abs(A)^2 integral_(-oo)^oo e^(-m omega x^2 slash planck.reduce) d x = abs(A)^2 sqrt((pi planck.reduce)/(m omega)) = 1 $
thus $A^2 = sqrt(m omega slash pi planck.reduce)$, so 
$ psi_0 (x) = ((m omega)/(pi planck.reduce))^(1 slash 4) e^(-(m omega)/(2 planck.reduce)x^2) $
so the interesting part here is that we don't exactly know if $A$ is positive or negative, but rather the $abs(A)^2$, this may be the philosophy of quantum mechanics: we cannot determine the wave function, but the probability distribution, i.e. $abs(psi)^2$

Anyway, with the ground state and ladder operator, we can obtain all eigenstates:
$ psi_n (x) = A_n (hat(a)_+)^n psi_0 (x), space "with" E_n = (n + 1/2)planck.reduce omega $ <states>
now we are going to calculate the normalization constant ($A_n$) algebraically, we have
$ hat(a)_+ psi_n = c_n psi_(n+1), space hat(a)_- psi_n = d_n psi_(n -1) $ <prop>
and operator $hat(a)_plus.minus$ is hermitian conjugate of $hat(a)_minus.plus$:
$ angle.l f bar.v hat(a)_plus.minus g angle.r &= 1/sqrt(2planck.reduce m omega) integral_(-oo)^oo f^* (minus.plus planck.reduce d/(d x) + m omega x)g d x\
&= 1/sqrt(2planck.reduce m omega) integral_(-oo)^oo  (m omega x f)^* g d x + 1/sqrt(2planck.reduce m omega) integral_(-oo)^oo f^* (minus.plus)planck.reduce (d g)/(d x) d x \
&= 1/sqrt(2planck.reduce m omega) integral_(-oo)^oo  (m omega x f)^* g d x + 1/sqrt(2planck.reduce m omega) f^* (minus.plus) planck.reduce g bar.v_(-oo)^oo - 1/sqrt(2 planck.reduce m omega) integral_(-oo)^oo (minus.plus)planck.reduce g d f^* \
&= 1/sqrt(2planck.reduce m omega)integral_(-oo)^oo [(plus.minus planck.reduce d/(d x) + m omega x)f]^* g d x \
&= angle.l hat(a)_minus.plus f bar.v g angle.r
$
invoking @ham and @states we have
$ planck.reduce omega (hat(a)_plus.minus hat(a)_minus.plus plus.minus 1/2) psi_n = [(n + 1/2)planck.reduce omega] psi_n\
hat(a)_+ hat(a)_- psi_n = n psi_n , space hat(a)_- hat(a)_+ psi_n = (n + 1)psi_n $ <32>
to utilize the hermitian conjugate, consider the integral
$ integral_(-oo)^oo (hat(a)_plus.minus psi_n)^* (hat(a)_plus.minus psi_n) d x = integral_(-oo)^oo (hat(a)_minus.plus hat(a)_plus.minus psi_n)^* psi_n d x $
with @32 and @prop
$ &integral_(-oo)^oo (hat(a)_+ psi_n)^* (hat(a)_+ psi_n) d x = abs(c_n)^2 integral abs(psi_(n+1))^2 d x = (n + 1)integral_(-oo)^oo abs(psi_n)^2 d x\
&integral_(-oo)^oo (hat(a)_- psi_n)^* (hat(a)_- psi_n)d x = abs(d_n)^2 integral abs(psi_(n-1))^2 d x = n integral_(oo)^oo abs(psi_n) ^2 d x $
but hey, we haven't make any restrictions on $psi$ this far, since we come from @prop, now we assume that $psi_n$ and $psi_(n plus.minus 1)$ are normalized, and we are finding the relation between the normalized constants, it follows that $ abs(c_n)^2  = n+1$ and $abs(d_n)^2 = n$, and hence
$ hat(a)_+ psi_n = sqrt(n+1)psi_(n+1) , space hat(a)_- psi_n = sqrt(n) psi_(n - 1) $
once again we still don't know whether to pick the positive or the negative sign LOL, whatever, finally here it comes
$ psi_n = 1/sqrt(n!)(hat(a)_+)^n psi_n $
Orthogonality
$ angle.l psi_m bar.v hat(a)_plus hat(a)_- psi_n angle.r &= n angle.l psi_m bar.v psi_n angle.r\
  &= angle.l hat(a)_- psi_m bar.v hat(a)_- psi_n angle.r = angle.l hat(a)_+ hat(a)_- psi_m bar.v psi_n angle.r\
  &= m angle.l psi_m bar.v psi_n angle.r
$
$ angle.l psi_m bar.v psi_n angle.r = delta_(m n) $

=== Analytic Method
Let's, at least, write the Schrodinger equation once
$ -planck.reduce^2/(2m)(d^2 psi)/(d x^2) + 1/2 m omega^2 x^2 psi = E psi $
and someone, somehow, introduced two variables
$ xi eq.triple sqrt((m omega)/planck.reduce)x, space K eq.triple (2E)/(planck.reduce omega) $
and how we got to know that we need to do this though, I have no idea, whatever, than the equation reads
$ (d^2 psi)/(d xi^2) = (xi^2 - K)psi $ <37>
when $xi$ goes infinity with constant energy
$ (d^2 psi)/(d xi^2) approx xi^2 psi $
so we have this
$ psi(xi) approx A e^(-xi^2 slash 2) + B e^(+xi^2 slash 2) $
clearly $B e^(+xi^2 slash 2)$ results in diverging, so $B$ must be 0; now, quite like the separation of variables, we separate the asymptotic part like this
$ psi(xi) = h(xi) e^(-xi^2 slash 2) $
in hopes that $h(xi)$ has a simpler form.
$ (d psi)/(d xi) = ((d h)/(d xi) - xi h) e^(-xi^2 slash 2)\
  (d^2 psi)/(d xi^2) = ((d^2 h)/(d xi^2) - 2xi (d h)/(d xi) + (xi^2 - 1)h)e^(-xi^2 slash 2) $
@37 than becomes _*Hermite equation*_
$ (d^2 h)/(d xi^2) - 2xi (d h)/(d xi)+ (K-1)h = 0 $ <46>
find the power series solution in terms of $xi$
$ h(xi) = sum_(j=0)^oo a_j xi^j\
  (d h)/(d xi) = sum_(j = 0)^oo (j+1) a_(j+1) xi^(j)\
  (d^2 h)/(d xi^2) = sum_(j=0)^oo (j+1)(j+2)a_(j+2)xi^j $
putting back into @46
$ sum_(j=0)^oo (j+1)(j+2)a_(j+2)xi^j - 2sum_(j=0)^oo (j+1)a_(j+1) xi^(j+1) + (K-1)sum_(j=0)^oo a_j xi^j = 0\
sum_(j = 0)^oo [(j+1)(j+2)a_(j+2) - 2 j a_j + (K-1)a_j]xi^j = 0\
a_(j+2) = (2j+1-K)/((j+1)(j+2))  a_j $ <48>
the upper asympotic behavior of this recursion is
$ a_(j+2) approx 2/j a_j\
  a_j approx C/((j slash 2)!)\
  h(xi) approx C sum 1/((j slash 2)!) xi^j approx C sum 1/(j!) xi^(2j)  approx C e^(xi^2) $
so this recursion must terminates, i.e. 
$ exists n, space "s.t." space K = 2n+1 $
again we derived this energy with a different method
$ E_n = (n+1/2)planck.reduce omega,  space n = 0, 1, 2, dots.h  $
the recursion formula @48 is now
$ a_(j+2)  =  (-2(n-j))/((j+1)(j+2))a_j $
$ psi_n(x) = ((m omega)/(pi planck.reduce))^(1/4) 1/sqrt(2^n n!) H_n (xi)e^(-xi^2 / 2) $
First few Hermite polynomials
$ &H_0(x) = 1\
  &H_1(x) = 2x\
  &H_2(x) = 4x^2 - 2\
  &H_3(x) = 8x^3 - 12x $
According to rodrigues formula
$ H_n (x) = (-1)^n e^(x^2) (d/(d x))^n e^(-x^2) $
generating function
$ e^(-t^2 + 2t x) = sum_(n = 0)^oo H_n (x) t^n / (n!) $
$ angle.l H_m bar.v e^(-x^2) H_n angle.r &= integral_(-oo)^oo (-1)^m (d/(d x))^m e^(-x^2) H_n d x\
  &= integral_(-oo)^oo (-1)^n H_n d (d/(d x))^(m - 1)e^(-x^2)\
  &= (-1)^n H_n (d/(d x))^(m - 1)e^(-x^2) bar.v_(-oo)^oo + integral_(-oo)^oo (-1)^(m-1) (d/(d x))^(m - 1) e^(-x^2) d H_n\
  &= integral_(-oo)^oo (-1)^(m-2) (d/(d x))^(m-2) e^(-x^2) d (d/(d x)) H_n\
  &= integral_(-oo)^oo (-1)^(m - k) (d/(d x))^(m - k) e^(-x^2) (d/(d x))^k H_n d x $
when $m = n = k$
$ angle.l H_m bar.v e^(-x^2) H_n angle.r &= integral_(-oo)^oo e^(-x^2) n! 2^n d x\
  &= 2^n sqrt(pi) n! $
so for the general case
$ angle.l H_m bar.v e^(-x^2) H_n angle.r = 2^n sqrt(pi) n! delta_(m n) $
$ (d H_n)/(d x) &= (-1)^n 2 x e^(x^2) (d/(d x))^n e^(-x^2) + (-1)^n e^(x^2) (d/(d x))^(n+1) e^(-x^2)\
  &= 2x H_n(x) - H_(n+1)(x) $

== Free Particle
$ - planck.reduce^2/(2 m) (d^2 psi)/(d x^2) = E psi\
  (d^2 psi)/(d x^2) = -k^2 psi, space k eq.triple sqrt(2 m E)/planck.reduce\
  psi(x) = A e^(i k x) + b e^(-i k x)\
  $
$ Psi(x, t) = A e^(i k(x - (planck.reduce k)/(2 m)t)) + B e^(-i k (x + (planck.reduce k)/ (2m))) $
a set of orthogonal eigenfunctions
$ Psi_k (x, t) = A e^(i k ( x - (planck.reduce k)/(2 m)t)) $
linear combination
$ Psi = 1/sqrt(2 pi)integral_(-oo)^oo phi.alt(k)e^(i k (x - (planck.reduce k)/(2 m)t)) d k $
$ cases(
  Psi(x, 0) = 1/sqrt(2 pi) integral_(-oo)^oo phi.alt(k)e^(i k x) d x"," space "position space",
  phi.alt(k) = 1/sqrt(2 pi)integral_(-oo)^oo Psi(x, 0) e^(-i k x) d x"," space "momentum space"
  )\
  T = (4 m pi)/(planck.reduce k^2)\
  lambda = (2 pi)/k $
$ hat(p)Psi_k = -i planck.reduce d/(d x) Psi_k = p A e^(i k(x - (planck.reduce k)/(2 m)t)) = planck.reduce k A e^(i k(x - (planck.reduce k)/(2 m)t))\
  p = planck.reduce k $
this is essentially the *de Broglie formula*

=== wave packet

$ Psi = 1/sqrt(2 pi) integral_(-oo)^oo phi.alt(k)e^(i k x - omega t) d k $
to get the classical velocity, consider a wave packet centered around some momentum $k_0$, so $phi.alt$ is negligible except in the vicinity of $k_0$, we can expand $omega (k)$ into taylor series and keep the leading terms:
$ omega(k) approx omega_0 + omega_0^' (k - k_0) $
to center the integral at $k_0$, let $s := k - k_0$
$ Psi & approx 1/(sqrt(2 pi)) integral_(-oo)^oo phi.alt(k_0+s)e^(i[(k_0 + s)x - (omega_0 + omega_0^' s)t]) d s\
  &= 1/sqrt(2 pi)e^(i (k_0 x - omega_0 t))integral_(-oo)^oo phi.alt(k_0 + s)e^(i s (x - omega_0^' t)) d s $
the phase velocity
$ v_"phase" = omega/k |_(k = k_0) $
and a well-defined group velocity (independent of $k$ inside the integral):
$ v_"group" = (d omega)/(d k) |_(k = k_0) $
for the free particle, $omega = (planck.reduce k^2 slash 2 m)$, $v_"classical" = v_"group" = 2 v_"phase"$
