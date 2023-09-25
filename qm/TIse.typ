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
== Infinite Square Well
$ V(x) = cases(
  0"," space 0 lt.eq x lt.eq a,
  oo"," space "otherwise"
  ) $
Inside the well, time-independent Schrodinger equation
$ -planck.reduce^2 /(2m)(d^2 psi)/(d x^2) = E psi $
let $k eq.triple sqrt(2m E)/planck.reduce$, we have
$ (d^2 psi)/(d x^2) = -k^2 psi $
where the general solution is
$ psi(x) = A sin k x + B cos k x $ <sol1>
_A_ and _B_ are, of course, *complex constants*. As for the continuity,
$ psi(0) = psi(a) = 0 $
plugging into @sol1,
$ A sin 0 + B cos 0 = B = 0 $
so $psi(x) = A sin k x $. And again
$ A sin k a = 0 $
which means either $A=0$ (discarded), or $sin k a = 0$, i.e.
$ k a = 0, plus.minus pi, plus.minus 2pi, plus.minus 3pi, dots.h $
So here comes the quantum number, let's say, $n$, we have
$ k_n = (n pi)/a, space "with" n in II $
and the energy
$ E_n = (planck.reduce^2 k_n^2)/(2m) = (n^2 pi^2 planck.reduce^2)/(2m a^2) $
at last, the normalization
$ integral_0^a abs(A)^2 sin^2 (k x) d  x = abs(A)^2 a/2 =  1, space space "so" abs(A)^2 = 2/a $
and the time-independent solutions are

== Harmonic Oscillator
consider the potential 
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
once again we still don't know whether to pick the positive or the negative sign LOL, but, whatever, finally here it comes
$ psi_n = 1/sqrt(n!)(hat(a)_+)^n psi_n $
are they orthogonal?
$ angle.l psi_m bar.v hat(a)_plus hat(a)_- psi_n angle.r &= n angle.l psi_m bar.v psi_n angle.r\
  &= angle.l hat(a)_- psi_m bar.v hat(a)_- psi_n angle.r = angle.l hat(a)_+ hat(a)_- psi_m bar.v psi_n angle.r\
  &= m angle.l psi_m bar.v psi_n angle.r
$
$ angle.l psi_m bar.v psi_n angle.r = delta_(m n) $
