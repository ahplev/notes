 #set text(
   font: "Noto Sans SignWriting Regular",
 )
#set math.equation(numbering: "(1)")
== Equation in 3 dimensions
Schrodinger’s equation:
$ i planck.reduce (diff Psi)/(diff t) = hat(H) Psi $ <seq>
Hamiltonian:
$ -(planck.reduce^2)/(2m) nabla^2 + V $
with time-independent potential $V$,
$ -(planck.reduce^2)/(2m)nabla^2 psi + V psi = E psi $
In spherical Coordinates
$ nabla^2 = 1/r^2(diff)/(diff r)(r^2(diff)/(diff r)) + 1/(r^2 sin theta)diff/(diff theta)(sin theta diff/(diff theta)) + 1/(r^2 sin^2 theta)(diff^2/(diff phi.alt^2)) $
along with separation of variables
$ psi(r, theta, phi.alt) = R(r)Y(theta, phi.alt) $
the @seq becomes
$ -planck.reduce^2/(2m)[Y/r^2 d/(d r)(r^2(d R)/(d r))+ R/(r^2 sin theta)diff/(diff theta)(sin theta (diff Y)/(diff theta))+R/(r^2 sin^2 theta)(diff^2 Y)/(diff phi.alt^2)] + V R Y = E R Y $
$ {1/R d/(d r)(r^2 (d R)/(d r)) - (2m r^2)/planck.reduce^2 [V(r) - E]} + 1/Y {1/(sin theta)diff/(diff theta)(sin theta (diff Y)/(diff theta)) + 1/(sin^2 theta)(diff^2 Y)/(diff phi.alt^2)} = 0 $
write the constant as $l(l+1)$:
$ 1/R d/(d r)(r^2 (d R)/(d r)) - (2m r^2)/planck.reduce^2 [V(r) - E] = l(l+1) $ <req>
and
$ sin theta diff/(diff theta)(sin theta (diff Y)/(diff theta)) + (diff^2 Y)/(diff phi.alt^2) = -l(l+1)Y sin^2 theta $
let $Y(theta, phi.alt) = Theta(theta)Phi(phi.alt)$, we have
$ {1/Theta [sin theta d/(d theta)(sin theta (d Theta)/(d theta))] + l(l+1)sin^2 theta} + 1/Phi (d^2 Phi)/(d phi.alt^2) = 0 $
Again but the separation constant is now $m^2$
$ 1/Theta [sin theta d/(d theta)(sin theta (d Theta)/(d theta))] + l(l+1)sin^2 theta = m^2 $ <aeq1>
$ 1/Phi (d^2 Phi)/(d phi.alt^2) =  -m^2 $ <aeq0>
solving @aeq0 and @aeq1
$ Phi(phi.alt) = e^(i m phi.alt) $
alas, in spherical coordinate, we have $Phi(phi.alt + 2 pi) = Phi(phi.alt)$, $m$ must be an $italic("integer")$, i.e.
$ m = 0, plus.minus 1, plus.minus 2, dots.h $
solution to @aeq1 is 
$ Theta(theta) = A P^m_l (cos theta) $
$ P^m_l(x) eq.triple (-1)^m (1-x^2)^(m slash 2)(d/(d x))^m P_l(x) $ <alf>
$ P_l(x) eq.triple 1/(2^l l!)(d/(d x))^l(x^2 - 1)^l $
according to @alf, $m$ should be no larger than $l$, there are $(2l + 1)$ possible values for $m$:
$ m = 0, plus.minus 1, dots.h, plus.minus l $
and of course $l$ should be a natural number.

The normalization condition:
$ integral abs(psi)^2 r^2 sin theta d r d theta d phi.alt = integral abs(R)^2 r^2 d r integral abs(Y)^2 d Omega = 1 $
IDK how to normalize this one, but we can do it separately
$ integral_0^oo abs(R)^2 r^2 d r = 1 and integral_0^pi integral_0^(2pi) abs(Y)^2 sin theta d theta d phi.alt = 1 $
The normalized angular wave functions are *spherical harmonics*
$ Y_l^m (theta, phi.alt) = sqrt((2l+1)/(4pi) ((l-m)!)/((l+m)!)) e^(i m phi.alt)P_l^m (cos theta) $
As for the radial equation @req, let $u(r) eq.triple r R(r)$,
$ R &= u/r\
  (d R) / (d r) &= [r (d u)/(d r) - u]/r^2\
  (d/(d r))[r^2 (d R)/(d r)] &= r (d 2 u)/(d r^2) $
the radia equation is then
$ -planck.reduce^2/(2m) (d^2 u)/(d r^2) + [V + planck.reduce^2/(2m)(l(l+1))/r^2]u = E u $
and the effective potential
$ V_"eff" = V + planck.reduce^2/(2m)(l(l+1))/r^2 $

== Hydrogen Atom
Potential
$ V(r) = -e^2/(4pi epsilon.alt_0) 1/r $
radial equation
$ -planck.reduce^2/(2m_e)(d^2 u)/(d r^2) + [-e^2/(4 pi epsilon.alt_0) 1/r + planck.reduce^2/(2 m_e)(l(l+1))/r^2]u = E u $
we only consider the bounding state, i.e. $E < 0$

Let $k eq.triple sqrt(-2m_e E)/planck.reduce$, which is real, we have
$ 1/kappa^2 (d^2 u)/(d r^2) = [1 - (m_e e^2)/(2pi epsilon.alt planck.reduce^2 kappa) 1/(kappa r) + (l(l+1))/((kappa r)^2)]u $
and again we change the variable and constant
$ rho eq.triple kappa r "and" rho_0 eq.triple (m_e e^2)/(2pi epsilon.alt planck.reduce^2 kappa) $ <28>
so that
$ (d^2 u)/(d rho^2) = [1 - rho_0/rho + (l(l+1))/rho^2]u $
as $rho -> oo$, $(d^2 u)/(d rho^2) = u$, $u(rho) = A e^(-rho) + B e^(rho)$, but _u_ have to converge, so $B = 0$, and the same holds for $rho -> 0$, where $u(rho) = C rho^(l+1)$, so we assume that
$ u(rho) = rho^(l+1)e^(-rho)v(rho) $
so
$ (d^2 u)/(d rho^2) = rho^l e^(-rho){[-2l-2+rho+(l(l+1))/rho]v + 2(l + 1 - rho)(d v)/ (d rho) + rho (d^2 v)/(d rho^2)} $
then the radial equation reads
$ rho (d^2 v)/(d rho^2) + 2(l + 1 - rho)(d v)/(d rho) + [rho_0 - 2(l + 1)]v = 0 $
using power series solution, we have
$ v(rho) = sum_(j = 0)^oo c_j rho^j $
$ c_(j+1) = [(2(j + l + 1) - rho_0)/((j+1)(j+2l+2))]c_j $ <recur>
for $v(rho)$ to converge, @recur must terminates, i.e. $exists N "such that" c_N = 0$, then we have
$ 2(N+l) - rho_0 = 0 $
according @28, we have
$ E = -(planck.reduce^2 kappa^2)/(2m) = -(m_e e^4)/(8 pi^2 epsilon.alt_0^2 planck.reduce^2 rho_0^2) $
Defining
$ n eq.triple N + l $ <principle>
we have $rho_0 = 2n$, and the famous result, *Bohr formula*
$ E_n = -[m_e/(2planck.reduce^2)(e^2/(4pi epsilon.alt_0))^2]1/n^2 = E_1/n^2, n = 1,2,3,dots.h $
since $E$ is function of $kappa$, it's obvious that
$ kappa = (m_e e^2)/(4 pi epsilon.alt_0 planck.reduce^2)1/n = 1/(a n) $
$a$ is so-called *Bohr radius*

The spatial wave function is labeled by three quantum nubers($n, l, m$)
$ psi_(n l m)(r, theta, phi.alt) = R_(n l)(r)Y_l^m(theta, phi.alt) $
notice that $l < n$, since if $N = 0$, $v(r)$ will vanish, i.e.
$ l = 0, 1, 2,dots.h, n - 1 $
and for each $l$ there are $2l+1$ possible vlaues of m, so the total degeneracy of the energy level $E_n$ is
$ d(n) = sum_(l = 0)^(n - 1)2l+1 = n^2 $
again, to formalize our solution
$ v(rho) = L_(n - l - 1)^(2l+1)(2rho) $
where
$ L_q^p(x) eq.triple (-1)^p (d/(d x))^p L_(p+q)(x) $
is an *associated Laguerre polynomial*, and the *Laguerre polynomial* is
$ L_q(x) eq.triple e^x/(q!)(d/(d x))^q (e^(-x) x^q) $
finally we give the normalized hydrogen wave function here
$ psi_(n l m) = sqrt((2/(n a))^3 ((n-l-1)!)/(2n(n+l)!))e^(-r slash n a) ((2r)/(n a))^l [L_(n-l-1)^(2l+1)(2r slash n a)]Y_l^m (theta, phi.alt) $
and the orthogonality
$ integral psi_(n l m)^* psi_(n' l' m') r^2 d r d Omega = delta_(n n') delta_(l l') delta_(m m') $

== Angular Momentum

$ bold(L) = bold(r) times bold(p) = vec(delim: "[", y p_z - z p_y, z p_x - x p_z, x p_y - y p_x) $
$ [L_i, L_j] = i planck.reduce epsilon.alt_(i j k) L_k $
$ sigma_(L_i)^2 sigma_(L_j)^2 gt.eq (1/(2i) angle.l i planck.reduce epsilon.alt_(i j k) L_k angle.r)^2 $
$ sigma_(L_i) sigma_(L_j) gt.eq planck.reduce/2 abs(angle.l epsilon.alt_(i j k) L_k angle.r) $
$ [L^2, bold(L)] = 0 $
we define a *ladder operator*
$ L_(plus.minus) eq.triple L_x plus.minus i L_y $
its commutator with $L_z$
$ [L_z, L_plus.minus] = [L_z, L_x] plus.minus i[L_z, L_y] = i planck.reduce L_y plus.minus i(-i planck.reduce L_x) = plus.minus planck.reduce L_plus.minus $
since $L^2$ and $L_plus.minus$ are commutative, we have
$ L^2 (L_plus.minus f) = L_plus.minus (L^2 f) = L_plus.minus (lambda f)  = lambda(L_plus.minus f) $
$L_plus.minus f$ is still an eigenfunction of $L^2$, with the same eigenvalue, but for $L_z$
$ L_z (L_plus.minus f) = (L_z L_plus.minus - L_plus.minus L_z) f + L_plus.minus L_z f = plus.minus planck.reduce L_plus.minus f + L_plus.minus (mu f) = (mu plus.minus planck.reduce)(L_plus.minus f) $
thus the operator changes the eigenvalue with an amount of $planck.reduce$ for $L_z$
since the _z_ component cannot exceed the total angular momentum, there must exists a top rung (maximal state) such that
$ L_+ f_t = 0 $
introduce some constants
$ L_z f_t = planck.reduce l f_t, space L^2 f_t = lambda f_t $
Then,
$ L_plus.minus L_minus.plus = (L_x plus.minus i L_y)(L_x minus.plus i L_y) = L_x^2 + L_y^2 minus.plus i[L_x, L_y] = L^2 - L_z^2 minus.plus i(i planck.reduce L_z) $
$ L^2 = L_plus.minus L_minus.plus + L_z^2 minus.plus planck.reduce L_z $
so we can say that
$ L^2 f_t = (L_- L_+ + L_z^2 + planck.reduce L_z)f_t = (0 + planck.reduce^2 l^2 + planck.reduce^2 l)f_t = planck.reduce^2 l(l+1)f_t $ <60>
and hence
$ lambda = planck.reduce^2 l(l+q) $
this indicates the eigenvalue of $L^2$ in terms of the _maximum_ eigenvalue of $L_z$

Of course there is the bottom rung where $L_- f_b = 0$, let's say
$ L_z f_b = planck.reduce macron(l)f_b, space L^2 f_b = lambda f_b $
using @60 again, we have $L^2 f_b = planck.reduce^2 macron(l)(macron(l) - 1)f_b$, and therefore
$ lambda = planck.reduce^2 macron(l)(macron(l)-1) $
