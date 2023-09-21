#set text(
 font: "Noto Sans SignWriting Regular",
)
#set math.equation(numbering: "(1)")
== Pade approximant
Pade approximant of order [$m slash n]$ is the rational function also denoted as $[m slash n]_f (x)$
$ R(x) = (sum_(j=0)^m a_j x^j)/(1+sum_(k=1)^n b_k x^k) = (a_0 + a_1 x+ a_2 x^2 + dots.h + a_m x^m)/(1 + b_1 x + b_2 x^2 + dots.h + b_n x^n) $
which agrees with $f(x)$ to the highest possible order, i.e.
$ f(0) &= R(0) \ f'(0) &= R'(0) \ f''(0) &= R''(0) \ &dots.v \ f^((m+n)) (0) &= R^((m+n))(0) $
To calculate the coefficients, we write the formula
$ [m slash n]_f (x) = f(0) + f'(0) x + (f''(0))/(k!) x^2 + dots.h +  (f^((m + n))(0))/((m+n)!) x^(m+n) $
this says the pade approximant is an approximation of function $f$ to the order of $m+n $, i.e.
$ f(x) = [m slash n]_f (x) + O(x^(m+n+1)) $
plugging in the definition and solving the undetermined coefficients
$ sum_(j=0)^m a_j x^j = (1 + sum_(k = 1)^n b_k x^k)sum_(i = 0)^(m+n) (f^((k))(0))/(k!)x^k $
