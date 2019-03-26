<TeXmacs|1.0.7.20>

<style|<tuple|generic|maxima>>

<\body>
  What follows is a derivation of the interpolation basis functions for
  interpolating a function known at 3 equally spaced points, separated by
  unit distances. The function and the first 3 derivatives of it must be
  known at these points to be useful. The c-basis functions of x multiply the
  appropriate derivatives of the function near the center point, the the
  d-basis functions multiply the derivatives to the right of the center, and
  the b-basis is identical to the d-basis evaluated at -x (by symmetry).

  <\session|maxima|default>
    <\output>
      \;

      Maxima 5.33.0 http://maxima.sourceforge.net

      using Lisp SBCL 1.1.18

      Distributed under the GNU Public License. See the file COPYING.

      Dedicated to the memory of William Schelter.

      The function bug_report() provides bug reporting information.

      <math|<with|math-display|true|<text|*** My very own personal
      maxima-init.mac has been loaded. *** >>>

      <\errput>
        STYLE-WARNING: redefining MAXIMA::MAIN-PROMPT in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-STRIPDOLLAR in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-MEXPT in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-CHOOSE in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-INT in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-SUM in DEFUN

        STYLE-WARNING: redefining MAXIMA::TEX-LSUM in DEFUN
      </errput>
    </output>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>1) >
    <|unfolded-io>
      polyorder : 3;
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      >>3>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>2) >
    <|unfolded-io>
      pint(x) := sum( a[i] * x^i / i!, i, 0, 3*(polyorder + 1)-1 );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o2>)
      >><math-up|pint><around*|(|x|)>\<assign\><math-up|sum><around*|(|<frac|a<rsub|i>*x<rsup|i>|i!>,i,0,3*<around*|(|<math-up|polyorder>+1|)>-1|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>3) >
    <|unfolded-io>
      ders(f) := makelist( diff(f(x), x, i), i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o3>)
      >><math-up|ders><around*|(|f|)>\<assign\><math-up|makelist><around*|(|<math-up|diff><around*|(|f<around*|(|x|)>,x,i|)>,i,0,<math-up|polyorder>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>4) >
    <|unfolded-io>
      dpint : ders(pint);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o4>)
      >><around*|[|<frac|a<rsub|11>*x<rsup|11>|39916800>+<frac|a<rsub|10>*x<rsup|10>|3628800>+<frac|a<rsub|9>*x<rsup|9>|362880>+<frac|a<rsub|8>*x<rsup|8>|40320>+<frac|a<rsub|7>*x<rsup|7>|5040>+<frac|a<rsub|6>*x<rsup|6>|720>+<frac|a<rsub|5>*x<rsup|5>|120>+<frac|a<rsub|4>*x<rsup|4>|24>+<frac|a<rsub|3>*x<rsup|3>|6>+<frac|a<rsub|2>*x<rsup|2>|2>+a<rsub|1>*x+a<rsub|0>,<frac|a<rsub|11>*x<rsup|10>|3628800>+<frac|a<rsub|10>*x<rsup|9>|362880>+<frac|a<rsub|9>*x<rsup|8>|40320>+<frac|a<rsub|8>*x<rsup|7>|5040>+<frac|a<rsub|7>*x<rsup|6>|720>+<frac|a<rsub|6>*x<rsup|5>|120>+<frac|a<rsub|5>*x<rsup|4>|24>+<frac|a<rsub|4>*x<rsup|3>|6>+<frac|a<rsub|3>*x<rsup|2>|2>+a<rsub|2>*x+a<rsub|1>,<frac|a<rsub|11>*x<rsup|9>|362880>+<frac|a<rsub|10>*x<rsup|8>|40320>+<frac|a<rsub|9>*x<rsup|7>|5040>+<frac|a<rsub|8>*x<rsup|6>|720>+<frac|a<rsub|7>*x<rsup|5>|120>+<frac|a<rsub|6>*x<rsup|4>|24>+<frac|a<rsub|5>*x<rsup|3>|6>+<frac|a<rsub|4>*x<rsup|2>|2>+a<rsub|3>*x+a<rsub|2>,<frac|a<rsub|11>*x<rsup|8>|40320>+<frac|a<rsub|10>*x<rsup|7>|5040>+<frac|a<rsub|9>*x<rsup|6>|720>+<frac|a<rsub|8>*x<rsup|5>|120>+<frac|a<rsub|7>*x<rsup|4>|24>+<frac|a<rsub|6>*x<rsup|3>|6>+<frac|a<rsub|5>*x<rsup|2>|2>+a<rsub|4>*x+a<rsub|3>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>5) >
    <|unfolded-io>
      eqns : makelist( ev(dpint[i], x=-h) = b[i-1], i, 1, polyorder+1);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o5>)
      >><around*|[|-<frac|a<rsub|11>*h<rsup|11>|39916800>+<frac|a<rsub|10>*h<rsup|10>|3628800>-<frac|a<rsub|9>*h<rsup|9>|362880>+<frac|a<rsub|8>*h<rsup|8>|40320>-<frac|a<rsub|7>*h<rsup|7>|5040>+<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|11>*h<rsup|10>|3628800>-<frac|a<rsub|10>*h<rsup|9>|362880>+<frac|a<rsub|9>*h<rsup|8>|40320>-<frac|a<rsub|8>*h<rsup|7>|5040>+<frac|a<rsub|7>*h<rsup|6>|720>-<frac|a<rsub|6>*h<rsup|5>|120>+<frac|a<rsub|5>*h<rsup|4>|24>-<frac|a<rsub|4>*h<rsup|3>|6>+<frac|a<rsub|3>*h<rsup|2>|2>-a<rsub|2>*h+a<rsub|1>=b<rsub|1>,-<frac|a<rsub|11>*h<rsup|9>|362880>+<frac|a<rsub|10>*h<rsup|8>|40320>-<frac|a<rsub|9>*h<rsup|7>|5040>+<frac|a<rsub|8>*h<rsup|6>|720>-<frac|a<rsub|7>*h<rsup|5>|120>+<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>,<frac|a<rsub|11>*h<rsup|8>|40320>-<frac|a<rsub|10>*h<rsup|7>|5040>+<frac|a<rsub|9>*h<rsup|6>|720>-<frac|a<rsub|8>*h<rsup|5>|120>+<frac|a<rsub|7>*h<rsup|4>|24>-<frac|a<rsub|6>*h<rsup|3>|6>+<frac|a<rsub|5>*h<rsup|2>|2>-a<rsub|4>*h+a<rsub|3>=b<rsub|3>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>6) >
    <|unfolded-io>
      eqns : append( eqns, makelist( ev(dpint[i], x=0) = c[i-1], i, 1,
      polyorder+1));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o6>)
      >><around*|[|-<frac|a<rsub|11>*h<rsup|11>|39916800>+<frac|a<rsub|10>*h<rsup|10>|3628800>-<frac|a<rsub|9>*h<rsup|9>|362880>+<frac|a<rsub|8>*h<rsup|8>|40320>-<frac|a<rsub|7>*h<rsup|7>|5040>+<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|11>*h<rsup|10>|3628800>-<frac|a<rsub|10>*h<rsup|9>|362880>+<frac|a<rsub|9>*h<rsup|8>|40320>-<frac|a<rsub|8>*h<rsup|7>|5040>+<frac|a<rsub|7>*h<rsup|6>|720>-<frac|a<rsub|6>*h<rsup|5>|120>+<frac|a<rsub|5>*h<rsup|4>|24>-<frac|a<rsub|4>*h<rsup|3>|6>+<frac|a<rsub|3>*h<rsup|2>|2>-a<rsub|2>*h+a<rsub|1>=b<rsub|1>,-<frac|a<rsub|11>*h<rsup|9>|362880>+<frac|a<rsub|10>*h<rsup|8>|40320>-<frac|a<rsub|9>*h<rsup|7>|5040>+<frac|a<rsub|8>*h<rsup|6>|720>-<frac|a<rsub|7>*h<rsup|5>|120>+<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>,<frac|a<rsub|11>*h<rsup|8>|40320>-<frac|a<rsub|10>*h<rsup|7>|5040>+<frac|a<rsub|9>*h<rsup|6>|720>-<frac|a<rsub|8>*h<rsup|5>|120>+<frac|a<rsub|7>*h<rsup|4>|24>-<frac|a<rsub|6>*h<rsup|3>|6>+<frac|a<rsub|5>*h<rsup|2>|2>-a<rsub|4>*h+a<rsub|3>=b<rsub|3>,a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>,a<rsub|3>=c<rsub|3>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>7) >
    <|unfolded-io>
      eqns : append( eqns, makelist( ev(dpint[i], x=h) = d[i-1], i, 1,
      polyorder+1));
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o7>)
      >><around*|[|-<frac|a<rsub|11>*h<rsup|11>|39916800>+<frac|a<rsub|10>*h<rsup|10>|3628800>-<frac|a<rsub|9>*h<rsup|9>|362880>+<frac|a<rsub|8>*h<rsup|8>|40320>-<frac|a<rsub|7>*h<rsup|7>|5040>+<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|11>*h<rsup|10>|3628800>-<frac|a<rsub|10>*h<rsup|9>|362880>+<frac|a<rsub|9>*h<rsup|8>|40320>-<frac|a<rsub|8>*h<rsup|7>|5040>+<frac|a<rsub|7>*h<rsup|6>|720>-<frac|a<rsub|6>*h<rsup|5>|120>+<frac|a<rsub|5>*h<rsup|4>|24>-<frac|a<rsub|4>*h<rsup|3>|6>+<frac|a<rsub|3>*h<rsup|2>|2>-a<rsub|2>*h+a<rsub|1>=b<rsub|1>,-<frac|a<rsub|11>*h<rsup|9>|362880>+<frac|a<rsub|10>*h<rsup|8>|40320>-<frac|a<rsub|9>*h<rsup|7>|5040>+<frac|a<rsub|8>*h<rsup|6>|720>-<frac|a<rsub|7>*h<rsup|5>|120>+<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>,<frac|a<rsub|11>*h<rsup|8>|40320>-<frac|a<rsub|10>*h<rsup|7>|5040>+<frac|a<rsub|9>*h<rsup|6>|720>-<frac|a<rsub|8>*h<rsup|5>|120>+<frac|a<rsub|7>*h<rsup|4>|24>-<frac|a<rsub|6>*h<rsup|3>|6>+<frac|a<rsub|5>*h<rsup|2>|2>-a<rsub|4>*h+a<rsub|3>=b<rsub|3>,a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>,a<rsub|3>=c<rsub|3>,<frac|a<rsub|11>*h<rsup|11>|39916800>+<frac|a<rsub|10>*h<rsup|10>|3628800>+<frac|a<rsub|9>*h<rsup|9>|362880>+<frac|a<rsub|8>*h<rsup|8>|40320>+<frac|a<rsub|7>*h<rsup|7>|5040>+<frac|a<rsub|6>*h<rsup|6>|720>+<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>+<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>+a<rsub|1>*h+a<rsub|0>=d<rsub|0>,<frac|a<rsub|11>*h<rsup|10>|3628800>+<frac|a<rsub|10>*h<rsup|9>|362880>+<frac|a<rsub|9>*h<rsup|8>|40320>+<frac|a<rsub|8>*h<rsup|7>|5040>+<frac|a<rsub|7>*h<rsup|6>|720>+<frac|a<rsub|6>*h<rsup|5>|120>+<frac|a<rsub|5>*h<rsup|4>|24>+<frac|a<rsub|4>*h<rsup|3>|6>+<frac|a<rsub|3>*h<rsup|2>|2>+a<rsub|2>*h+a<rsub|1>=d<rsub|1>,<frac|a<rsub|11>*h<rsup|9>|362880>+<frac|a<rsub|10>*h<rsup|8>|40320>+<frac|a<rsub|9>*h<rsup|7>|5040>+<frac|a<rsub|8>*h<rsup|6>|720>+<frac|a<rsub|7>*h<rsup|5>|120>+<frac|a<rsub|6>*h<rsup|4>|24>+<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>+a<rsub|3>*h+a<rsub|2>=d<rsub|2>,<frac|a<rsub|11>*h<rsup|8>|40320>+<frac|a<rsub|10>*h<rsup|7>|5040>+<frac|a<rsub|9>*h<rsup|6>|720>+<frac|a<rsub|8>*h<rsup|5>|120>+<frac|a<rsub|7>*h<rsup|4>|24>+<frac|a<rsub|6>*h<rsup|3>|6>+<frac|a<rsub|5>*h<rsup|2>|2>+a<rsub|4>*h+a<rsub|3>=d<rsub|3>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>8) >
    <|unfolded-io>
      vars : makelist( a[i], i, 0, 3*(polyorder+1)-1 );
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o8>)
      >><around*|[|a<rsub|0>,a<rsub|1>,a<rsub|2>,a<rsub|3>,a<rsub|4>,a<rsub|5>,a<rsub|6>,a<rsub|7>,a<rsub|8>,a<rsub|9>,a<rsub|10>,a<rsub|11>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>9) >
    <|unfolded-io>
      soln : solve(eqns, vars);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o9>)
      >><around*|[|<around*|[|a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>,a<rsub|3>=c<rsub|3>,a<rsub|4>=<frac|<around*|(|b<rsub|3>-d<rsub|3>|)>*h<rsup|3>+<around*|(|21*d<rsub|2>-192*c<rsub|2>+21*b<rsub|2>|)>*h<rsup|2>+<around*|(|165*b<rsub|1>-165*d<rsub|1>|)>*h+480*d<rsub|0>-960*c<rsub|0>+480*b<rsub|0>|4*h<rsup|4>>,a<rsub|5>=-<frac|<around*|(|5*d<rsub|3>+320*c<rsub|3>+5*b<rsub|3>|)>*h<rsup|3>+<around*|(|120*b<rsub|2>-120*d<rsub|2>|)>*h<rsup|2>+<around*|(|1065*d<rsub|1>+4800*c<rsub|1>+1065*b<rsub|1>|)>*h-3465*d<rsub|0>+3465*b<rsub|0>|4*h<rsup|5>>,a<rsub|6>=-<frac|<around*|(|45*b<rsub|3>-45*d<rsub|3>|)>*h<rsup|3>+<around*|(|855*d<rsub|2>-4320*c<rsub|2>+855*b<rsub|2>|)>*h<rsup|2>+<around*|(|5895*b<rsub|1>-5895*d<rsub|1>|)>*h+14400*d<rsub|0>-28800*c<rsub|0>+14400*b<rsub|0>|2*h<rsup|6>>,a<rsub|7>=<frac|<around*|(|315*d<rsub|3>+10080*c<rsub|3>+315*b<rsub|3>|)>*h<rsup|3>+<around*|(|6930*b<rsub|2>-6930*d<rsub|2>|)>*h<rsup|2>+<around*|(|55125*d<rsub|1>+201600*c<rsub|1>+55125*b<rsub|1>|)>*h-155925*d<rsub|0>+155925*b<rsub|0>|2*h<rsup|7>>,a<rsub|8>=<frac|<around*|(|1260*b<rsub|3>-1260*d<rsub|3>|)>*h<rsup|3>+<around*|(|21420*d<rsub|2>-80640*c<rsub|2>+21420*b<rsub|2>|)>*h<rsup|2>+<around*|(|132300*b<rsub|1>-132300*d<rsub|1>|)>*h+302400*d<rsub|0>-604800*c<rsub|0>+302400*b<rsub|0>|h<rsup|8>>,a<rsub|9>=-<frac|<around*|(|11340*d<rsub|3>+241920*c<rsub|3>+11340*b<rsub|3>|)>*h<rsup|3>+<around*|(|226800*b<rsub|2>-226800*d<rsub|2>|)>*h<rsup|2>+<around*|(|1644300*d<rsub|1>+5443200*c<rsub|1>+1644300*b<rsub|1>|)>*h-4365900*d<rsub|0>+4365900*b<rsub|0>|h<rsup|9>>,a<rsub|10>=-<frac|<around*|(|37800*b<rsub|3>-37800*d<rsub|3>|)>*h<rsup|3>+<around*|(|567000*d<rsub|2>-1814400*c<rsub|2>+567000*b<rsub|2>|)>*h<rsup|2>+<around*|(|3288600*b<rsub|1>-3288600*d<rsub|1>|)>*h+7257600*d<rsub|0>-14515200*c<rsub|0>+7257600*b<rsub|0>|h<rsup|10>>,a<rsub|11>=<frac|<around*|(|415800*d<rsub|3>+6652800*c<rsub|3>+415800*b<rsub|3>|)>*h<rsup|3>+<around*|(|7484400*b<rsub|2>-7484400*d<rsub|2>|)>*h<rsup|2>+<around*|(|51143400*d<rsub|1>+159667200*c<rsub|1>+51143400*b<rsub|1>|)>*h-130977000*d<rsub|0>+130977000*b<rsub|0>|h<rsup|11>>|]>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>10) >
    <|unfolded-io>
      avec : transpose(matrix( makelist( rhs(soln[1][i]), i, 1,
      3*(polyorder+1))));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o10>)
      >><matrix|<tformat|<table|<row|<cell|c<rsub|0>>>|<row|<cell|c<rsub|1>>>|<row|<cell|c<rsub|2>>>|<row|<cell|c<rsub|3>>>|<row|<cell|<frac|<around*|(|b<rsub|3>-d<rsub|3>|)>*h<rsup|3>+<around*|(|21*d<rsub|2>-192*c<rsub|2>+21*b<rsub|2>|)>*h<rsup|2>+<around*|(|165*b<rsub|1>-165*d<rsub|1>|)>*h+480*d<rsub|0>-960*c<rsub|0>+480*b<rsub|0>|4*h<rsup|4>>>>|<row|<cell|-<frac|<around*|(|5*d<rsub|3>+320*c<rsub|3>+5*b<rsub|3>|)>*h<rsup|3>+<around*|(|120*b<rsub|2>-120*d<rsub|2>|)>*h<rsup|2>+<around*|(|1065*d<rsub|1>+4800*c<rsub|1>+1065*b<rsub|1>|)>*h-3465*d<rsub|0>+3465*b<rsub|0>|4*h<rsup|5>>>>|<row|<cell|-<frac|<around*|(|45*b<rsub|3>-45*d<rsub|3>|)>*h<rsup|3>+<around*|(|855*d<rsub|2>-4320*c<rsub|2>+855*b<rsub|2>|)>*h<rsup|2>+<around*|(|5895*b<rsub|1>-5895*d<rsub|1>|)>*h+14400*d<rsub|0>-28800*c<rsub|0>+14400*b<rsub|0>|2*h<rsup|6>>>>|<row|<cell|<frac|<around*|(|315*d<rsub|3>+10080*c<rsub|3>+315*b<rsub|3>|)>*h<rsup|3>+<around*|(|6930*b<rsub|2>-6930*d<rsub|2>|)>*h<rsup|2>+<around*|(|55125*d<rsub|1>+201600*c<rsub|1>+55125*b<rsub|1>|)>*h-155925*d<rsub|0>+155925*b<rsub|0>|2*h<rsup|7>>>>|<row|<cell|<frac|<around*|(|1260*b<rsub|3>-1260*d<rsub|3>|)>*h<rsup|3>+<around*|(|21420*d<rsub|2>-80640*c<rsub|2>+21420*b<rsub|2>|)>*h<rsup|2>+<around*|(|132300*b<rsub|1>-132300*d<rsub|1>|)>*h+302400*d<rsub|0>-604800*c<rsub|0>+302400*b<rsub|0>|h<rsup|8>>>>|<row|<cell|-<frac|<around*|(|11340*d<rsub|3>+241920*c<rsub|3>+11340*b<rsub|3>|)>*h<rsup|3>+<around*|(|226800*b<rsub|2>-226800*d<rsub|2>|)>*h<rsup|2>+<around*|(|1644300*d<rsub|1>+5443200*c<rsub|1>+1644300*b<rsub|1>|)>*h-4365900*d<rsub|0>+4365900*b<rsub|0>|h<rsup|9>>>>|<row|<cell|-<frac|<around*|(|37800*b<rsub|3>-37800*d<rsub|3>|)>*h<rsup|3>+<around*|(|567000*d<rsub|2>-1814400*c<rsub|2>+567000*b<rsub|2>|)>*h<rsup|2>+<around*|(|3288600*b<rsub|1>-3288600*d<rsub|1>|)>*h+7257600*d<rsub|0>-14515200*c<rsub|0>+7257600*b<rsub|0>|h<rsup|10>>>>|<row|<cell|<frac|<around*|(|415800*d<rsub|3>+6652800*c<rsub|3>+415800*b<rsub|3>|)>*h<rsup|3>+<around*|(|7484400*b<rsub|2>-7484400*d<rsub|2>|)>*h<rsup|2>+<around*|(|51143400*d<rsub|1>+159667200*c<rsub|1>+51143400*b<rsub|1>|)>*h-130977000*d<rsub|0>+130977000*b<rsub|0>|h<rsup|11>>>>>>>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>11) >
    <|unfolded-io>
      xvec : transpose(matrix(makelist( x^n/n!, n, 0, 3*(polyorder+1)-1)));
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o11>)
      >><matrix|<tformat|<table|<row|<cell|1>>|<row|<cell|x>>|<row|<cell|<frac|x<rsup|2>|2>>>|<row|<cell|<frac|x<rsup|3>|6>>>|<row|<cell|<frac|x<rsup|4>|24>>>|<row|<cell|<frac|x<rsup|5>|120>>>|<row|<cell|<frac|x<rsup|6>|720>>>|<row|<cell|<frac|x<rsup|7>|5040>>>|<row|<cell|<frac|x<rsup|8>|40320>>>|<row|<cell|<frac|x<rsup|9>|362880>>>|<row|<cell|<frac|x<rsup|10>|3628800>>>|<row|<cell|<frac|x<rsup|11>|39916800>>>>>>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>12) >
    <|unfolded-io>
      intpoly : transpose(xvec).avec;
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o12>)
      >><frac|<around*|(|<around*|(|415800*d<rsub|3>+6652800*c<rsub|3>+415800*b<rsub|3>|)>*h<rsup|3>+<around*|(|7484400*b<rsub|2>-7484400*d<rsub|2>|)>*h<rsup|2>+<around*|(|51143400*d<rsub|1>+159667200*c<rsub|1>+51143400*b<rsub|1>|)>*h-130977000*d<rsub|0>+130977000*b<rsub|0>|)>*x<rsup|11>|39916800*h<rsup|11>>-<frac|<around*|(|<around*|(|37800*b<rsub|3>-37800*d<rsub|3>|)>*h<rsup|3>+<around*|(|567000*d<rsub|2>-1814400*c<rsub|2>+567000*b<rsub|2>|)>*h<rsup|2>+<around*|(|3288600*b<rsub|1>-3288600*d<rsub|1>|)>*h+7257600*d<rsub|0>-14515200*c<rsub|0>+7257600*b<rsub|0>|)>*x<rsup|10>|3628800*h<rsup|10>>-<frac|<around*|(|<around*|(|11340*d<rsub|3>+241920*c<rsub|3>+11340*b<rsub|3>|)>*h<rsup|3>+<around*|(|226800*b<rsub|2>-226800*d<rsub|2>|)>*h<rsup|2>+<around*|(|1644300*d<rsub|1>+5443200*c<rsub|1>+1644300*b<rsub|1>|)>*h-4365900*d<rsub|0>+4365900*b<rsub|0>|)>*x<rsup|9>|362880*h<rsup|9>>+<frac|<around*|(|<around*|(|1260*b<rsub|3>-1260*d<rsub|3>|)>*h<rsup|3>+<around*|(|21420*d<rsub|2>-80640*c<rsub|2>+21420*b<rsub|2>|)>*h<rsup|2>+<around*|(|132300*b<rsub|1>-132300*d<rsub|1>|)>*h+302400*d<rsub|0>-604800*c<rsub|0>+302400*b<rsub|0>|)>*x<rsup|8>|40320*h<rsup|8>>+<frac|<around*|(|<around*|(|315*d<rsub|3>+10080*c<rsub|3>+315*b<rsub|3>|)>*h<rsup|3>+<around*|(|6930*b<rsub|2>-6930*d<rsub|2>|)>*h<rsup|2>+<around*|(|55125*d<rsub|1>+201600*c<rsub|1>+55125*b<rsub|1>|)>*h-155925*d<rsub|0>+155925*b<rsub|0>|)>*x<rsup|7>|10080*h<rsup|7>>-<frac|<around*|(|<around*|(|45*b<rsub|3>-45*d<rsub|3>|)>*h<rsup|3>+<around*|(|855*d<rsub|2>-4320*c<rsub|2>+855*b<rsub|2>|)>*h<rsup|2>+<around*|(|5895*b<rsub|1>-5895*d<rsub|1>|)>*h+14400*d<rsub|0>-28800*c<rsub|0>+14400*b<rsub|0>|)>*x<rsup|6>|1440*h<rsup|6>>-<frac|<around*|(|<around*|(|5*d<rsub|3>+320*c<rsub|3>+5*b<rsub|3>|)>*h<rsup|3>+<around*|(|120*b<rsub|2>-120*d<rsub|2>|)>*h<rsup|2>+<around*|(|1065*d<rsub|1>+4800*c<rsub|1>+1065*b<rsub|1>|)>*h-3465*d<rsub|0>+3465*b<rsub|0>|)>*x<rsup|5>|480*h<rsup|5>>+<frac|<around*|(|<around*|(|b<rsub|3>-d<rsub|3>|)>*h<rsup|3>+<around*|(|21*d<rsub|2>-192*c<rsub|2>+21*b<rsub|2>|)>*h<rsup|2>+<around*|(|165*b<rsub|1>-165*d<rsub|1>|)>*h+480*d<rsub|0>-960*c<rsub|0>+480*b<rsub|0>|)>*x<rsup|4>|96*h<rsup|4>>+<frac|c<rsub|3>*x<rsup|3>|6>+<frac|c<rsub|2>*x<rsup|2>|2>+c<rsub|1>*x+c<rsub|0>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>13) >
    <|unfolded-io>
      bbases : makelist( diff(intpoly, b[i]), i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o13>)
      >><around*|[|<frac|105*x<rsup|11>|32*h<rsup|11>>-<frac|2*x<rsup|10>|h<rsup|10>>-<frac|385*x<rsup|9>|32*h<rsup|9>>+<frac|15*x<rsup|8>|2*h<rsup|8>>+<frac|495*x<rsup|7>|32*h<rsup|7>>-<frac|10*x<rsup|6>|h<rsup|6>>-<frac|231*x<rsup|5>|32*h<rsup|5>>+<frac|5*x<rsup|4>|h<rsup|4>>,<frac|41*x<rsup|11>|32*h<rsup|10>>-<frac|29*x<rsup|10>|32*h<rsup|9>>-<frac|145*x<rsup|9>|32*h<rsup|8>>+<frac|105*x<rsup|8>|32*h<rsup|7>>+<frac|175*x<rsup|7>|32*h<rsup|6>>-<frac|131*x<rsup|6>|32*h<rsup|5>>-<frac|71*x<rsup|5>|32*h<rsup|4>>+<frac|55*x<rsup|4>|32*h<rsup|3>>,<frac|3*x<rsup|11>|16*h<rsup|9>>-<frac|5*x<rsup|10>|32*h<rsup|8>>-<frac|5*x<rsup|9>|8*h<rsup|7>>+<frac|17*x<rsup|8>|32*h<rsup|6>>+<frac|11*x<rsup|7>|16*h<rsup|5>>-<frac|19*x<rsup|6>|32*h<rsup|4>>-<frac|x<rsup|5>|4*h<rsup|3>>+<frac|7*x<rsup|4>|32*h<rsup|2>>,<frac|x<rsup|11>|96*h<rsup|8>>-<frac|x<rsup|10>|96*h<rsup|7>>-<frac|x<rsup|9>|32*h<rsup|6>>+<frac|x<rsup|8>|32*h<rsup|5>>+<frac|x<rsup|7>|32*h<rsup|4>>-<frac|x<rsup|6>|32*h<rsup|3>>-<frac|x<rsup|5>|96*h<rsup|2>>+<frac|x<rsup|4>|96*h>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>14) >
    <|unfolded-io>
      cbases : makelist( diff(intpoly, c[i]), i, 0, polyorder);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o14>)
      >><around*|[|<frac|4*x<rsup|10>|h<rsup|10>>-<frac|15*x<rsup|8>|h<rsup|8>>+<frac|20*x<rsup|6>|h<rsup|6>>-<frac|10*x<rsup|4>|h<rsup|4>>+1,<frac|4*x<rsup|11>|h<rsup|10>>-<frac|15*x<rsup|9>|h<rsup|8>>+<frac|20*x<rsup|7>|h<rsup|6>>-<frac|10*x<rsup|5>|h<rsup|4>>+x,<frac|x<rsup|10>|2*h<rsup|8>>-<frac|2*x<rsup|8>|h<rsup|6>>+<frac|3*x<rsup|6>|h<rsup|4>>-<frac|2*x<rsup|4>|h<rsup|2>>+<frac|x<rsup|2>|2>,<frac|x<rsup|11>|6*h<rsup|8>>-<frac|2*x<rsup|9>|3*h<rsup|6>>+<frac|x<rsup|7>|h<rsup|4>>-<frac|2*x<rsup|5>|3*h<rsup|2>>+<frac|x<rsup|3>|6>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>15) >
    <|unfolded-io>
      dbases : makelist( diff(intpoly, d[i]), i, 0, polyorder);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o15>)
      >><around*|[|-<frac|105*x<rsup|11>|32*h<rsup|11>>-<frac|2*x<rsup|10>|h<rsup|10>>+<frac|385*x<rsup|9>|32*h<rsup|9>>+<frac|15*x<rsup|8>|2*h<rsup|8>>-<frac|495*x<rsup|7>|32*h<rsup|7>>-<frac|10*x<rsup|6>|h<rsup|6>>+<frac|231*x<rsup|5>|32*h<rsup|5>>+<frac|5*x<rsup|4>|h<rsup|4>>,<frac|41*x<rsup|11>|32*h<rsup|10>>+<frac|29*x<rsup|10>|32*h<rsup|9>>-<frac|145*x<rsup|9>|32*h<rsup|8>>-<frac|105*x<rsup|8>|32*h<rsup|7>>+<frac|175*x<rsup|7>|32*h<rsup|6>>+<frac|131*x<rsup|6>|32*h<rsup|5>>-<frac|71*x<rsup|5>|32*h<rsup|4>>-<frac|55*x<rsup|4>|32*h<rsup|3>>,-<frac|3*x<rsup|11>|16*h<rsup|9>>-<frac|5*x<rsup|10>|32*h<rsup|8>>+<frac|5*x<rsup|9>|8*h<rsup|7>>+<frac|17*x<rsup|8>|32*h<rsup|6>>-<frac|11*x<rsup|7>|16*h<rsup|5>>-<frac|19*x<rsup|6>|32*h<rsup|4>>+<frac|x<rsup|5>|4*h<rsup|3>>+<frac|7*x<rsup|4>|32*h<rsup|2>>,<frac|x<rsup|11>|96*h<rsup|8>>+<frac|x<rsup|10>|96*h<rsup|7>>-<frac|x<rsup|9>|32*h<rsup|6>>-<frac|x<rsup|8>|32*h<rsup|5>>+<frac|x<rsup|7>|32*h<rsup|4>>+<frac|x<rsup|6>|32*h<rsup|3>>-<frac|x<rsup|5>|96*h<rsup|2>>-<frac|x<rsup|4>|96*h>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>21) >
    <|unfolded-io>
      bbases - ev(dbases, x = -x)*makelist((-1)^i, i, 0, length(dbases)-1);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o21>)
      >><around*|[|0,0,0,0|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>20) >
    <|unfolded-io>
      makelist((-1)^i, i, 1, length(dbases)-1);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o20>)
      >><around*|[|-1,1,-1|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>18) >
    <|unfolded-io>
      %, factor;
    <|unfolded-io>
      \;

      rat: replaced 0.3102264404296875 by 20331/65536 = 0.3102264404296875

      \;

      rat: replaced -0.0902252197265625 by -5913/65536 = -0.0902252197265625

      \;

      rat: replaced 0.0098876953125 by 81/8192 = 0.0098876953125

      \;

      rat: replaced - \ \ \ 4.119873046875e-4 by -27/65536 = -
      \ \ \ 4.119873046875e-4

      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o18>)
      >><around*|[|<frac|3<rsup|4>*251|2<rsup|16>>,-<frac|3<rsup|4>*73|2<rsup|16>>,<frac|3<rsup|4>|2<rsup|13>>,-<frac|3<rsup|3>|2<rsup|16>>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>19) >
    <|unfolded-io>
      meanPoly : integrate(intpoly, x, -h, h)/(2*h);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o19>)
      >><frac|-<frac|<around*|(|205*d<rsub|3>-2464*c<rsub|3>-51*b<rsub|3>|)>*h<rsup|4>+<around*|(|-5638*d<rsub|2>-16384*c<rsub|2>-1018*b<rsub|2>|)>*h<rsup|3>+<around*|(|65355*d<rsub|1>-147840*c<rsub|1>-7605*b<rsub|1>|)>*h<rsup|2>+<around*|(|-374475*d<rsub|0>-491520*c<rsub|0>-21045*b<rsub|0>|)>*h|887040>-<frac|<around*|(|51*d<rsub|3>+2464*c<rsub|3>-205*b<rsub|3>|)>*h<rsup|4>+<around*|(|-1018*d<rsub|2>-16384*c<rsub|2>-5638*b<rsub|2>|)>*h<rsup|3>+<around*|(|7605*d<rsub|1>+147840*c<rsub|1>-65355*b<rsub|1>|)>*h<rsup|2>+<around*|(|-21045*d<rsub|0>-491520*c<rsub|0>-374475*b<rsub|0>|)>*h|887040>|2*h>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>20) >
    <|unfolded-io>
      meanPoly: expand(meanPoly);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o20>)
      >>-<frac|d<rsub|3>*h<rsup|3>|6930>+<frac|b<rsub|3>*h<rsup|3>|6930>+<frac|13*d<rsub|2>*h<rsup|2>|3465>+<frac|64*c<rsub|2>*h<rsup|2>|3465>+<frac|13*b<rsub|2>*h<rsup|2>|3465>-<frac|19*d<rsub|1>*h|462>+<frac|19*b<rsub|1>*h|462>+<frac|103*d<rsub|0>|462>+<frac|128*c<rsub|0>|231>+<frac|103*b<rsub|0>|462>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>21) >
    <|unfolded-io>
      meanPoly;
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o21>)
      >>-<frac|d<rsub|3>*h<rsup|3>|6930>+<frac|b<rsub|3>*h<rsup|3>|6930>+<frac|13*d<rsub|2>*h<rsup|2>|3465>+<frac|64*c<rsub|2>*h<rsup|2>|3465>+<frac|13*b<rsub|2>*h<rsup|2>|3465>-<frac|19*d<rsub|1>*h|462>+<frac|19*b<rsub|1>*h|462>+<frac|103*d<rsub|0>|462>+<frac|128*c<rsub|0>|231>+<frac|103*b<rsub|0>|462>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>22) >
    <|unfolded-io>
      meanPoly : fullratsimp(meanPoly);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o22>)
      >>-<frac|<around*|(|d<rsub|3>-b<rsub|3>|)>*h<rsup|3>+<around*|(|-26*d<rsub|2>-128*c<rsub|2>-26*b<rsub|2>|)>*h<rsup|2>+<around*|(|285*d<rsub|1>-285*b<rsub|1>|)>*h-1545*d<rsub|0>-3840*c<rsub|0>-1545*b<rsub|0>|6930>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>23) >
    <|unfolded-io>
      meanCoeffs : makelist(diff(meanPoly, h, i)/i!, i, 0, polyorder);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o23>)
      >><around*|[|-<frac|<around*|(|d<rsub|3>-b<rsub|3>|)>*h<rsup|3>+<around*|(|-26*d<rsub|2>-128*c<rsub|2>-26*b<rsub|2>|)>*h<rsup|2>+<around*|(|285*d<rsub|1>-285*b<rsub|1>|)>*h-1545*d<rsub|0>-3840*c<rsub|0>-1545*b<rsub|0>|6930>,-<frac|3*<around*|(|d<rsub|3>-b<rsub|3>|)>*h<rsup|2>+2*<around*|(|-26*d<rsub|2>-128*c<rsub|2>-26*b<rsub|2>|)>*h+285*d<rsub|1>-285*b<rsub|1>|6930>,-<frac|6*<around*|(|d<rsub|3>-b<rsub|3>|)>*h+2*<around*|(|-26*d<rsub|2>-128*c<rsub|2>-26*b<rsub|2>|)>|13860>,-<frac|d<rsub|3>-b<rsub|3>|6930>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>24) >
    <|unfolded-io>
      meanCoeffs : ev(meanCoeffs, h = 0);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o24>)
      >><around*|[|-<frac|-1545*d<rsub|0>-3840*c<rsub|0>-1545*b<rsub|0>|6930>,-<frac|285*d<rsub|1>-285*b<rsub|1>|6930>,-<frac|-26*d<rsub|2>-128*c<rsub|2>-26*b<rsub|2>|6930>,-<frac|d<rsub|3>-b<rsub|3>|6930>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>25) >
    <|unfolded-io>
      meanCoeffs : factor(meanCoeffs);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o25>)
      >><around*|[|<frac|103*d<rsub|0>+256*c<rsub|0>+103*b<rsub|0>|462>,-<frac|19*<around*|(|d<rsub|1>-b<rsub|1>|)>|462>,<frac|13*d<rsub|2>+64*c<rsub|2>+13*b<rsub|2>|3465>,-<frac|d<rsub|3>-b<rsub|3>|6930>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>26) >
    <|unfolded-io>
      meanCoeffs, expand, numer;
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o26>)
      >><around*|[|0.2229437229437229*d<rsub|0>+0.5541125541125541*c<rsub|0>+0.2229437229437229*b<rsub|0>,0.04112554112554113*b<rsub|1>-0.04112554112554113*d<rsub|1>,0.003751803751803752*d<rsub|2>+0.01847041847041847*c<rsub|2>+0.003751803751803752*b<rsub|2>,1.4430014430014432\<times\>10<rsup|-4>*b<rsub|3>-1.4430014430014432\<times\>10<rsup|-4>*d<rsub|3>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>27) >
    <|unfolded-io>
      kill(all);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o0>)
      >><math-bf|done>>>
    </unfolded-io>

    <\input>
      <with|color|red|(<with|math-font-family|rm|%i>1) >
    <|input>
      \;
    </input>
  </session>
</body>