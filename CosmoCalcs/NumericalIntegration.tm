<TeXmacs|1.0.7.20>

<style|<tuple|generic|maxima>>

<\body>
  This file shows (though it does not prove) that using only even order
  derivatives on the edges to constrain the terms of an interpolating
  polynomial results in a series in h, the sample spacing, that is purely
  even.

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
      polyorder :2;
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      >>2>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>2) >
    <|unfolded-io>
      pint(x) := sum( a[i] * x^(i) / (i)!, i, 0, (polyorder + 1)-1 +
      2*ceiling((polyorder+1)/2) );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o2>)
      >><math-up|pint><around*|(|x|)>\<assign\><math-up|sum><around*|(|<frac|a<rsub|i>*x<rsup|i>|i!>,i,0,<math-up|polyorder>+1-1+2*<around*|\<lceil\>|<frac|<math-up|polyorder>+1|2>|\<rceil\>>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>3) >
    <|unfolded-io>
      pm(x) := sum( b[i] * x^(i) / (i)!, i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o3>)
      >><math-up|pm><around*|(|x|)>\<assign\><math-up|sum><around*|(|<frac|b<rsub|i>*x<rsup|i>|i!>,i,0,<math-up|polyorder>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>4) >
    <|unfolded-io>
      p0(x) := sum( c[i] * x^i / i!, i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o4>)
      >><with|math-font-family|rm|p0><around*|(|x|)>\<assign\><math-up|sum><around*|(|<frac|c<rsub|i>*x<rsup|i>|i!>,i,0,<math-up|polyorder>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>5) >
    <|unfolded-io>
      pp(x) := sum(d[i] * x^(i) / (i)!, i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o5>)
      >><math-up|pp><around*|(|x|)>\<assign\><math-up|sum><around*|(|<frac|d<rsub|i>*x<rsup|i>|i!>,i,0,<math-up|polyorder>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>6) >
    <|unfolded-io>
      ders(f) := makelist( diff(f(x), x, i), i, 0, polyorder );
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o6>)
      >><math-up|ders><around*|(|f|)>\<assign\><math-up|makelist><around*|(|<math-up|diff><around*|(|f<around*|(|x|)>,x,i|)>,i,0,<math-up|polyorder>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>7) >
    <|unfolded-io>
      derseven(f) := makelist( diff(f(x), x, 2*i), i, 0,
      ceiling((polyorder-1)/2));
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o7>)
      >><math-up|derseven><around*|(|f|)>\<assign\><math-up|makelist><around*|(|<math-up|diff><around*|(|f<around*|(|x|)>,x,2*i|)>,i,0,<around*|\<lceil\>|<frac|<math-up|polyorder>-1|2>|\<rceil\>>|)>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>8) >
    <|unfolded-io>
      dpint : ders(pint);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o8>)
      >><around*|[|<frac|a<rsub|6>*x<rsup|6>|720>+<frac|a<rsub|5>*x<rsup|5>|120>+<frac|a<rsub|4>*x<rsup|4>|24>+<frac|a<rsub|3>*x<rsup|3>|6>+<frac|a<rsub|2>*x<rsup|2>|2>+a<rsub|1>*x+a<rsub|0>,<frac|a<rsub|6>*x<rsup|5>|120>+<frac|a<rsub|5>*x<rsup|4>|24>+<frac|a<rsub|4>*x<rsup|3>|6>+<frac|a<rsub|3>*x<rsup|2>|2>+a<rsub|2>*x+a<rsub|1>,<frac|a<rsub|6>*x<rsup|4>|24>+<frac|a<rsub|5>*x<rsup|3>|6>+<frac|a<rsub|4>*x<rsup|2>|2>+a<rsub|3>*x+a<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>9) >
    <|unfolded-io>
      dpinteven : derseven(pint);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o9>)
      >><around*|[|<frac|a<rsub|6>*x<rsup|6>|720>+<frac|a<rsub|5>*x<rsup|5>|120>+<frac|a<rsub|4>*x<rsup|4>|24>+<frac|a<rsub|3>*x<rsup|3>|6>+<frac|a<rsub|2>*x<rsup|2>|2>+a<rsub|1>*x+a<rsub|0>,<frac|a<rsub|6>*x<rsup|4>|24>+<frac|a<rsub|5>*x<rsup|3>|6>+<frac|a<rsub|4>*x<rsup|2>|2>+a<rsub|3>*x+a<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>10) >
    <|unfolded-io>
      dpm : derseven(pm);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o10>)
      >><around*|[|<frac|b<rsub|2>*x<rsup|2>|2>+b<rsub|1>*x+b<rsub|0>,b<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>11) >
    <|unfolded-io>
      dp0 : ders(p0);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o11>)
      >><around*|[|<frac|c<rsub|2>*x<rsup|2>|2>+c<rsub|1>*x+c<rsub|0>,c<rsub|2>*x+c<rsub|1>,c<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>12) >
    <|unfolded-io>
      dpp : derseven(pp);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o12>)
      >><around*|[|<frac|d<rsub|2>*x<rsup|2>|2>+d<rsub|1>*x+d<rsub|0>,d<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>13) >
    <|unfolded-io>
      eqns : makelist( ev(dpinteven[i], x=-h) = ev(dpm[i], x=0), i, 1,
      length(dpm));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o13>)
      >><around*|[|<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>14) >
    <|unfolded-io>
      eqns : append( eqns, makelist( ev(dpint[i], x=0) = ev(dp0[i], x=0), i,
      1, length(dp0)));
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o14>)
      >><around*|[|<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>,a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>15) >
    <|unfolded-io>
      eqns : append( eqns, makelist( ev(dpinteven[i], x=h) = ev(dpp[i], x=0),
      i, 1, length(dpp)));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o15>)
      >><around*|[|<frac|a<rsub|6>*h<rsup|6>|720>-<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>-<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>-a<rsub|1>*h+a<rsub|0>=b<rsub|0>,<frac|a<rsub|6>*h<rsup|4>|24>-<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>-a<rsub|3>*h+a<rsub|2>=b<rsub|2>,a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>,<frac|a<rsub|6>*h<rsup|6>|720>+<frac|a<rsub|5>*h<rsup|5>|120>+<frac|a<rsub|4>*h<rsup|4>|24>+<frac|a<rsub|3>*h<rsup|3>|6>+<frac|a<rsub|2>*h<rsup|2>|2>+a<rsub|1>*h+a<rsub|0>=d<rsub|0>,<frac|a<rsub|6>*h<rsup|4>|24>+<frac|a<rsub|5>*h<rsup|3>|6>+<frac|a<rsub|4>*h<rsup|2>|2>+a<rsub|3>*h+a<rsub|2>=d<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>16) >
    <|unfolded-io>
      vars : makelist( a[i], i, 0, (polyorder+1)-1 +
      2*ceiling((polyorder+1)/2));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o16>)
      >><around*|[|a<rsub|0>,a<rsub|1>,a<rsub|2>,a<rsub|3>,a<rsub|4>,a<rsub|5>,a<rsub|6>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>17) >
    <|unfolded-io>
      soln : solve(eqns, vars);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o17>)
      >><around*|[|<around*|[|a<rsub|0>=c<rsub|0>,a<rsub|1>=c<rsub|1>,a<rsub|2>=c<rsub|2>,a<rsub|3>=<frac|<around*|(|3*b<rsub|2>-3*d<rsub|2>|)>*h<rsup|2>-120*c<rsub|1>*h+60*d<rsub|0>-60*b<rsub|0>|14*h<rsup|3>>,a<rsub|4>=-<frac|<around*|(|2*d<rsub|2>+56*c<rsub|2>+2*b<rsub|2>|)>*h<rsup|2>-60*d<rsub|0>+120*c<rsub|0>-60*b<rsub|0>|3*h<rsup|4>>,a<rsub|5>=-<frac|<around*|(|30*b<rsub|2>-30*d<rsub|2>|)>*h<rsup|2>-360*c<rsub|1>*h+180*d<rsub|0>-180*b<rsub|0>|7*h<rsup|5>>,a<rsub|6>=<frac|<around*|(|20*d<rsub|2>+200*c<rsub|2>+20*b<rsub|2>|)>*h<rsup|2>-240*d<rsub|0>+480*c<rsub|0>-240*b<rsub|0>|h<rsup|6>>|]>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>18) >
    <|unfolded-io>
      avec : transpose(matrix( makelist( rhs(soln[1][i]), i, 1, (polyorder+1)
      + 2*ceiling((polyorder+1)/2))));
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o18>)
      >><matrix|<tformat|<table|<row|<cell|c<rsub|0>>>|<row|<cell|c<rsub|1>>>|<row|<cell|c<rsub|2>>>|<row|<cell|<frac|<around*|(|3*b<rsub|2>-3*d<rsub|2>|)>*h<rsup|2>-120*c<rsub|1>*h+60*d<rsub|0>-60*b<rsub|0>|14*h<rsup|3>>>>|<row|<cell|-<frac|<around*|(|2*d<rsub|2>+56*c<rsub|2>+2*b<rsub|2>|)>*h<rsup|2>-60*d<rsub|0>+120*c<rsub|0>-60*b<rsub|0>|3*h<rsup|4>>>>|<row|<cell|-<frac|<around*|(|30*b<rsub|2>-30*d<rsub|2>|)>*h<rsup|2>-360*c<rsub|1>*h+180*d<rsub|0>-180*b<rsub|0>|7*h<rsup|5>>>>|<row|<cell|<frac|<around*|(|20*d<rsub|2>+200*c<rsub|2>+20*b<rsub|2>|)>*h<rsup|2>-240*d<rsub|0>+480*c<rsub|0>-240*b<rsub|0>|h<rsup|6>>>>>>>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>19) >
    <|unfolded-io>
      xvec : transpose(matrix(makelist( x^n/n!, n, 0, (polyorder+1)-1 +
      2*ceiling((polyorder+1)/2))));
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o19>)
      >><matrix|<tformat|<table|<row|<cell|1>>|<row|<cell|x>>|<row|<cell|<frac|x<rsup|2>|2>>>|<row|<cell|<frac|x<rsup|3>|6>>>|<row|<cell|<frac|x<rsup|4>|24>>>|<row|<cell|<frac|x<rsup|5>|120>>>|<row|<cell|<frac|x<rsup|6>|720>>>>>>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>20) >
    <|unfolded-io>
      intpoly : transpose(xvec).avec;
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o20>)
      >><frac|<around*|(|<around*|(|20*d<rsub|2>+200*c<rsub|2>+20*b<rsub|2>|)>*h<rsup|2>-240*d<rsub|0>+480*c<rsub|0>-240*b<rsub|0>|)>*x<rsup|6>|720*h<rsup|6>>-<frac|<around*|(|<around*|(|30*b<rsub|2>-30*d<rsub|2>|)>*h<rsup|2>-360*c<rsub|1>*h+180*d<rsub|0>-180*b<rsub|0>|)>*x<rsup|5>|840*h<rsup|5>>-<frac|<around*|(|<around*|(|2*d<rsub|2>+56*c<rsub|2>+2*b<rsub|2>|)>*h<rsup|2>-60*d<rsub|0>+120*c<rsub|0>-60*b<rsub|0>|)>*x<rsup|4>|72*h<rsup|4>>+<frac|<around*|(|<around*|(|3*b<rsub|2>-3*d<rsub|2>|)>*h<rsup|2>-120*c<rsub|1>*h+60*d<rsub|0>-60*b<rsub|0>|)>*x<rsup|3>|84*h<rsup|3>>+<frac|c<rsub|2>*x<rsup|2>|2>+c<rsub|1>*x+c<rsub|0>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>21) >
    <|unfolded-io>
      meanPoly : integrate(intpoly, x, -h, h)/(2*h);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o21>)
      >><frac|<frac|<around*|(|7*d<rsub|2>+256*c<rsub|2>-23*b<rsub|2>|)>*h<rsup|3>-1080*c<rsub|1>*h<rsup|2>+<around*|(|-120*d<rsub|0>+3840*c<rsub|0>+1320*b<rsub|0>|)>*h|5040>-<frac|<around*|(|23*d<rsub|2>-256*c<rsub|2>-7*b<rsub|2>|)>*h<rsup|3>-1080*c<rsub|1>*h<rsup|2>+<around*|(|-1320*d<rsub|0>-3840*c<rsub|0>+120*b<rsub|0>|)>*h|5040>|2*h>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>22) >
    <|unfolded-io>
      meanPoly : expand(meanPoly);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o22>)
      >>-<frac|d<rsub|2>*h<rsup|2>|630>+<frac|16*c<rsub|2>*h<rsup|2>|315>-<frac|b<rsub|2>*h<rsup|2>|630>+<frac|5*d<rsub|0>|42>+<frac|16*c<rsub|0>|21>+<frac|5*b<rsub|0>|42>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>23) >
    <|unfolded-io>
      meanPoly : fullratsimp(meanPoly);
    <|unfolded-io>
      \;

      \ <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o23>)
      >>-<frac|<around*|(|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|)>*h<rsup|2>-75*d<rsub|0>-480*c<rsub|0>-75*b<rsub|0>|630>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>59) >
    <|unfolded-io>
      meanCoeffs : makelist(diff(meanPoly, h, i)/i!, i, 0, polyorder);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o59>)
      >><around*|[|-<frac|<around*|(|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|)>*h<rsup|2>-75*d<rsub|0>-480*c<rsub|0>-75*b<rsub|0>|630>,-<frac|<around*|(|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|)>*h|315>,-<frac|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|630>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>60) >
    <|unfolded-io>
      meanCoeffs : ev(meanCoeffs, h = 0);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o60>)
      >><around*|[|-<frac|-75*d<rsub|0>-480*c<rsub|0>-75*b<rsub|0>|630>,0,-<frac|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|630>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>61) >
    <|unfolded-io>
      meanCoeffs : factor(meanCoeffs);
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o61>)
      >><around*|[|<frac|5*d<rsub|0>+32*c<rsub|0>+5*b<rsub|0>|42>,0,-<frac|d<rsub|2>-32*c<rsub|2>+b<rsub|2>|630>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>62) >
    <|unfolded-io>
      meanCoeffs, expand, numer;
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o62>)
      >><around*|[|0.119047619047619*d<rsub|0>+0.7619047619047619*c<rsub|0>+0.119047619047619*b<rsub|0>,0,-0.001587301587301587*d<rsub|2>+0.05079365079365079*c<rsub|2>-0.001587301587301587*b<rsub|2>|]>>>
    </unfolded-io>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>28) >
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