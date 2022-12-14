/*!glpk.js - v4.49.0
 * https://github.com/hgourvest/glpk.js
 * Copyright (c) 2013 Henri Gourvest; Licensed GPLv2*/
(function (exports) {
  var s = Number.MAX_VALUE,
    aa = Number.MIN_VALUE;
  function w(a) {
    throw Error(a);
  }
  function x() {}
  exports.glp_get_print_func = function () {
    return x;
  };
  exports.glp_set_print_func = function (a) {
    x = a;
  };
  function da(a, b) {
    for (var c in b) a[c] = b[c];
  }
  function ga(a, b, c, d, e) {
    for (; 0 < e; b++, d++, e--) a[b] = c[d];
  }
  function ha(a, b, c, d) {
    for (; 0 < d; b++, d--) a[b] = c;
  }
  function ia(a, b, c) {
    for (; 0 < c; b++, c--) a[b] = {};
  }
  function ja() {
    return new Date().getTime();
  }
  function la(a) {
    return (ja() - a) / 1e3;
  }
  function ma(a, b, c) {
    var d = Array(b);
    ga(d, 0, a, 1, b);
    d.sort(c);
    ga(a, 1, d, 0, b);
  }
  var na = {},
    ra = (exports.glp_version = function () {
      return pa + "." + qa;
    });
  function ta(a) {
    a = "string" == typeof a ? a.charCodeAt(0) : -1;
    return (0 <= a && 31 >= a) || 127 == a;
  }
  function ua(a) {
    a = "string" == typeof a ? a.charCodeAt(0) : -1;
    return (65 <= a && 90 >= a) || (97 <= a && 122 >= a);
  }
  function va(a) {
    a = "string" == typeof a ? a.charCodeAt(0) : -1;
    return (
      (65 <= a && 90 >= a) || (97 <= a && 122 >= a) || (48 <= a && 57 >= a)
    );
  }
  function wa(a) {
    a = "string" == typeof a ? a.charCodeAt(0) : -1;
    return 48 <= a && 57 >= a;
  }
  function xa() {
    function a(a, d, e, k, l, p, m) {
      a >>>= 0;
      e = (e && a && { 2: "0b", 8: "0", 16: "0x" }[d]) || "";
      a = e + c(a.toString(d), p || 0, "0", !1);
      return b(a, e, k, l, m);
    }
    function b(a, b, d, e, l, p) {
      var m = e - a.length;
      0 < m &&
        (a =
          d || !l
            ? c(a, e, p, d)
            : a.slice(0, b.length) + c("", m, "0", !0) + a.slice(b.length));
      return a;
    }
    function c(a, b, c, d) {
      c || (c = " ");
      b = a.length >= b ? "" : Array((1 + b - a.length) >>> 0).join(c);
      return d ? a + b : b + a;
    }
    var d = arguments,
      e = 0;
    return d[e++].replace(
      /%%|%(\d+\$)?([-+\'#0 ]*)(\*\d+\$|\*|\d+)?(\.(\*\d+\$|\*|\d+))?([scboxXuideEfFgG])/g,
      function (f, g, h, k, l, p, m) {
        var q, r;
        if ("%%" == f) return "%";
        var n = !1;
        r = "";
        var t = (l = !1);
        q = " ";
        for (var y = h.length, E = 0; h && E < y; E++)
          switch (h.charAt(E)) {
            case " ":
              r = " ";
              break;
            case "+":
              r = "+";
              break;
            case "-":
              n = !0;
              break;
            case "'":
              q = h.charAt(E + 1);
              break;
            case "0":
              l = !0;
              break;
            case "#":
              t = !0;
          }
        k = k
          ? "*" == k
            ? +d[e++]
            : "*" == k.charAt(0)
            ? +d[k.slice(1, -1)]
            : +k
          : 0;
        0 > k && ((k = -k), (n = !0));
        if (!isFinite(k))
          throw Error("sprintf: (minimum-)width must be finite");
        p = p
          ? "*" == p
            ? +d[e++]
            : "*" == p.charAt(0)
            ? +d[p.slice(1, -1)]
            : +p
          : -1 < "fFeE".indexOf(m)
          ? 6
          : "d" == m
          ? 0
          : void 0;
        g = g ? d[g.slice(0, -1)] : d[e++];
        switch (m) {
          case "s":
            return (
              (m = String(g)),
              null != p && (m = m.slice(0, p)),
              b(m, "", n, k, l, q)
            );
          case "c":
            return (
              (m = String.fromCharCode(+g)),
              null != p && (m = m.slice(0, p)),
              b(m, "", n, k, l, void 0)
            );
          case "b":
            return a(g, 2, t, n, k, p, l);
          case "o":
            return a(g, 8, t, n, k, p, l);
          case "x":
            return a(g, 16, t, n, k, p, l);
          case "X":
            return a(g, 16, t, n, k, p, l).toUpperCase();
          case "u":
            return a(g, 10, t, n, k, p, l);
          case "i":
          case "d":
            return (
              (q = +g || 0),
              (q = Math.round(q - (q % 1))),
              (f = 0 > q ? "-" : r),
              (g = f + c(String(Math.abs(q)), p, "0", !1)),
              b(g, f, n, k, l)
            );
          case "e":
          case "E":
          case "f":
          case "F":
          case "g":
          case "G":
            return (
              (q = +g),
              (f = 0 > q ? "-" : r),
              (r = ["toExponential", "toFixed", "toPrecision"][
                "efg".indexOf(m.toLowerCase())
              ]),
              (m = ["toString", "toUpperCase"]["eEfFgG".indexOf(m) % 2]),
              (g = f + Math.abs(q)[r](p)),
              b(g, f, n, k, l)[m]()
            );
          default:
            return f;
        }
      }
    );
  }
  function ya(a) {
    a.Dd = 3621377730;
    a.ke = null;
    a.V = null;
    a.name = null;
    a.eb = null;
    a.dir = za;
    a.ha = 0;
    a.hb = 100;
    a.K = 200;
    a.g = a.i = 0;
    a.L = 0;
    a.n = Array(1 + a.hb);
    a.f = Array(1 + a.K);
    a.fc = {};
    a.Lc = {};
    a.valid = 0;
    a.head = new Int32Array(1 + a.hb);
    a.Qd = null;
    a.U = null;
    a.na = a.sa = Aa;
    a.aa = 0;
    a.$ = 0;
    a.some = 0;
    a.df = Aa;
    a.ae = 0;
    a.za = Aa;
    a.ta = 0;
  }
  var Ba = (exports.glp_create_prob = function () {
      var a = {};
      ya(a);
      return a;
    }),
    Ca = (exports.glp_set_prob_name = function (a, b) {
      var c = a.V;
      null != c &&
        0 != c.reason &&
        w("glp_set_prob_name: operation not allowed");
      a.name = b;
    }),
    Da = (exports.glp_set_obj_name = function (a, b) {
      var c = a.V;
      null != c &&
        0 != c.reason &&
        w("glp_set_obj_name: operation not allowed");
      a.eb = b;
    }),
    Fa = (exports.glp_set_obj_dir = function (a, b) {
      var c = a.V;
      null != c && 0 != c.reason && w("glp_set_obj_dir: operation not allowed");
      b != za &&
        b != Ea &&
        w("glp_set_obj_dir: dir = " + b + "; invalid direction flag");
      a.dir = b;
    }),
    La = (exports.glp_add_rows = function (a, b) {
      var c = a.V,
        d;
      1 > b && w("glp_add_rows: nrs = " + b + "; invalid number of rows");
      b > 1e8 - a.g && w("glp_add_rows: nrs = " + b + "; too many rows");
      var e = a.g + b;
      if (a.hb < e) {
        for (; a.hb < e; ) a.hb += a.hb;
        a.n.length = 1 + a.hb;
        a.head = new Int32Array(1 + a.hb);
      }
      for (var f = a.g + 1; f <= e; f++) {
        a.n[f] = d = {};
        d.ea = f;
        d.name = null;
        d.rb = null;
        d.La = 0;
        d.origin = 0;
        d.qc = 0;
        if (null != c)
          switch (c.reason) {
            case Ga:
              d.La = c.N.La;
              d.origin = Ha;
              break;
            case Ia:
              (d.La = c.N.La), (d.origin = Ja);
          }
        d.type = Ka;
        d.c = d.d = 0;
        d.k = null;
        d.ma = 1;
        d.m = A;
        d.bind = 0;
        d.r = d.J = 0;
        d.dc = d.mc = 0;
        d.Sa = 0;
      }
      a.g = e;
      a.valid = 0;
      null != c && 0 != c.reason && (c.pe = 1);
      return e - b + 1;
    }),
    Oa = (exports.glp_add_cols = function (a, b) {
      var c = a.V;
      null != c && 0 != c.reason && w("glp_add_cols: operation not allowed");
      1 > b && w("glp_add_cols: ncs = " + b + "; invalid number of columns");
      b > 1e8 - a.i && w("glp_add_cols: ncs = " + b + "; too many columns");
      var d = a.i + b;
      if (a.K < d) {
        for (; a.K < d; ) a.K += a.K;
        a.f.length = 1 + a.K;
      }
      for (var e = a.i + 1; e <= d; e++)
        (a.f[e] = c = {}),
          (c.C = e),
          (c.name = null),
          (c.rb = null),
          (c.kind = Ma),
          (c.type = B),
          (c.c = c.d = 0),
          (c.u = 0),
          (c.k = null),
          (c.va = 1),
          (c.m = Na),
          (c.bind = 0),
          (c.r = c.J = 0),
          (c.dc = c.mc = 0),
          (c.Sa = 0);
      a.i = d;
      return d - b + 1;
    }),
    Pa = (exports.glp_set_row_name = function (a, b, c) {
      (1 <= b && b <= a.g) ||
        w("glp_set_row_name: i = " + b + "; row number out of range");
      b = a.n[b];
      null != b.name && (delete a.fc[b.name], (b.name = null));
      null != c && ((b.name = c), (a.fc[b.name] = b));
    }),
    Qa = (exports.glp_set_col_name = function (a, b, c) {
      var d = a.V;
      null != d &&
        0 != d.reason &&
        w("glp_set_col_name: operation not allowed");
      (1 <= b && b <= a.i) ||
        w("glp_set_col_name: j = " + b + "; column number out of range");
      b = a.f[b];
      null != b.name && (delete a.Lc[b.name], (b.name = null));
      null != c && ((b.name = c), (a.Lc[b.name] = b));
    }),
    Va = (exports.glp_set_row_bnds = function (a, b, c, d, e) {
      (1 <= b && b <= a.g) ||
        w("glp_set_row_bnds: i = " + b + "; row number out of range");
      a = a.n[b];
      a.type = c;
      switch (c) {
        case Ka:
          a.c = a.d = 0;
          a.m != A && (a.m = Ra);
          break;
        case Sa:
          a.c = d;
          a.d = 0;
          a.m != A && (a.m = G);
          break;
        case Ta:
          a.c = 0;
          a.d = e;
          a.m != A && (a.m = Ua);
          break;
        case I:
          a.c = d;
          a.d = e;
          a.m != A &&
            a.m != G &&
            a.m != Ua &&
            (a.m = Math.abs(d) <= Math.abs(e) ? G : Ua);
          break;
        case B:
          a.c = a.d = d;
          a.m != A && (a.m = Na);
          break;
        default:
          w(
            "glp_set_row_bnds: i = " +
              b +
              "; type = " +
              c +
              "; invalid row type"
          );
      }
    }),
    Wa = (exports.glp_set_col_bnds = function (a, b, c, d, e) {
      (1 <= b && b <= a.i) ||
        w("glp_set_col_bnds: j = " + b + "; column number out of range");
      a = a.f[b];
      a.type = c;
      switch (c) {
        case Ka:
          a.c = a.d = 0;
          a.m != A && (a.m = Ra);
          break;
        case Sa:
          a.c = d;
          a.d = 0;
          a.m != A && (a.m = G);
          break;
        case Ta:
          a.c = 0;
          a.d = e;
          a.m != A && (a.m = Ua);
          break;
        case I:
          a.c = d;
          a.d = e;
          a.m != A &&
            a.m != G &&
            a.m != Ua &&
            (a.m = Math.abs(d) <= Math.abs(e) ? G : Ua);
          break;
        case B:
          a.c = a.d = d;
          a.m != A && (a.m = Na);
          break;
        default:
          w(
            "glp_set_col_bnds: j = " +
              b +
              "; type = " +
              c +
              "; invalid column type"
          );
      }
    }),
    Xa = (exports.glp_set_obj_coef = function (a, b, c) {
      var d = a.V;
      null != d &&
        0 != d.reason &&
        w("glp_set_obj_coef: operation not allowed");
      (0 <= b && b <= a.i) ||
        w("glp_set_obj_coef: j = " + b + "; column number out of range");
      0 == b ? (a.ha = c) : (a.f[b].u = c);
    }),
    Ya = (exports.glp_set_mat_row = function (a, b, c, d, e) {
      var f, g, h;
      (1 <= b && b <= a.g) ||
        w("glp_set_mat_row: i = " + b + "; row number out of range");
      for (var k = a.n[b]; null != k.k; )
        (g = k.k),
          (k.k = g.B),
          (f = g.f),
          null == g.ra ? (f.k = g.I) : (g.ra.I = g.I),
          null != g.I && (g.I.ra = g.ra),
          a.L--,
          f.m == A && (a.valid = 0);
      (0 <= c && c <= a.i) ||
        w(
          "glp_set_mat_row: i = " + b + "; len = " + c + "; invalid row length "
        );
      c > 5e8 - a.L &&
        w(
          "glp_set_mat_row: i = " +
            b +
            "; len = " +
            c +
            "; too many constraint coefficients"
        );
      for (h = 1; h <= c; h++)
        (g = d[h]),
          (1 <= g && g <= a.i) ||
            w(
              "glp_set_mat_row: i = " +
                b +
                "; ind[" +
                h +
                "] = " +
                g +
                "; column index out of range"
            ),
          (f = a.f[g]),
          null != f.k &&
            f.k.n.ea == b &&
            w(
              "glp_set_mat_row: i = " +
                b +
                "; ind[" +
                h +
                "] = " +
                g +
                "; duplicate column indices not allowed"
            ),
          (g = {}),
          a.L++,
          (g.n = k),
          (g.f = f),
          (g.j = e[h]),
          (g.ua = null),
          (g.B = k.k),
          (g.ra = null),
          (g.I = f.k),
          null != g.B && (g.B.ua = g),
          null != g.I && (g.I.ra = g),
          (k.k = f.k = g),
          f.m == A && 0 != g.j && (a.valid = 0);
      for (g = k.k; null != g; g = b)
        (b = g.B),
          0 == g.j &&
            (null == g.ua ? (k.k = b) : (g.ua.B = b),
            null != b && (b.ua = g.ua),
            (g.f.k = g.I),
            null != g.I && (g.I.ra = null),
            a.L--);
    }),
    Za = (exports.glp_set_mat_col = function (a, b, c, d, e) {
      var f = a.V,
        g,
        h,
        k;
      null != f && 0 != f.reason && w("glp_set_mat_col: operation not allowed");
      (1 <= b && b <= a.i) ||
        w("glp_set_mat_col: j = " + b + "; column number out of range");
      for (f = a.f[b]; null != f.k; )
        (h = f.k),
          (f.k = h.I),
          (g = h.n),
          null == h.ua ? (g.k = h.B) : (h.ua.B = h.B),
          null != h.B && (h.B.ua = h.ua),
          a.L--;
      (0 <= c && c <= a.g) ||
        w(
          "glp_set_mat_col: j = " +
            b +
            "; len = " +
            c +
            "; invalid column length"
        );
      c > 5e8 - a.L &&
        w(
          "glp_set_mat_col: j = " +
            b +
            "; len = " +
            c +
            "; too many constraint coefficients"
        );
      for (k = 1; k <= c; k++)
        (h = d[k]),
          (1 <= h && h <= a.g) ||
            w(
              "glp_set_mat_col: j = " +
                b +
                "; ind[" +
                k +
                "] = " +
                h +
                "; row index out of range"
            ),
          (g = a.n[h]),
          null != g.k &&
            g.k.f.C == b &&
            w(
              "glp_set_mat_col: j = " +
                b +
                "; ind[" +
                k +
                "] = " +
                h +
                "; duplicate row indices not allowed"
            ),
          (h = {}),
          a.L++,
          (h.n = g),
          (h.f = f),
          (h.j = e[k]),
          (h.ua = null),
          (h.B = g.k),
          (h.ra = null),
          (h.I = f.k),
          null != h.B && (h.B.ua = h),
          null != h.I && (h.I.ra = h),
          (g.k = f.k = h);
      for (h = f.k; null != h; h = b)
        (b = h.I),
          0 == h.j &&
            ((h.n.k = h.B),
            null != h.B && (h.B.ua = null),
            null == h.ra ? (f.k = b) : (h.ra.I = b),
            null != b && (b.ra = h.ra),
            a.L--);
      f.m == A && (a.valid = 0);
    });
  exports.glp_load_matrix = function (a, b, c, d, e) {
    var f = a.V,
      g,
      h,
      k,
      l;
    null != f && 0 != f.reason && w("glp_load_matrix: operation not allowed");
    for (k = 1; k <= a.g; k++)
      for (f = a.n[k]; null != f.k; ) (h = f.k), (f.k = h.B), a.L--;
    for (h = 1; h <= a.i; h++) a.f[h].k = null;
    0 > b &&
      w(
        "glp_load_matrix: ne = " +
          b +
          "; invalid number of constraint coefficients"
      );
    5e8 < b &&
      w("glp_load_matrix: ne = " + b + "; too many constraint coefficients");
    for (l = 1; l <= b; l++)
      (k = c[l]),
        (h = d[l]),
        (1 <= k && k <= a.g) ||
          w(
            "glp_load_matrix: ia[" + l + "] = " + k + "; row index out of range"
          ),
        (f = a.n[k]),
        (1 <= h && h <= a.i) ||
          w(
            "glp_load_matrix: ja[" +
              l +
              "] = " +
              h +
              "; column index out of range"
          ),
        (g = a.f[h]),
        (h = {}),
        a.L++,
        (h.n = f),
        (h.f = g),
        (h.j = e[l]),
        (h.ua = null),
        (h.B = f.k),
        null != h.B && (h.B.ua = h),
        (f.k = h);
    for (k = 1; k <= a.g; k++)
      for (h = a.n[k].k; null != h; h = h.B) {
        g = h.f;
        if (null != g.k && g.k.n.ea == k) {
          for (l = 1; l <= b && (c[l] != k || d[l] != g.C); l++);
          w(
            "glp_load_mat: ia[" +
              l +
              "] = " +
              k +
              "; ja[" +
              l +
              "] = " +
              g.C +
              "; duplicate indices not allowed"
          );
        }
        h.ra = null;
        h.I = g.k;
        null != h.I && (h.I.ra = h);
        g.k = h;
      }
    for (k = 1; k <= a.g; k++)
      for (f = a.n[k], h = f.k; null != h; h = b)
        (b = h.B),
          0 == h.j &&
            (null == h.ua ? (f.k = b) : (h.ua.B = b),
            null != b && (b.ua = h.ua),
            null == h.ra ? (h.f.k = h.I) : (h.ra.I = h.I),
            null != h.I && (h.I.ra = h.ra),
            a.L--);
    a.valid = 0;
  };
  exports.glp_check_dup = function (a, b, c, d, e) {
    var f, g, h, k, l;
    0 > a && w("glp_check_dup: m = %d; invalid parameter");
    0 > b && w("glp_check_dup: n = %d; invalid parameter");
    0 > c && w("glp_check_dup: ne = %d; invalid parameter");
    0 < c && null == d && w("glp_check_dup: ia = " + d + "; invalid parameter");
    0 < c && null == e && w("glp_check_dup: ja = " + e + "; invalid parameter");
    for (h = 1; h <= c; h++)
      if (((f = d[h]), (g = e[h]), !(1 <= f && f <= a && 1 <= g && g <= b)))
        return (a = -h);
    if (0 == a || 0 == b) return 0;
    k = new Int32Array(1 + a);
    l = new Int32Array(1 + c);
    b = new Int8Array(1 + b);
    for (h = 1; h <= c; h++) (f = d[h]), (l[h] = k[f]), (k[f] = h);
    for (f = 1; f <= a; f++) {
      for (h = k[f]; 0 != h; h = l[h]) {
        g = e[h];
        if (b[g]) {
          for (h = 1; h <= c && (d[h] != f || e[h] != g); h++);
          for (h++; h <= c && (d[h] != f || e[h] != g); h++);
          return (a = +h);
        }
        b[g] = 1;
      }
      for (h = k[f]; 0 != h; h = l[h]) b[e[h]] = 0;
    }
    return 0;
  };
  var $a = (exports.glp_sort_matrix = function (a) {
      var b, c, d;
      (null != a && 3621377730 == a.Dd) ||
        w("glp_sort_matrix: P = " + a + "; invalid problem object");
      for (c = a.g; 1 <= c; c--) a.n[c].k = null;
      for (d = a.i; 1 <= d; d--)
        for (b = a.f[d].k; null != b; b = b.I)
          (c = b.n.ea),
            (b.ua = null),
            (b.B = a.n[c].k),
            null != b.B && (b.B.ua = b),
            (a.n[c].k = b);
      for (d = a.i; 1 <= d; d--) a.f[d].k = null;
      for (c = a.g; 1 <= c; c--)
        for (b = a.n[c].k; null != b; b = b.B)
          (d = b.f.C),
            (b.ra = null),
            (b.I = a.f[d].k),
            null != b.I && (b.I.ra = b),
            (a.f[d].k = b);
    }),
    ab = (exports.glp_del_rows = function (a, b, c) {
      var d = a.V,
        e,
        f,
        g;
      (1 <= b && b <= a.g) ||
        w("glp_del_rows: nrs = " + b + "; invalid number of rows");
      for (g = 1; g <= b; g++)
        (f = c[g]),
          (1 <= f && f <= a.g) ||
            w(
              "glp_del_rows: num[" +
                g +
                "] = " +
                f +
                "; row number out of range"
            ),
          (e = a.n[f]),
          null != d &&
            0 != d.reason &&
            (d.reason != Ga &&
              d.reason != Ia &&
              w("glp_del_rows: operation not allowed"),
            e.La != d.N.La &&
              w(
                "glp_del_rows: num[" +
                  g +
                  "] = " +
                  f +
                  "; invalid attempt to delete row created not in current subproblem"
              ),
            e.m != A &&
              w(
                "glp_del_rows: num[" +
                  g +
                  "] = " +
                  f +
                  "; invalid attempt to delete active row (constraint)"
              ),
            (d.tf = 1)),
          0 == e.ea &&
            w(
              "glp_del_rows: num[" +
                g +
                "] = " +
                f +
                "; duplicate row numbers not allowed"
            ),
          Pa(a, f, null),
          Ya(a, f, 0, null, null),
          (e.ea = 0);
      b = 0;
      for (f = 1; f <= a.g; f++)
        (e = a.n[f]), 0 != e.ea && ((e.ea = ++b), (a.n[e.ea] = e));
      a.g = b;
      a.valid = 0;
    });
  exports.glp_del_cols = function (a, b, c) {
    var d = a.V,
      e,
      f;
    null != d && 0 != d.reason && w("glp_del_cols: operation not allowed");
    (1 <= b && b <= a.i) ||
      w("glp_del_cols: ncs = " + b + "; invalid number of columns");
    for (f = 1; f <= b; f++)
      (d = c[f]),
        (1 <= d && d <= a.i) ||
          w(
            "glp_del_cols: num[" +
              f +
              "] = " +
              d +
              "; column number out of range"
          ),
        (e = a.f[d]),
        0 == e.C &&
          w(
            "glp_del_cols: num[" +
              f +
              "] = " +
              d +
              "; duplicate column numbers not allowed"
          ),
        Qa(a, d, null),
        Za(a, d, 0, null, null),
        (e.C = 0),
        e.m == A && (a.valid = 0);
    b = 0;
    for (d = 1; d <= a.i; d++)
      (e = a.f[d]), 0 != e.C && ((e.C = ++b), (a.f[e.C] = e));
    a.i = b;
    if (a.valid)
      for (c = a.g, e = a.head, d = 1; d <= b; d++)
        (f = a.f[d].bind), 0 != f && (e[f] = c + d);
  };
  var hb = (exports.glp_copy_prob = function (a, b, c) {
      var d = a.V,
        e = {},
        f,
        g,
        h,
        k;
      null != d && 0 != d.reason && w("glp_copy_prob: operation not allowed");
      a == b &&
        w("glp_copy_prob: copying problem object to itself not allowed");
      c != bb &&
        c != cb &&
        w("glp_copy_prob: names = " + c + "; invalid parameter");
      db(a);
      c && null != b.name && Ca(a, b.name);
      c && null != b.eb && Da(a, b.eb);
      a.dir = b.dir;
      a.ha = b.ha;
      0 < b.g && La(a, b.g);
      0 < b.i && Oa(a, b.i);
      eb(b, e);
      fb(a, e);
      a.na = b.na;
      a.sa = b.sa;
      a.aa = b.aa;
      a.some = b.some;
      a.df = b.df;
      a.ae = b.ae;
      a.za = b.za;
      a.ta = b.ta;
      for (f = 1; f <= b.g; f++)
        (d = a.n[f]),
          (e = b.n[f]),
          c && null != e.name && Pa(a, f, e.name),
          (d.type = e.type),
          (d.c = e.c),
          (d.d = e.d),
          (d.ma = e.ma),
          (d.m = e.m),
          (d.r = e.r),
          (d.J = e.J),
          (d.dc = e.dc),
          (d.mc = e.mc),
          (d.Sa = e.Sa);
      h = new Int32Array(1 + b.g);
      k = new Float64Array(1 + b.g);
      for (f = 1; f <= b.i; f++)
        (d = a.f[f]),
          (e = b.f[f]),
          c && null != e.name && Qa(a, f, e.name),
          (d.kind = e.kind),
          (d.type = e.type),
          (d.c = e.c),
          (d.d = e.d),
          (d.u = e.u),
          (g = gb(b, f, h, k)),
          Za(a, f, g, h, k),
          (d.va = e.va),
          (d.m = e.m),
          (d.r = e.r),
          (d.J = e.J),
          (d.dc = e.dc),
          (d.mc = e.mc),
          (d.Sa = e.Sa);
    }),
    db = (exports.glp_erase_prob = function (a) {
      var b = a.V;
      null != b && 0 != b.reason && w("glp_erase_prob: operation not allowed");
      a.Dd = 1061109567;
      a.ke = null;
      a.n = null;
      a.f = null;
      a.fc = null;
      a.Lc = null;
      a.head = null;
      a.Qd = null;
      a.U = null;
      ya(a);
    });
  exports.glp_get_prob_name = function (a) {
    return a.name;
  };
  var ib = (exports.glp_get_obj_name = function (a) {
      return a.eb;
    }),
    jb = (exports.glp_get_obj_dir = function (a) {
      return a.dir;
    }),
    kb = (exports.glp_get_num_rows = function (a) {
      return a.g;
    }),
    lb = (exports.glp_get_num_cols = function (a) {
      return a.i;
    }),
    mb = (exports.glp_get_row_name = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_row_name: i = " + b + "; row number out of range");
      return a.n[b].name;
    }),
    nb = (exports.glp_get_col_name = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_name: j = " + b + "; column number out of range");
      return a.f[b].name;
    }),
    ob = (exports.glp_get_row_type = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_row_type: i = " + b + "; row number out of range");
      return a.n[b].type;
    }),
    pb = (exports.glp_get_row_lb = function (a, b) {
      var c;
      (1 <= b && b <= a.g) ||
        w("glp_get_row_lb: i = " + b + "; row number out of range");
      switch (a.n[b].type) {
        case Ka:
        case Ta:
          c = -s;
          break;
        case Sa:
        case I:
        case B:
          c = a.n[b].c;
      }
      return c;
    }),
    qb = (exports.glp_get_row_ub = function (a, b) {
      var c;
      (1 <= b && b <= a.g) ||
        w("glp_get_row_ub: i = " + b + "; row number out of range");
      switch (a.n[b].type) {
        case Ka:
        case Sa:
          c = +s;
          break;
        case Ta:
        case I:
        case B:
          c = a.n[b].d;
      }
      return c;
    }),
    rb = (exports.glp_get_col_type = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_type: j = " + b + "; column number out of range");
      return a.f[b].type;
    }),
    sb = (exports.glp_get_col_lb = function (a, b) {
      var c;
      (1 <= b && b <= a.i) ||
        w("glp_get_col_lb: j = " + b + "; column number out of range");
      switch (a.f[b].type) {
        case Ka:
        case Ta:
          c = -s;
          break;
        case Sa:
        case I:
        case B:
          c = a.f[b].c;
      }
      return c;
    }),
    tb = (exports.glp_get_col_ub = function (a, b) {
      var c;
      (1 <= b && b <= a.i) ||
        w("glp_get_col_ub: j = " + b + "; column number out of range");
      switch (a.f[b].type) {
        case Ka:
        case Sa:
          c = +s;
          break;
        case Ta:
        case I:
        case B:
          c = a.f[b].d;
      }
      return c;
    });
  exports.glp_get_obj_coef = function (a, b) {
    (0 <= b && b <= a.i) ||
      w("glp_get_obj_coef: j = " + b + "; column number out of range");
    return 0 == b ? a.ha : a.f[b].u;
  };
  exports.glp_get_num_nz = function (a) {
    return a.L;
  };
  var ub = (exports.glp_get_mat_row = function (a, b, c, d) {
      var e;
      (1 <= b && b <= a.g) ||
        w("glp_get_mat_row: i = " + b + "; row number out of range");
      e = 0;
      for (a = a.n[b].k; null != a; a = a.B)
        e++, null != c && (c[e] = a.f.C), null != d && (d[e] = a.j);
      return e;
    }),
    gb = (exports.glp_get_mat_col = function (a, b, c, d) {
      var e;
      (1 <= b && b <= a.i) ||
        w("glp_get_mat_col: j = " + b + "; column number out of range");
      e = 0;
      for (a = a.f[b].k; null != a; a = a.I)
        e++, null != c && (c[e] = a.n.ea), null != d && (d[e] = a.j);
      return e;
    }),
    vb = (exports.glp_create_index = function (a) {
      var b, c;
      if (null == a.fc)
        for (a.fc = {}, c = 1; c <= a.g; c++)
          (b = a.n[c]), null != b.name && (a.fc[b.name] = b);
      if (null == a.Lc)
        for (a.Lc = {}, c = 1; c <= a.i; c++)
          (b = a.f[c]), null != b.name && (a.Lc[b.name] = b);
    }),
    wb = (exports.glp_find_row = function (a, b) {
      var c = 0;
      null == a.fc && w("glp_find_row: row name index does not exist");
      var d = a.fc[b];
      d && (c = d.ea);
      return c;
    }),
    xb = (exports.glp_find_col = function (a, b) {
      var c = 0;
      null == a.Lc && w("glp_find_col: column name index does not exist");
      var d = a.Lc[b];
      d && (c = d.C);
      return c;
    }),
    yb = (exports.glp_delete_index = function (a) {
      a.fc = null;
      a.fc = null;
    }),
    zb = (exports.glp_set_rii = function (a, b, c) {
      (1 <= b && b <= a.g) ||
        w("glp_set_rii: i = " + b + "; row number out of range");
      0 >= c &&
        w("glp_set_rii: i = " + b + "; rii = " + c + "; invalid scale factor");
      if (a.valid && a.n[b].ma != c)
        for (var d = a.n[b].k; null != d; d = d.B)
          if (d.f.m == A) {
            a.valid = 0;
            break;
          }
      a.n[b].ma = c;
    }),
    Ab = (exports.glp_set_sjj = function (a, b, c) {
      (1 <= b && b <= a.i) ||
        w("glp_set_sjj: j = " + b + "; column number out of range");
      0 >= c &&
        w("glp_set_sjj: j = " + b + "; sjj = " + c + "; invalid scale factor");
      a.valid && a.f[b].va != c && a.f[b].m == A && (a.valid = 0);
      a.f[b].va = c;
    }),
    Cb = (exports.glp_get_rii = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_rii: i = " + b + "; row number out of range");
      return a.n[b].ma;
    }),
    Db = (exports.glp_get_sjj = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_sjj: j = " + b + "; column number out of range");
      return a.f[b].va;
    }),
    Eb = (exports.glp_unscale_prob = function (a) {
      var b = kb(a),
        c = lb(a),
        d;
      for (d = 1; d <= b; d++) zb(a, d, 1);
      for (b = 1; b <= c; b++) Ab(a, b, 1);
    }),
    Fb = (exports.glp_set_row_stat = function (a, b, c) {
      (1 <= b && b <= a.g) ||
        w("glp_set_row_stat: i = " + b + "; row number out of range");
      c != A &&
        c != G &&
        c != Ua &&
        c != Ra &&
        c != Na &&
        w("glp_set_row_stat: i = " + b + "; stat = " + c + "; invalid status");
      b = a.n[b];
      if (c != A)
        switch (b.type) {
          case Ka:
            c = Ra;
            break;
          case Sa:
            c = G;
            break;
          case Ta:
            c = Ua;
            break;
          case I:
            c != Ua && (c = G);
            break;
          case B:
            c = Na;
        }
      if ((b.m == A && c != A) || (b.m != A && c == A)) a.valid = 0;
      b.m = c;
    }),
    Gb = (exports.glp_set_col_stat = function (a, b, c) {
      (1 <= b && b <= a.i) ||
        w("glp_set_col_stat: j = " + b + "; column number out of range");
      c != A &&
        c != G &&
        c != Ua &&
        c != Ra &&
        c != Na &&
        w("glp_set_col_stat: j = " + b + "; stat = " + c + "; invalid status");
      b = a.f[b];
      if (c != A)
        switch (b.type) {
          case Ka:
            c = Ra;
            break;
          case Sa:
            c = G;
            break;
          case Ta:
            c = Ua;
            break;
          case I:
            c != Ua && (c = G);
            break;
          case B:
            c = Na;
        }
      if ((b.m == A && c != A) || (b.m != A && c == A)) a.valid = 0;
      b.m = c;
    }),
    Hb = (exports.glp_std_basis = function (a) {
      var b;
      for (b = 1; b <= a.g; b++) Fb(a, b, A);
      for (b = 1; b <= a.i; b++) {
        var c = a.f[b];
        c.type == I && Math.abs(c.c) > Math.abs(c.d)
          ? Gb(a, b, Ua)
          : Gb(a, b, G);
      }
    }),
    sc = (exports.glp_simplex = function (a, b) {
      function c(a, b) {
        var c;
        if (
          !Ib(a) &&
          ((c = Jb(a)),
          0 != c &&
            (c == Kb
              ? b.o >= Lb && x("glp_simplex: initial basis is invalid")
              : c == Mb
              ? b.o >= Lb && x("glp_simplex: initial basis is singular")
              : c == Nb &&
                b.o >= Lb &&
                x("glp_simplex: initial basis is ill-conditioned")),
          0 != c)
        )
          return c;
        b.cb == Ob
          ? (c = Pb(a, b))
          : b.cb == Qb
          ? ((c = Rb(a, b)), c == Sb && a.valid && (c = Pb(a, b)))
          : b.cb == Tb && (c = Rb(a, b));
        return c;
      }
      function d(a, b) {
        function d() {
          Ub(e, f);
          f = null;
          Vb(e, a);
          return (r = 0);
        }
        var e,
          f = null,
          g = {},
          r;
        b.o >= Wb && x("Preprocessing...");
        e = Xb();
        Yb(e, a, Zb);
        r = $b(e, 0);
        0 != r &&
          (r == ac
            ? b.o >= Wb && x("PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION")
            : r == bc &&
              b.o >= Wb &&
              x("PROBLEM HAS NO DUAL FEASIBLE SOLUTION"));
        if (0 != r) return r;
        f = Ba();
        cc(e, f);
        if (0 == f.g && 0 == f.i)
          return (
            (f.na = f.sa = dc),
            (f.aa = f.ha),
            b.o >= fc &&
              0 == b.fb &&
              x(a.$ + ": obj = " + f.aa + "  infeas = 0.0"),
            b.o >= Wb && x("OPTIMAL SOLUTION FOUND BY LP PREPROCESSOR"),
            d()
          );
        b.o >= Wb &&
          x(
            f.g +
              " row" +
              (1 == f.g ? "" : "s") +
              ", " +
              f.i +
              " column" +
              (1 == f.i ? "" : "s") +
              ", " +
              f.L +
              " non-zero" +
              (1 == f.L ? "" : "s") +
              ""
          );
        eb(a, g);
        fb(f, g);
        var g = na,
          n = g.Fb;
        g.Fb = !n || b.o < Wb ? cb : bb;
        gc(f, hc);
        g.Fb = n;
        g = na;
        n = g.Fb;
        g.Fb = !n || b.o < Wb ? cb : bb;
        ic(f);
        g.Fb = n;
        f.$ = a.$;
        r = c(f, b);
        a.$ = f.$;
        return 0 != r || f.na != dc || f.sa != dc
          ? (b.o >= Lb &&
              x(
                "glp_simplex: unable to recover undefined or non-optimal solution"
              ),
            0 == r && (f.na == jc ? (r = ac) : f.sa == jc && (r = bc)),
            r)
          : d();
      }
      function e(a, b) {
        function c() {
          f.m = G;
          f.r = f.c;
        }
        function d() {
          f.m = Ua;
          f.r = f.d;
        }
        var e, f, g, n, t;
        a.valid = 0;
        a.na = a.sa = dc;
        a.aa = a.ha;
        n = t = a.some = 0;
        for (g = 1; g <= a.g; g++) {
          e = a.n[g];
          e.m = A;
          e.r = e.J = 0;
          if (e.type == Sa || e.type == I || e.type == B)
            e.c > +b.Gb &&
              ((a.na = jc), 0 == a.some && b.cb != Ob && (a.some = g)),
              n < +e.c && (n = +e.c);
          if (e.type == Ta || e.type == I || e.type == B)
            e.d < -b.Gb &&
              ((a.na = jc), 0 == a.some && b.cb != Ob && (a.some = g)),
              n < -e.d && (n = -e.d);
        }
        for (e = g = 1; e <= a.i; e++)
          (f = a.f[e]), g < Math.abs(f.u) && (g = Math.abs(f.u));
        g = (a.dir == za ? 1 : -1) / g;
        for (e = 1; e <= a.i; e++) {
          f = a.f[e];
          f.type == Ka
            ? ((f.m = Ra), (f.r = 0))
            : f.type == Sa
            ? c()
            : f.type == Ta
            ? d()
            : f.type == I
            ? 0 < g * f.u
              ? c()
              : 0 > g * f.u
              ? d()
              : Math.abs(f.c) <= Math.abs(f.d)
              ? c()
              : d()
            : f.type == B && ((f.m = Na), (f.r = f.c));
          f.J = f.u;
          a.aa += f.u * f.r;
          if (f.type == Ka || f.type == Sa)
            g * f.J < -b.tb &&
              ((a.sa = jc), 0 == a.some && b.cb == Ob && (a.some = a.g + e)),
              t < -g * f.J && (t = -g * f.J);
          if (f.type == Ka || f.type == Ta)
            g * f.J > +b.tb &&
              ((a.sa = jc), 0 == a.some && b.cb == Ob && (a.some = a.g + e)),
              t < +g * f.J && (t = +g * f.J);
        }
        b.o >= fc &&
          0 == b.fb &&
          x(
            "~" +
              a.$ +
              ": obj = " +
              a.aa +
              "  infeas = " +
              (b.cb == Ob ? n : t) +
              ""
          );
        b.o >= Wb &&
          0 == b.fb &&
          (a.na == dc && a.sa == dc
            ? x("OPTIMAL SOLUTION FOUND")
            : a.na == jc
            ? x("PROBLEM HAS NO FEASIBLE SOLUTION")
            : b.cb == Ob
            ? x("PROBLEM HAS UNBOUNDED SOLUTION")
            : x("PROBLEM HAS NO DUAL FEASIBLE SOLUTION"));
      }
      var f;
      (null != a && 3621377730 == a.Dd) ||
        w("glp_simplex: P = " + a + "; invalid problem object");
      null != a.V && 0 != a.V.reason && w("glp_simplex: operation not allowed");
      null == b && (b = new kc());
      b.o != lc &&
        b.o != Lb &&
        b.o != fc &&
        b.o != Wb &&
        b.o != mc &&
        w("glp_simplex: msg_lev = " + b.o + "; invalid parameter");
      b.cb != Ob &&
        b.cb != Qb &&
        b.cb != Tb &&
        w("glp_simplex: meth = " + b.cb + "; invalid parameter");
      b.fd != nc &&
        b.fd != oc &&
        w("glp_simplex: pricing = " + b.fd + "; invalid parameter");
      b.ne != pc &&
        b.ne != qc &&
        w("glp_simplex: r_test = " + b.ne + "; invalid parameter");
      (0 < b.Gb && 1 > b.Gb) ||
        w("glp_simplex: tol_bnd = " + b.Gb + "; invalid parameter");
      (0 < b.tb && 1 > b.tb) ||
        w("glp_simplex: tol_dj = " + b.tb + "; invalid parameter");
      (0 < b.xe && 1 > b.xe) ||
        w("glp_simplex: tol_piv = " + b.xe + "; invalid parameter");
      0 > b.oc && w("glp_simplex: it_lim = " + b.oc + "; invalid parameter");
      0 > b.sb && w("glp_simplex: tm_lim = " + b.sb + "; invalid parameter");
      1 > b.bc && w("glp_simplex: out_frq = " + b.bc + "; invalid parameter");
      0 > b.fb && w("glp_simplex: out_dly = " + b.fb + "; invalid parameter");
      b.yc != bb &&
        b.yc != cb &&
        w("glp_simplex: presolve = " + b.yc + "; invalid parameter");
      a.na = a.sa = Aa;
      a.aa = 0;
      a.some = 0;
      for (f = 1; f <= a.g; f++) {
        var g = a.n[f];
        if (g.type == I && g.c >= g.d)
          return (
            b.o >= Lb &&
              x(
                "glp_simplex: row " +
                  f +
                  ": lb = " +
                  g.c +
                  ", ub = " +
                  g.d +
                  "; incorrect bounds"
              ),
            (f = rc)
          );
      }
      for (f = 1; f <= a.i; f++)
        if (((g = a.f[f]), g.type == I && g.c >= g.d))
          return (
            b.o >= Lb &&
              x(
                "glp_simplex: column " +
                  f +
                  ": lb = " +
                  g.c +
                  ", ub = " +
                  g.d +
                  "; incorrect bounds"
              ),
            (f = rc)
          );
      b.o >= Wb &&
        (x("GLPK Simplex Optimizer, v" + ra() + ""),
        x(
          a.g +
            " row" +
            (1 == a.g ? "" : "s") +
            ", " +
            a.i +
            " column" +
            (1 == a.i ? "" : "s") +
            ", " +
            a.L +
            " non-zero" +
            (1 == a.L ? "" : "s") +
            ""
        ));
      0 == a.L ? (e(a, b), (f = 0)) : (f = b.yc ? d(a, b) : c(a, b));
      return f;
    }),
    kc = (exports.SMCP = function (a) {
      a = a || {};
      this.o = a.msg_lev || Wb;
      this.cb = a.meth || Ob;
      this.fd = a.pricing || oc;
      this.ne = a.r_test || qc;
      this.Gb = a.tol_bnd || 1e-7;
      this.tb = a.tol_dj || 1e-7;
      this.xe = a.tol_piv || 1e-10;
      this.hf = a.obj_ll || -s;
      this.jf = a.obj_ul || +s;
      this.oc = a.it_lim || 2147483647;
      this.sb = a.tm_lim || 2147483647;
      this.bc = a.out_frq || 500;
      this.fb = a.out_dly || 0;
      this.yc = a.presolve || cb;
    }),
    xc = (exports.glp_get_status = function (a) {
      var b;
      b = tc(a);
      switch (b) {
        case dc:
          switch (uc(a)) {
            case dc:
              b = vc;
              break;
            case jc:
              b = wc;
          }
      }
      return b;
    }),
    tc = (exports.glp_get_prim_stat = function (a) {
      return a.na;
    }),
    uc = (exports.glp_get_dual_stat = function (a) {
      return a.sa;
    }),
    yc = (exports.glp_get_obj_val = function (a) {
      return a.aa;
    }),
    zc = (exports.glp_get_row_stat = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_row_stat: i = " + b + "; row number out of range");
      return a.n[b].m;
    }),
    Ac = (exports.glp_get_row_prim = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_row_prim: i = " + b + "; row number out of range");
      return a.n[b].r;
    }),
    Bc = (exports.glp_get_row_dual = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_get_row_dual: i = " + b + "; row number out of range");
      return a.n[b].J;
    }),
    Cc = (exports.glp_get_col_stat = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_stat: j = " + b + "; column number out of range");
      return a.f[b].m;
    }),
    Dc = (exports.glp_get_col_prim = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_prim: j = " + b + "; column number out of range");
      return a.f[b].r;
    }),
    Ec = (exports.glp_get_col_dual = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_dual: j = " + b + "; column number out of range");
      return a.f[b].J;
    });
  exports.glp_get_unbnd_ray = function (a) {
    var b = a.some;
    b > a.g + a.i && (b = 0);
    return b;
  };
  var Hc = (exports.glp_set_col_kind = function (a, b, c) {
      (1 <= b && b <= a.i) ||
        w("glp_set_col_kind: j = " + b + "; column number out of range");
      var d = a.f[b];
      switch (c) {
        case Ma:
          d.kind = Ma;
          break;
        case Fc:
          d.kind = Fc;
          break;
        case Gc:
          d.kind = Fc;
          (d.type == I && 0 == d.c && 1 == d.d) || Wa(a, b, I, 0, 1);
          break;
        default:
          w(
            "glp_set_col_kind: j = " +
              b +
              "; kind = " +
              c +
              "; invalid column kind"
          );
      }
    }),
    Ic = (exports.glp_get_col_kind = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_get_col_kind: j = " + b + "; column number out of range");
      var c = a.f[b],
        d = c.kind;
      switch (d) {
        case Fc:
          c.type == I && 0 == c.c && 1 == c.d && (d = Gc);
      }
      return d;
    }),
    Jc = (exports.glp_get_num_int = function (a) {
      for (var b, c = 0, d = 1; d <= a.i; d++)
        (b = a.f[d]), b.kind == Fc && c++;
      return c;
    }),
    Kc = (exports.glp_get_num_bin = function (a) {
      for (var b, c = 0, d = 1; d <= a.i; d++)
        (b = a.f[d]),
          b.kind == Fc && b.type == I && 0 == b.c && 1 == b.d && c++;
      return c;
    });
  exports.glp_intopt = function (a, b) {
    function c(a, b) {
      var c;
      if (xc(a) != vc)
        return (
          b.o >= Lb &&
            x(
              "glp_intopt: optimal basis to initial LP relaxation not provided"
            ),
          (c = Lc)
        );
      b.o >= Wb && x("Integer optimization begins...");
      var d = a.g;
      c = a.i;
      var e, f;
      a.V = e = {};
      e.i = c;
      e.wc = d;
      e.ac = new Int8Array(1 + d + c);
      e.bd = new Float64Array(1 + d + c);
      e.cd = new Float64Array(1 + d + c);
      e.mf = new Int8Array(1 + d + c);
      e.lf = new Float64Array(1 + d + c);
      e.kf = new Float64Array(1 + d + c);
      for (f = 1; f <= d; f++) {
        var q = a.n[f];
        e.ac[f] = q.type;
        e.bd[f] = q.c;
        e.cd[f] = q.d;
        e.mf[f] = q.m;
        e.lf[f] = q.r;
        e.kf[f] = q.J;
      }
      for (f = 1; f <= c; f++)
        (q = a.f[f]),
          (e.ac[d + f] = q.type),
          (e.bd[d + f] = q.c),
          (e.cd[d + f] = q.d),
          (e.mf[d + f] = q.m),
          (e.lf[d + f] = q.r),
          (e.kf[d + f] = q.J);
      e.ih = a.aa;
      e.fe = 0;
      e.Sc = 0;
      e.ya = null;
      e.head = e.Xa = null;
      e.Pd = e.Zf = e.Jg = 0;
      e.Ig = 0;
      e.se = null;
      e.qe = e.te = null;
      e.re = null;
      e.N = null;
      e.A = a;
      e.ad = new Int8Array(1 + c);
      e.Fg = e.Gg = 0;
      e.qf = null;
      e.of = e.rf = null;
      e.pf = null;
      d = { size: 0 };
      d.head = d.Xa = null;
      d.fh = 0;
      d.N = null;
      e.Bd = d;
      e.Xf = null;
      e.Ne = null;
      e.Fd = null;
      e.Bg = new Int32Array(1 + c);
      e.Tg = new Float64Array(1 + c);
      e.p = b;
      e.hc = ja();
      e.Lg = 0;
      e.qh = 0;
      e.reason = 0;
      e.pe = 0;
      e.tf = 0;
      e.Tc = 0;
      e.Jf = 0;
      e.ud = 0;
      e.gf = 0;
      e.stop = 0;
      Mc(e, null);
      c = Nc(e);
      var d = e.A,
        r = d.g;
      f = d.i;
      if (r != e.wc) {
        var n,
          r = r - e.wc;
        n = new Int32Array(1 + r);
        for (q = 1; q <= r; q++) n[q] = e.wc + q;
        ab(d, r, n);
      }
      r = e.wc;
      for (q = 1; q <= r; q++)
        Va(d, q, e.ac[q], e.bd[q], e.cd[q]),
          Fb(d, q, e.mf[q]),
          (d.n[q].r = e.lf[q]),
          (d.n[q].J = e.kf[q]);
      for (q = 1; q <= f; q++)
        Wa(d, q, e.ac[r + q], e.bd[r + q], e.cd[r + q]),
          Gb(d, q, e.mf[r + q]),
          (d.f[q].r = e.lf[r + q]),
          (d.f[q].J = e.kf[r + q]);
      d.na = d.sa = dc;
      d.aa = e.ih;
      Oc(e.Bd);
      d.V = null;
      0 == c
        ? a.za == dc
          ? (b.o >= Wb && x("INTEGER OPTIMAL SOLUTION FOUND"), (a.za = vc))
          : (b.o >= Wb && x("PROBLEM HAS NO INTEGER FEASIBLE SOLUTION"),
            (a.za = jc))
        : c == Pc
        ? b.o >= Wb &&
          x("RELATIVE MIP GAP TOLERANCE REACHED; SEARCH TERMINATED")
        : c == Qc
        ? b.o >= Wb && x("TIME LIMIT EXCEEDED; SEARCH TERMINATED")
        : c == Sb
        ? b.o >= Lb && x("glp_intopt: cannot solve current LP relaxation")
        : c == Rc && b.o >= Wb && x("SEARCH TERMINATED BY APPLICATION");
      return c;
    }
    function d(a, b) {
      function d() {
        Ub(m, q);
        q = null;
        Vb(m, a);
        return n;
      }
      var e = na,
        f = e.Fb,
        m,
        q = null,
        r = {},
        n;
      b.o >= Wb && x("Preprocessing...");
      m = Xb();
      Yb(m, a, Sc);
      e.Fb = !f || b.o < Wb ? cb : bb;
      n = Tc(m, b);
      e.Fb = f;
      0 != n &&
        (n == ac
          ? b.o >= Wb && x("PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION")
          : n == bc &&
            b.o >= Wb &&
            x("LP RELAXATION HAS NO DUAL FEASIBLE SOLUTION"));
      if (0 != n) return n;
      q = Ba();
      cc(m, q);
      if (0 == q.g && 0 == q.i)
        return (
          (q.za = vc),
          (q.ta = q.ha),
          b.o >= Wb &&
            (x("Objective value = " + q.ta + ""),
            x("INTEGER OPTIMAL SOLUTION FOUND BY MIP PREPROCESSOR")),
          d()
        );
      if (b.o >= Wb) {
        var t = Jc(q),
          y = Kc(q);
        x(
          q.g +
            " row" +
            (1 == q.g ? "" : "s") +
            ", " +
            q.i +
            " column" +
            (1 == q.i ? "" : "s") +
            ", " +
            q.L +
            " non-zero" +
            (1 == q.L ? "" : "s") +
            ""
        );
        x(
          t +
            " integer variable" +
            (1 == t ? "" : "s") +
            ", " +
            (0 == y
              ? "none of"
              : 1 == t && 1 == y
              ? ""
              : 1 == y
              ? "one of"
              : y == t
              ? "all of"
              : y + " of") +
            " which " +
            (1 == y ? "is" : "are") +
            " binary"
        );
      }
      eb(a, r);
      fb(q, r);
      e.Fb = !f || b.o < Wb ? cb : bb;
      gc(q, Uc | Vc | Wc | Xc);
      e.Fb = f;
      e.Fb = !f || b.o < Wb ? cb : bb;
      ic(q);
      e.Fb = f;
      b.o >= Wb && x("Solving LP relaxation...");
      e = new kc();
      e.o = b.o;
      q.$ = a.$;
      n = sc(q, e);
      a.$ = q.$;
      if (0 != n)
        return (
          b.o >= Lb && x("glp_intopt: cannot solve LP relaxation"), (n = Sb)
        );
      n = xc(q);
      n == vc ? (n = 0) : n == jc ? (n = ac) : n == wc && (n = bc);
      if (0 != n) return n;
      q.$ = a.$;
      n = c(q, b);
      a.$ = q.$;
      return q.za != vc && q.za != dc ? ((a.za = q.za), n) : d();
    }
    var e, f;
    (null != a && 3621377730 == a.Dd) ||
      w("glp_intopt: P = " + a + "; invalid problem object");
    null != a.V && w("glp_intopt: operation not allowed");
    null == b && (b = new Yc());
    b.o != lc &&
      b.o != Lb &&
      b.o != fc &&
      b.o != Wb &&
      b.o != mc &&
      w("glp_intopt: msg_lev = " + b.o + "; invalid parameter");
    b.Jb != Zc &&
      b.Jb != $c &&
      b.Jb != ad &&
      b.Jb != bd &&
      b.Jb != cd &&
      w("glp_intopt: br_tech = " + b.Jb + "; invalid parameter");
    b.kc != dd &&
      b.kc != ed &&
      b.kc != fd &&
      b.kc != gd &&
      w("glp_intopt: bt_tech = " + b.kc + "; invalid parameter");
    (0 < b.Ub && 1 > b.Ub) ||
      w("glp_intopt: tol_int = " + b.Ub + "; invalid parameter");
    (0 < b.we && 1 > b.we) ||
      w("glp_intopt: tol_obj = " + b.we + "; invalid parameter");
    0 > b.sb && w("glp_intopt: tm_lim = " + b.sb + "; invalid parameter");
    0 > b.bc && w("glp_intopt: out_frq = " + b.bc + "; invalid parameter");
    0 > b.fb && w("glp_intopt: out_dly = " + b.fb + "; invalid parameter");
    (0 <= b.Me && 256 >= b.Me) ||
      w("glp_intopt: cb_size = " + b.Me + "; invalid parameter");
    b.ed != hd &&
      b.ed != id &&
      b.ed != jd &&
      w("glp_intopt: pp_tech = " + b.ed + "; invalid parameter");
    0 > b.ce && w("glp_intopt: mip_gap = " + b.ce + "; invalid parameter");
    b.Ed != bb &&
      b.Ed != cb &&
      w("glp_intopt: mir_cuts = " + b.Ed + "; invalid parameter");
    b.Ad != bb &&
      b.Ad != cb &&
      w("glp_intopt: gmi_cuts = " + b.Ad + "; invalid parameter");
    b.xd != bb &&
      b.xd != cb &&
      w("glp_intopt: cov_cuts = " + b.xd + "; invalid parameter");
    b.vd != bb &&
      b.vd != cb &&
      w("glp_intopt: clq_cuts = " + b.vd + "; invalid parameter");
    b.yc != bb &&
      b.yc != cb &&
      w("glp_intopt: presolve = " + b.yc + "; invalid parameter");
    b.sd != bb &&
      b.sd != cb &&
      w("glp_intopt: binarize = " + b.sd + "; invalid parameter");
    b.Xe != bb &&
      b.Xe != cb &&
      w("glp_intopt: fp_heur = " + b.Xe + "; invalid parameter");
    a.za = Aa;
    a.ta = 0;
    for (e = 1; e <= a.g; e++)
      if (((f = a.n[e]), f.type == I && f.c >= f.d))
        return (
          b.o >= Lb &&
            x(
              "glp_intopt: row " +
                e +
                ": lb = " +
                f.c +
                ", ub = " +
                f.d +
                "; incorrect bounds"
            ),
          (e = rc)
        );
    for (e = 1; e <= a.i; e++)
      if (((f = a.f[e]), f.type == I && f.c >= f.d))
        return (
          b.o >= Lb &&
            x(
              "glp_intopt: column " +
                e +
                ": lb = " +
                f.c +
                ", ub = " +
                f.d +
                "; incorrect bounds"
            ),
          (e = rc)
        );
    for (e = 1; e <= a.i; e++)
      if (((f = a.f[e]), f.kind == Fc)) {
        if ((f.type == Sa || f.type == I) && f.c != Math.floor(f.c))
          return (
            b.o >= Lb &&
              x(
                "glp_intopt: integer column " +
                  e +
                  " has non-integer lower bound " +
                  f.c +
                  ""
              ),
            (e = rc)
          );
        if ((f.type == Ta || f.type == I) && f.d != Math.floor(f.d))
          return (
            b.o >= Lb &&
              x(
                "glp_intopt: integer column " +
                  e +
                  " has non-integer upper bound " +
                  f.d +
                  ""
              ),
            (e = rc)
          );
        if (f.type == B && f.c != Math.floor(f.c))
          return (
            b.o >= Lb &&
              x(
                "glp_intopt: integer column " +
                  e +
                  " has non-integer fixed value " +
                  f.c +
                  ""
              ),
            (e = rc)
          );
      }
    b.o >= Wb &&
      ((e = Jc(a)),
      (f = Kc(a)),
      x("GLPK Integer Optimizer, v" + ra() + ""),
      x(
        a.g +
          " row" +
          (1 == a.g ? "" : "s") +
          ", " +
          a.i +
          " column" +
          (1 == a.i ? "" : "s") +
          ", " +
          a.L +
          " non-zero" +
          (1 == a.L ? "" : "s") +
          ""
      ),
      x(
        e +
          " integer variable" +
          (1 == e ? "" : "s") +
          ", " +
          (0 == f
            ? "none of"
            : 1 == e && 1 == f
            ? ""
            : 1 == f
            ? "one of"
            : f == e
            ? "all of"
            : f + " of") +
          " which " +
          (1 == f ? "is" : "are") +
          " binary"
      ));
    return (e = b.yc ? d(a, b) : c(a, b));
  };
  var Yc = (exports.IOCP = function (a) {
    a = a || {};
    this.o = a.msg_lev || Wb;
    this.Jb = a.br_tech || bd;
    this.kc = a.bt_tech || fd;
    this.Ub = a.tol_int || 1e-5;
    this.we = a.tol_obj || 1e-7;
    this.sb = a.tm_lim || 2147483647;
    this.bc = a.out_frq || 5e3;
    this.fb = a.out_dly || 1e4;
    this.ob = a.cb_func || null;
    this.Uc = a.cb_info || null;
    this.Me = a.cb_size || 0;
    this.ed = a.pp_tech || jd;
    this.ce = a.mip_gap || 0;
    this.Ed = a.mir_cuts || cb;
    this.Ad = a.gmi_cuts || cb;
    this.xd = a.cov_cuts || cb;
    this.vd = a.clq_cuts || cb;
    this.yc = a.presolve || cb;
    this.sd = a.binarize || cb;
    this.Xe = a.fp_heur || cb;
  });
  exports.glp_mip_status = function (a) {
    return a.za;
  };
  exports.glp_mip_obj_val = function (a) {
    return a.ta;
  };
  var kd = (exports.glp_mip_row_val = function (a, b) {
      (1 <= b && b <= a.g) ||
        w("glp_mip_row_val: i = " + b + "; row number out of range");
      return a.n[b].Sa;
    }),
    ld = (exports.glp_mip_col_val = function (a, b) {
      (1 <= b && b <= a.i) ||
        w("glp_mip_col_val: j = " + b + "; column number out of range");
      return a.f[b].Sa;
    }),
    Ib = (exports.glp_bf_exists = function (a) {
      return 0 == a.g || a.valid;
    }),
    Jb = (exports.glp_factorize = function (a) {
      function b(a, b, c, d) {
        var e = a.g,
          f;
        f = a.head[b];
        if (f <= e) (b = 1), (c[1] = f), (d[1] = 1);
        else
          for (b = 0, a = a.f[f - e].k; null != a; a = a.I)
            b++, (c[b] = a.n.ea), (d[b] = -a.n.ma * a.j * a.f.va);
        return b;
      }
      var c = a.g,
        d = a.i,
        e = a.n,
        f = a.f,
        g = a.head,
        h,
        k,
        l;
      h = a.valid = 0;
      for (k = 1; k <= c + d; k++)
        if (
          (k <= c
            ? ((l = e[k].m), (e[k].bind = 0))
            : ((l = f[k - c].m), (f[k - c].bind = 0)),
          l == A)
        ) {
          h++;
          if (h > c) return (a = Kb);
          g[h] = k;
          k <= c ? (e[k].bind = h) : (f[k - c].bind = h);
        }
      if (h < c) return (a = Kb);
      if (0 < c) {
        null == a.U &&
          ((a.U = {
            valid: 0,
            type: md,
            gb: null,
            Za: null,
            Cd: 0,
            cc: 0.1,
            xc: 4,
            gc: 1,
            Mb: 1e-15,
            sc: 1e10,
            $c: 100,
            ic: 1e-6,
            vc: 100,
            kd: 1e3,
            Ch: -1,
            hg: 0,
          }),
          nd(a));
        switch (od(a.U, c, b, a)) {
          case pd:
            return (a = Mb);
          case qd:
            return (a = Nb);
        }
        a.valid = 1;
      }
      return 0;
    });
  exports.glp_bf_updated = function (a) {
    0 == a.g ||
      a.valid ||
      w("glp_bf_update: basis factorization does not exist");
    return 0 == a.g ? 0 : a.U.hg;
  };
  var eb = (exports.glp_get_bfcp = function (a, b) {
    var c = a.Qd;
    null == c
      ? ((b.type = md),
        (b.Cd = 0),
        (b.cc = 0.1),
        (b.xc = 4),
        (b.gc = bb),
        (b.Mb = 1e-15),
        (b.sc = 1e10),
        (b.$c = 100),
        (b.ic = 1e-6),
        (b.vc = 100),
        (b.kd = 0))
      : da(b, c);
  });
  function nd(a) {
    var b = {};
    eb(a, b);
    a = a.U;
    a.type = b.type;
    a.Cd = b.Cd;
    a.cc = b.cc;
    a.xc = b.xc;
    a.gc = b.gc;
    a.Mb = b.Mb;
    a.sc = b.sc;
    a.$c = b.$c;
    a.ic = b.ic;
    a.vc = b.vc;
    a.kd = b.kd;
  }
  var fb = (exports.glp_set_bfcp = function (a, b) {
      var c = a.Qd;
      null == b
        ? null != c && (a.Qd = null)
        : (null == c && (c = a.Qd = {}),
          da(c, b),
          c.type != md &&
            c.type != rd &&
            c.type != sd &&
            w("glp_set_bfcp: type = " + c.type + "; invalid parameter"),
          0 > c.Cd &&
            w("glp_set_bfcp: lu_size = " + c.Cd + "; invalid parameter"),
          (0 < c.cc && 1 > c.cc) ||
            w("glp_set_bfcp: piv_tol = " + c.cc + "; invalid parameter"),
          1 > c.xc &&
            w("glp_set_bfcp: piv_lim = " + c.xc + "; invalid parameter"),
          c.gc != bb &&
            c.gc != cb &&
            w("glp_set_bfcp: suhl = " + c.gc + "; invalid parameter"),
          (0 <= c.Mb && 1e-6 >= c.Mb) ||
            w("glp_set_bfcp: eps_tol = " + c.Mb + "; invalid parameter"),
          1 > c.sc &&
            w("glp_set_bfcp: max_gro = " + c.sc + "; invalid parameter"),
          (1 <= c.$c && 32767 >= c.$c) ||
            w("glp_set_bfcp: nfs_max = " + c.$c + "; invalid parameter"),
          (0 < c.ic && 1 > c.ic) ||
            w("glp_set_bfcp: upd_tol = " + c.ic + "; invalid parameter"),
          (1 <= c.vc && 32767 >= c.vc) ||
            w("glp_set_bfcp: nrs_max = " + c.vc + "; invalid parameter"),
          0 > c.kd &&
            w("glp_set_bfcp: rs_size = " + c.vc + "; invalid parameter"),
          0 == c.kd && (c.kd = 20 * c.vc));
      null != a.U && nd(a);
    }),
    td = (exports.glp_get_bhead = function (a, b) {
      0 == a.g ||
        a.valid ||
        w("glp_get_bhead: basis factorization does not exist");
      (1 <= b && b <= a.g) ||
        w("glp_get_bhead: k = " + b + "; index out of range");
      return a.head[b];
    }),
    ud = (exports.glp_get_row_bind = function (a, b) {
      0 == a.g ||
        a.valid ||
        w("glp_get_row_bind: basis factorization does not exist");
      (1 <= b && b <= a.g) ||
        w("glp_get_row_bind: i = " + b + "; row number out of range");
      return a.n[b].bind;
    }),
    vd = (exports.glp_get_col_bind = function (a, b) {
      0 == a.g ||
        a.valid ||
        w("glp_get_col_bind: basis factorization does not exist");
      (1 <= b && b <= a.i) ||
        w("glp_get_col_bind: j = " + b + "; column number out of range");
      return a.f[b].bind;
    }),
    xd = (exports.glp_ftran = function (a, b) {
      var c = a.g,
        d = a.n,
        e = a.f,
        f,
        g;
      0 == c || a.valid || w("glp_ftran: basis factorization does not exist");
      for (f = 1; f <= c; f++) b[f] *= d[f].ma;
      0 < c && wd(a.U, b);
      for (f = 1; f <= c; f++)
        (g = a.head[f]), (b[f] = g <= c ? b[f] / d[g].ma : b[f] * e[g - c].va);
    }),
    zd = (exports.glp_btran = function (a, b) {
      var c = a.g,
        d = a.n,
        e = a.f,
        f,
        g;
      0 == c || a.valid || w("glp_btran: basis factorization does not exist");
      for (f = 1; f <= c; f++)
        (g = a.head[f]), (b[f] = g <= c ? b[f] / d[g].ma : b[f] * e[g - c].va);
      0 < c && yd(a.U, b);
      for (f = 1; f <= c; f++) b[f] *= d[f].ma;
    });
  exports.glp_warm_up = function (a) {
    var b, c, d, e, f;
    a.na = a.sa = Aa;
    a.aa = 0;
    a.some = 0;
    for (d = 1; d <= a.g; d++) (b = a.n[d]), (b.r = b.J = 0);
    for (d = 1; d <= a.i; d++) (b = a.f[d]), (b.r = b.J = 0);
    if (!Ib(a) && ((e = Jb(a)), 0 != e)) return e;
    e = new Float64Array(1 + a.g);
    for (d = 1; d <= a.g; d++)
      (b = a.n[d]),
        b.m != A &&
          (b.m == G
            ? (b.r = b.c)
            : b.m == Ua
            ? (b.r = b.d)
            : b.m == Ra
            ? (b.r = 0)
            : b.m == Na && (b.r = b.c),
          (e[d] -= b.r));
    for (d = 1; d <= a.i; d++)
      if (
        ((b = a.f[d]),
        b.m != A &&
          (b.m == G
            ? (b.r = b.c)
            : b.m == Ua
            ? (b.r = b.d)
            : b.m == Ra
            ? (b.r = 0)
            : b.m == Na && (b.r = b.c),
          0 != b.r))
      )
        for (c = b.k; null != c; c = c.I) e[c.n.ea] += c.j * b.r;
    xd(a, e);
    a.na = dc;
    for (d = 1; d <= a.g; d++)
      if (((b = a.n[d]), b.m == A)) {
        b.r = e[b.bind];
        c = b.type;
        if (c == Sa || c == I || c == B)
          (f = 1e-6 + 1e-9 * Math.abs(b.c)), b.r < b.c - f && (a.na = Ad);
        if (c == Ta || c == I || c == B)
          (f = 1e-6 + 1e-9 * Math.abs(b.d)), b.r > b.d + f && (a.na = Ad);
      }
    for (d = 1; d <= a.i; d++)
      if (((b = a.f[d]), b.m == A)) {
        b.r = e[b.bind];
        c = b.type;
        if (c == Sa || c == I || c == B)
          (f = 1e-6 + 1e-9 * Math.abs(b.c)), b.r < b.c - f && (a.na = Ad);
        if (c == Ta || c == I || c == B)
          (f = 1e-6 + 1e-9 * Math.abs(b.d)), b.r > b.d + f && (a.na = Ad);
      }
    a.aa = a.ha;
    for (d = 1; d <= a.i; d++) (b = a.f[d]), (a.aa += b.u * b.r);
    for (d = 1; d <= a.g; d++) e[d] = 0;
    for (d = 1; d <= a.i; d++) (b = a.f[d]), b.m == A && (e[b.bind] = b.u);
    zd(a, e);
    a.sa = dc;
    for (d = 1; d <= a.g; d++)
      if (((b = a.n[d]), b.m == A)) b.J = 0;
      else if (
        ((b.J = -e[d]),
        (c = b.m),
        (b = a.dir == za ? +b.J : -b.J),
        ((c == Ra || c == G) && -1e-5 > b) ||
          ((c == Ra || c == Ua) && 1e-5 < b))
      )
        a.sa = Ad;
    for (d = 1; d <= a.i; d++)
      if (((b = a.f[d]), b.m == A)) b.J = 0;
      else {
        b.J = b.u;
        for (c = b.k; null != c; c = c.I) b.J += c.j * e[c.n.ea];
        c = b.m;
        b = a.dir == za ? +b.J : -b.J;
        if (
          ((c == Ra || c == G) && -1e-5 > b) ||
          ((c == Ra || c == Ua) && 1e-5 < b)
        )
          a.sa = Ad;
      }
    return 0;
  };
  var Bd = (exports.glp_eval_tab_row = function (a, b, c, d) {
      var e = a.g,
        f = a.i,
        g,
        h,
        k,
        l,
        p,
        m,
        q;
      0 == e ||
        a.valid ||
        w("glp_eval_tab_row: basis factorization does not exist");
      (1 <= b && b <= e + f) ||
        w("glp_eval_tab_row: k = " + b + "; variable number out of range");
      g = b <= e ? ud(a, b) : vd(a, b - e);
      0 == g && w("glp_eval_tab_row: k = " + b + "; variable must be basic");
      m = new Float64Array(1 + e);
      l = new Int32Array(1 + e);
      q = new Float64Array(1 + e);
      m[g] = 1;
      zd(a, m);
      h = 0;
      for (b = 1; b <= e + f; b++) {
        if (b <= e) {
          if (zc(a, b) == A) continue;
          p = -m[b];
        } else {
          if (Cc(a, b - e) == A) continue;
          k = gb(a, b - e, l, q);
          p = 0;
          for (g = 1; g <= k; g++) p += m[l[g]] * q[g];
        }
        0 != p && (h++, (c[h] = b), (d[h] = p));
      }
      return h;
    }),
    Cd = (exports.glp_eval_tab_col = function (a, b, c, d) {
      var e = a.g,
        f = a.i,
        g;
      0 == e ||
        a.valid ||
        w("glp_eval_tab_col: basis factorization does not exist");
      (1 <= b && b <= e + f) ||
        w("glp_eval_tab_col: k = " + b + "; variable number out of range");
      (b <= e ? zc(a, b) : Cc(a, b - e)) == A &&
        w("glp_eval_tab_col: k = " + b + "; variable must be non-basic");
      f = new Float64Array(1 + e);
      if (b <= e) f[b] = -1;
      else for (g = gb(a, b - e, c, d), b = 1; b <= g; b++) f[c[b]] = d[b];
      xd(a, f);
      g = 0;
      for (b = 1; b <= e; b++)
        0 != f[b] && (g++, (c[g] = td(a, b)), (d[g] = f[b]));
      return g;
    }),
    Dd = (exports.glp_transform_row = function (a, b, c, d) {
      var e, f, g, h, k, l, p, m, q, r;
      Ib(a) || w("glp_transform_row: basis factorization does not exist ");
      f = kb(a);
      g = lb(a);
      m = new Float64Array(1 + g);
      (0 <= b && b <= g) ||
        w("glp_transform_row: len = " + b + "; invalid row length");
      for (h = 1; h <= b; h++)
        (e = c[h]),
          (1 <= e && e <= g) ||
            w(
              "glp_transform_row: ind[" +
                h +
                "] = " +
                e +
                "; column index out of range"
            ),
          0 == d[h] &&
            w(
              "glp_transform_row: val[" +
                h +
                "] = 0; zero coefficient not allowed"
            ),
          0 != m[e] &&
            w(
              "glp_transform_row: ind[" +
                h +
                "] = " +
                e +
                "; duplicate column indices not allowed"
            ),
          (m[e] = d[h]);
      q = new Float64Array(1 + f);
      for (e = 1; e <= f; e++) (b = td(a, e)), (q[e] = b <= f ? 0 : m[b - f]);
      zd(a, q);
      b = 0;
      for (e = 1; e <= f; e++)
        zc(a, e) != A && ((p = -q[e]), 0 != p && (b++, (c[b] = e), (d[b] = p)));
      l = new Int32Array(1 + f);
      r = new Float64Array(1 + f);
      for (e = 1; e <= g; e++)
        if (Cc(a, e) != A) {
          p = m[e];
          k = gb(a, e, l, r);
          for (h = 1; h <= k; h++) p += r[h] * q[l[h]];
          0 != p && (b++, (c[b] = f + e), (d[b] = p));
        }
      return b;
    });
  exports.glp_transform_col = function (a, b, c, d) {
    var e, f, g, h;
    Ib(a) || w("glp_transform_col: basis factorization does not exist ");
    f = kb(a);
    h = new Float64Array(1 + f);
    (0 <= b && b <= f) ||
      w("glp_transform_col: len = " + b + "; invalid column length");
    for (g = 1; g <= b; g++)
      (e = c[g]),
        (1 <= e && e <= f) ||
          w(
            "glp_transform_col: ind[" +
              g +
              "] = " +
              e +
              "; row index out of range"
          ),
        0 == d[g] &&
          w(
            "glp_transform_col: val[" +
              g +
              "] = 0; zero coefficient not allowed"
          ),
        0 != h[e] &&
          w(
            "glp_transform_col: ind[" +
              g +
              "] = " +
              e +
              "; duplicate row indices not allowed"
          ),
        (h[e] = d[g]);
    xd(a, h);
    b = 0;
    for (e = 1; e <= f; e++)
      0 != h[e] && (b++, (c[b] = td(a, e)), (d[b] = h[e]));
    return b;
  };
  var Ed = (exports.glp_prim_rtest = function (a, b, c, d, e, f) {
      var g, h, k, l, p, m, q, r, n, t, y, E, C;
      tc(a) != dc &&
        w("glp_prim_rtest: basic solution is not primal feasible ");
      1 != e &&
        -1 != e &&
        w("glp_prim_rtest: dir = " + e + "; invalid parameter");
      (0 < f && 1 > f) ||
        w("glp_prim_rtest: eps = " + f + "; invalid parameter");
      h = kb(a);
      k = lb(a);
      l = 0;
      C = s;
      r = 0;
      for (p = 1; p <= b; p++)
        if (
          ((g = c[p]),
          (1 <= g && g <= h + k) ||
            w(
              "glp_prim_rtest: ind[" +
                p +
                "] = " +
                g +
                "; variable number out of range"
            ),
          g <= h
            ? ((m = ob(a, g)),
              (t = pb(a, g)),
              (y = qb(a, g)),
              (q = zc(a, g)),
              (n = Ac(a, g)))
            : ((m = rb(a, g - h)),
              (t = sb(a, g - h)),
              (y = tb(a, g - h)),
              (q = Cc(a, g - h)),
              (n = Dc(a, g - h))),
          q != A &&
            w(
              "glp_prim_rtest: ind[" +
                p +
                "] = " +
                g +
                "; non-basic variable not allowed"
            ),
          (g = 0 < e ? +d[p] : -d[p]),
          m != Ka)
        ) {
          if (m == Sa) {
            if (g > -f) continue;
            E = (t - n) / g;
          } else if (m == Ta) {
            if (g < +f) continue;
            E = (y - n) / g;
          } else if (m == I)
            if (0 > g) {
              if (g > -f) continue;
              E = (t - n) / g;
            } else {
              if (g < +f) continue;
              E = (y - n) / g;
            }
          else if (m == B) {
            if (-f < g && g < +f) continue;
            E = 0;
          }
          0 > E && (E = 0);
          if (C > E || (C == E && r < Math.abs(g)))
            (l = p), (C = E), (r = Math.abs(g));
        }
      return l;
    }),
    Fd = (exports.glp_dual_rtest = function (a, b, c, d, e, f) {
      var g, h, k, l, p, m, q, r, n, t, y;
      uc(a) != dc && w("glp_dual_rtest: basic solution is not dual feasible");
      1 != e &&
        -1 != e &&
        w("glp_dual_rtest: dir = " + e + "; invalid parameter");
      (0 < f && 1 > f) ||
        w("glp_dual_rtest: eps = " + f + "; invalid parameter");
      h = kb(a);
      k = lb(a);
      n = jb(a) == za ? 1 : -1;
      l = 0;
      y = s;
      q = 0;
      for (p = 1; p <= b; p++) {
        g = c[p];
        (1 <= g && g <= h + k) ||
          w(
            "glp_dual_rtest: ind[" +
              p +
              "] = " +
              g +
              "; variable number out of range"
          );
        g <= h
          ? ((m = zc(a, g)), (r = Bc(a, g)))
          : ((m = Cc(a, g - h)), (r = Ec(a, g - h)));
        m == A &&
          w(
            "glp_dual_rtest: ind[" +
              p +
              "] = " +
              g +
              "; basic variable not allowed"
          );
        g = 0 < e ? +d[p] : -d[p];
        if (m == G) {
          if (g < +f) continue;
          t = (n * r) / g;
        } else if (m == Ua) {
          if (g > -f) continue;
          t = (n * r) / g;
        } else if (m == Ra) {
          if (-f < g && g < +f) continue;
          t = 0;
        } else if (m == Na) continue;
        0 > t && (t = 0);
        if (y > t || (y == t && q < Math.abs(g)))
          (l = p), (y = t), (q = Math.abs(g));
      }
      return l;
    });
  function Gd(a, b, c, d, e, f, g) {
    var h,
      k,
      l,
      p = 0,
      m,
      q;
    a.na == Aa &&
      w("glp_analyze_row: primal basic solution components are undefined");
    a.sa != dc && w("glp_analyze_row: basic solution is not dual feasible");
    (0 <= b && b <= a.i) ||
      w("glp_analyze_row: len = " + b + "; invalid row length");
    q = 0;
    for (h = 1; h <= b; h++)
      (k = c[h]),
        (1 <= k && k <= a.g + a.i) ||
          w(
            "glp_analyze_row: ind[" +
              h +
              "] = " +
              k +
              "; row/column index out of range"
          ),
        k <= a.g
          ? (a.n[k].m == A &&
              w(
                "glp_analyze_row: ind[" +
                  h +
                  "] = " +
                  k +
                  "; basic auxiliary variable is not allowed"
              ),
            (m = a.n[k].r))
          : (a.f[k - a.g].m == A &&
              w(
                "glp_analyze_row: ind[" +
                  h +
                  "] = " +
                  k +
                  "; basic structural variable is not allowed"
              ),
            (m = a.f[k - a.g].r)),
        (q += d[h] * m);
    if (e == Sa) {
      if (q >= f) return 1;
      l = 1;
    } else if (e == Ta) {
      if (q <= f) return 1;
      l = -1;
    } else w("glp_analyze_row: type = " + e + "; invalid parameter");
    e = f - q;
    b = Fd(a, b, c, d, l, 1e-9);
    if (0 == b) return 2;
    k = c[b];
    m = k <= a.g ? a.n[k].r : a.f[k - a.g].r;
    c = e / d[b];
    g(b, m, c, q, e, k <= a.g ? a.n[k].J * c : a.f[k - a.g].J * c);
    return p;
  }
  exports.glp_analyze_bound = function (a, b, c) {
    var d, e, f, g, h, k, l, p, m, q, r, n, t, y;
    r = n = t = y = null;
    (null != a && 3621377730 == a.Dd) ||
      w("glp_analyze_bound: P = " + a + "; invalid problem object");
    e = a.g;
    f = a.i;
    (a.na == dc && a.sa == dc) ||
      w("glp_analyze_bound: optimal basic solution required");
    0 == e || a.valid || w("glp_analyze_bound: basis factorization required");
    (1 <= b && b <= e + f) ||
      w("glp_analyze_bound: k = " + b + "; variable number out of range");
    d = b <= e ? a.n[b] : a.f[b - e];
    g = d.m;
    f = d.r;
    g == A &&
      w("glp_analyze_bound: k = " + b + "; basic variable not allowed ");
    g = new Int32Array(1 + e);
    q = new Float64Array(1 + e);
    k = Cd(a, b, g, q);
    for (b = -1; 1 >= b; b += 2)
      (l = Ed(a, k, g, q, b, 1e-9)),
        0 == l
          ? ((h = 0), (l = 0 > b ? -s : +s))
          : ((h = g[l]),
            h <= e
              ? ((d = a.n[h]), (p = pb(a, d.ea)), (m = qb(a, d.ea)))
              : ((d = a.f[h - e]), (p = sb(a, d.C)), (m = tb(a, d.C))),
            (d = d.r),
            (d = (0 > b && 0 < q[l]) || (0 < b && 0 > q[l]) ? p - d : m - d),
            (l = f + d / q[l])),
        0 > b ? ((r = l), (n = h)) : ((t = l), (y = h));
    c(r, n, t, y);
  };
  exports.glp_analyze_coef = function (a, b, c) {
    var d,
      e,
      f,
      g,
      h,
      k,
      l,
      p,
      m,
      q,
      r,
      n,
      t,
      y,
      E,
      C,
      D,
      H,
      R,
      V,
      O = null,
      Q = null,
      F = null,
      W = null,
      X = null,
      ca = null;
    (null != a && 3621377730 == a.Dd) ||
      w("glp_analyze_coef: P = " + a + "; invalid problem object");
    e = a.g;
    f = a.i;
    (a.na == dc && a.sa == dc) ||
      w("glp_analyze_coef: optimal basic solution required");
    0 == e || a.valid || w("glp_analyze_coef: basis factorization required");
    (1 <= b && b <= e + f) ||
      w("glp_analyze_coef: k = " + b + "; variable number out of range");
    b <= e
      ? ((d = a.n[b]), (g = d.type), (n = d.c), (t = d.d), (y = 0))
      : ((d = a.f[b - e]), (g = d.type), (n = d.c), (t = d.d), (y = d.u));
    h = d.m;
    E = d.r;
    h != A &&
      w("glp_analyze_coef: k = " + b + "; non-basic variable not allowed");
    h = new Int32Array(1 + e);
    V = new Float64Array(1 + e);
    r = new Int32Array(1 + f);
    R = new Float64Array(1 + f);
    m = Bd(a, b, r, R);
    for (f = -1; 1 >= f; f += 2)
      a.dir == za ? (l = -f) : a.dir == Ea && (l = +f),
        (q = Fd(a, m, r, R, l, 1e-9)),
        0 == q
          ? ((C = 0 > f ? -s : +s), (k = 0), (q = E))
          : ((k = r[q]),
            (d = k <= e ? a.n[k] : a.f[k - e]),
            (l = d.J),
            (d = -l / R[q]),
            (C = y + d),
            (l = (0 > f && 0 < R[q]) || (0 < f && 0 > R[q]) ? 1 : -1),
            a.dir == Ea && (l = -l),
            (p = Cd(a, k, h, V)),
            (d = b <= e ? a.n[b] : a.f[b - e]),
            (d.type = Ka),
            (d.c = d.d = 0),
            (p = Ed(a, p, h, V, l, 1e-9)),
            (d = b <= e ? a.n[b] : a.f[b - e]),
            (d.type = g),
            (d.c = n),
            (d.d = t),
            0 == p
              ? (q = (0 > l && 0 < R[q]) || (0 < l && 0 > R[q]) ? -s : +s)
              : ((d = h[p]),
                d <= e
                  ? ((d = a.n[d]), (D = pb(a, d.ea)), (H = qb(a, d.ea)))
                  : ((d = a.f[d - e]), (D = sb(a, d.C)), (H = tb(a, d.C))),
                (d = d.r),
                (d =
                  (0 > l && 0 < V[p]) || (0 < l && 0 > V[p]) ? D - d : H - d),
                (q = E + (R[q] / V[p]) * d))),
        0 > f ? ((O = C), (Q = k), (F = q)) : ((W = C), (X = k), (ca = q));
    c(O, Q, F, W, X, ca);
  };
  exports.glp_ios_reason = function (a) {
    return a.reason;
  };
  exports.glp_ios_get_prob = function (a) {
    return a.A;
  };
  function Hd(a) {
    a.reason != Ia && w("glp_ios_pool_size: operation not allowed");
    return a.Bd.size;
  }
  function Id(a, b, c, d, e, f, g) {
    a.reason != Ia && w("glp_ios_add_row: operation not allowed");
    var h = a.Bd,
      k,
      l;
    k = { name: null };
    (0 <= b && 255 >= b) ||
      w("glp_ios_add_row: klass = " + b + "; invalid cut class");
    k.qc = b;
    k.k = null;
    (0 <= c && c <= a.i) ||
      w("glp_ios_add_row: len = " + c + "; invalid cut length");
    for (l = 1; l <= c; l++)
      (b = {}),
        (1 <= d[l] && d[l] <= a.i) ||
          w(
            "glp_ios_add_row: ind[" +
              l +
              "] = " +
              d[l] +
              "; column index out of range"
          ),
        (b.C = d[l]),
        (b.j = e[l]),
        (b.e = k.k),
        (k.k = b);
    f != Sa &&
      f != Ta &&
      f != B &&
      w("glp_ios_add_row: type = " + f + "; invalid cut type");
    k.type = f;
    k.cg = g;
    k.ca = h.Xa;
    k.e = null;
    null == k.ca ? (h.head = k) : (k.ca.e = k);
    h.Xa = k;
    h.size++;
  }
  function Jd(a, b) {
    (1 <= b && b <= a.A.i) ||
      w("glp_ios_can_branch: j = " + b + "; column number out of range");
    return a.ad[b];
  }
  function Kd(a, b) {
    var c = a.A,
      d = a.wc,
      e = a.i,
      f,
      g;
    g = c.ha;
    for (f = 1; f <= e; f++) {
      var h = c.f[f];
      if (h.kind == Fc && b[f] != Math.floor(b[f])) return 1;
      g += h.u * b[f];
    }
    if (c.za == dc)
      switch (c.dir) {
        case za:
          if (g >= a.A.ta) return 1;
          break;
        case Ea:
          if (g <= a.A.ta) return 1;
      }
    a.p.o >= fc && x("Solution found by heuristic: " + g + "");
    c.za = dc;
    c.ta = g;
    for (f = 1; f <= e; f++) c.f[f].Sa = b[f];
    for (e = 1; e <= d; e++)
      for (f = c.n[e], f.Sa = 0, g = f.k; null != g; g = g.B)
        f.Sa += g.j * g.f.Sa;
    return 0;
  }
  exports.glp_mpl_alloc_wksp = function () {
    return Ld();
  };
  exports._glp_mpl_init_rand = function (a, b) {
    0 != a.D && w("glp_mpl_init_rand: invalid call sequence\n");
    Md(a.Gd, b);
  };
  var Od = (exports.glp_mpl_read_model = function (a, b, c, d) {
    0 != a.D && w("glp_mpl_read_model: invalid call sequence");
    a = Nd(a, b, c, d);
    1 == a || 2 == a ? (a = 0) : 4 == a && (a = 1);
    return a;
  });
  exports.glp_mpl_read_model_from_string = function (a, b, c, d) {
    var e = 0;
    return Od(
      a,
      b,
      function () {
        return e < c.length ? c[e++] : -1;
      },
      d
    );
  };
  var Qd = (exports.glp_mpl_read_data = function (a, b, c) {
    1 != a.D && 2 != a.D && w("glp_mpl_read_data: invalid call sequence");
    a = Pd(a, b, c);
    2 == a ? (a = 0) : 4 == a && (a = 1);
    return a;
  });
  exports.glp_mpl_read_data_from_string = function (a, b, c) {
    var d = 0;
    return Qd(a, b, function () {
      return d < c.length ? c[d++] : -1;
    });
  };
  exports.glp_mpl_generate = function (a, b, c, d) {
    1 != a.D && 2 != a.D && w("glp_mpl_generate: invalid call sequence\n");
    a = Rd(a, b, c, d);
    3 == a ? (a = 0) : 4 == a && (a = 1);
    return a;
  };
  exports.glp_mpl_build_prob = function (a, b) {
    var c, d, e, f, g, h, k;
    3 != a.D && w("glp_mpl_build_prob: invalid call sequence\n");
    db(b);
    Ca(b, Sd(a));
    c = Td(a);
    0 < c && La(b, c);
    for (d = 1; d <= c; d++) {
      Pa(b, d, Ud(a, d));
      g = Vd(a, d, function (a, b) {
        h = a;
        k = b;
      });
      switch (g) {
        case Wd:
          g = Ka;
          break;
        case Xd:
          g = Sa;
          break;
        case Yd:
          g = Ta;
          break;
        case Zd:
          g = I;
          break;
        case $d:
          g = B;
      }
      g == I &&
        Math.abs(h - k) < 1e-9 * (1 + Math.abs(h)) &&
        ((g = B), Math.abs(h) <= Math.abs(k) ? (k = h) : (h = k));
      Va(b, d, g, h, k);
      0 != ae(a, d) &&
        x(
          "glp_mpl_build_prob: row " +
            Ud(a, d) +
            "; constant term " +
            ae(a, d) +
            " ignored"
        );
    }
    d = be(a);
    0 < d && Oa(b, d);
    for (e = 1; e <= d; e++) {
      Qa(b, e, ce(a, e));
      f = de(a, e);
      switch (f) {
        case ee:
        case fe:
          Hc(b, e, Fc);
      }
      g = ge(a, e, function (a, b) {
        h = a;
        k = b;
      });
      switch (g) {
        case Wd:
          g = Ka;
          break;
        case Xd:
          g = Sa;
          break;
        case Yd:
          g = Ta;
          break;
        case Zd:
          g = I;
          break;
        case $d:
          g = B;
      }
      if (f == fe) {
        if (g == Ka || g == Ta || 0 > h) h = 0;
        if (g == Ka || g == Sa || 1 < k) k = 1;
        g = I;
      }
      g == I &&
        Math.abs(h - k) < 1e-9 * (1 + Math.abs(h)) &&
        ((g = B), Math.abs(h) <= Math.abs(k) ? (k = h) : (h = k));
      Wa(b, e, g, h, k);
    }
    g = new Int32Array(1 + d);
    e = new Float64Array(1 + d);
    for (d = 1; d <= c; d++) (f = he(a, d, g, e)), Ya(b, d, f, g, e);
    for (d = 1; d <= c; d++)
      if (((f = ie(a, d)), f == je || f == ke)) {
        Da(b, Ud(a, d));
        Fa(b, f == je ? za : Ea);
        Xa(b, 0, ae(a, d));
        f = he(a, d, g, e);
        for (c = 1; c <= f; c++) Xa(b, g[c], e[c]);
        break;
      }
  };
  exports.glp_mpl_postsolve = function (a, b, c) {
    var d, e, f, g, h, k;
    (3 != a.D || a.Of) && w("glp_mpl_postsolve: invalid call sequence");
    c != Zb &&
      c != le &&
      c != Sc &&
      w("glp_mpl_postsolve: sol = " + c + "; invalid parameter");
    e = Td(a);
    f = be(a);
    (e == kb(b) && f == lb(b)) ||
      w("glp_mpl_postsolve: wrong problem object\n");
    if (!me(a)) return 0;
    for (d = 1; d <= e; d++)
      c == Zb
        ? ((g = zc(b, d)), (h = Ac(b, d)), (k = Bc(b, d)))
        : c == le
        ? ((g = 0), (h = glp_ipt_row_prim(b, d)), (k = glp_ipt_row_dual(b, d)))
        : c == Sc && ((g = 0), (h = kd(b, d)), (k = 0)),
        1e-9 > Math.abs(h) && (h = 0),
        1e-9 > Math.abs(k) && (k = 0),
        ne(a, d, g, h, k);
    for (d = 1; d <= f; d++)
      c == Zb
        ? ((g = Cc(b, d)), (h = Dc(b, d)), (k = Ec(b, d)))
        : c == le
        ? ((g = 0), (h = glp_ipt_col_prim(b, d)), (k = glp_ipt_col_dual(b, d)))
        : c == Sc && ((g = 0), (h = ld(b, d)), (k = 0)),
        1e-9 > Math.abs(h) && (h = 0),
        1e-9 > Math.abs(k) && (k = 0),
        oe(a, d, g, h, k);
    a = pe(a);
    3 == a ? (a = 0) : 4 == a && (a = 1);
    return a;
  };
  function qe(a, b) {
    var c, d, e;
    c = null;
    for (d = a.root; null != d; )
      (c = d),
        0 >= a.yg(a.info, b, c.key)
          ? ((e = 0), (d = c.left), c.pa++)
          : ((e = 1), (d = c.right));
    d = {};
    d.key = b;
    d.type = 0;
    d.link = null;
    d.pa = 1;
    d.R = c;
    d.ba = null == c ? 0 : e;
    d.wa = 0;
    d.left = null;
    d.right = null;
    a.size++;
    for (
      null == c ? (a.root = d) : 0 == e ? (c.left = d) : (c.right = d);
      null != c;

    ) {
      if (0 == e) {
        if (0 < c.wa) {
          c.wa = 0;
          break;
        }
        if (0 > c.wa) {
          re(a, c);
          break;
        }
        c.wa = -1;
      } else {
        if (0 > c.wa) {
          c.wa = 0;
          break;
        }
        if (0 < c.wa) {
          re(a, c);
          break;
        }
        c.wa = 1;
      }
      e = c.ba;
      c = c.R;
    }
    null == c && a.height++;
    return d;
  }
  function re(a, b) {
    var c, d, e, f, g;
    0 > b.wa
      ? ((c = b.R),
        (d = b.left),
        (e = d.right),
        0 >= d.wa
          ? (null == c
              ? (a.root = d)
              : 0 == b.ba
              ? (c.left = d)
              : (c.right = d),
            (b.pa -= d.pa),
            (d.R = c),
            (d.ba = b.ba),
            d.wa++,
            (d.right = b),
            (b.R = d),
            (b.ba = 1),
            (b.wa = -d.wa),
            (b.left = e),
            null != e && ((e.R = b), (e.ba = 0)))
          : ((f = e.left),
            (g = e.right),
            null == c ? (a.root = e) : 0 == b.ba ? (c.left = e) : (c.right = e),
            (b.pa -= d.pa + e.pa),
            (e.pa += d.pa),
            (b.wa = 0 <= e.wa ? 0 : 1),
            (d.wa = 0 >= e.wa ? 0 : -1),
            (e.R = c),
            (e.ba = b.ba),
            (e.wa = 0),
            (e.left = d),
            (e.right = b),
            (b.R = e),
            (b.ba = 1),
            (b.left = g),
            (d.R = e),
            (d.ba = 0),
            (d.right = f),
            null != f && ((f.R = d), (f.ba = 1)),
            null != g && ((g.R = b), (g.ba = 0))))
      : ((c = b.R),
        (d = b.right),
        (e = d.left),
        0 <= d.wa
          ? (null == c
              ? (a.root = d)
              : 0 == b.ba
              ? (c.left = d)
              : (c.right = d),
            (d.pa += b.pa),
            (d.R = c),
            (d.ba = b.ba),
            d.wa--,
            (d.left = b),
            (b.R = d),
            (b.ba = 0),
            (b.wa = -d.wa),
            (b.right = e),
            null != e && ((e.R = b), (e.ba = 1)))
          : ((f = e.left),
            (g = e.right),
            null == c ? (a.root = e) : 0 == b.ba ? (c.left = e) : (c.right = e),
            (d.pa -= e.pa),
            (e.pa += b.pa),
            (b.wa = 0 >= e.wa ? 0 : -1),
            (d.wa = 0 <= e.wa ? 0 : 1),
            (e.R = c),
            (e.ba = b.ba),
            (e.wa = 0),
            (e.left = b),
            (e.right = d),
            (b.R = e),
            (b.ba = 0),
            (b.right = f),
            (d.R = e),
            (d.ba = 1),
            (d.left = g),
            null != f && ((f.R = b), (f.ba = 1)),
            null != g && ((g.R = d), (g.ba = 0))));
  }
  var pd = 1,
    qd = 2,
    se = 3,
    te = 4,
    ue = 5;
  function od(a, b, c, d) {
    var e, f;
    f = a.valid = 0;
    switch (a.type) {
      case md:
        a.Za = null;
        null == a.gb &&
          ((f = {}),
          (f.hb = f.g = 0),
          (f.valid = 0),
          (f.ia = ve()),
          (f.$d = 50),
          (f.Pb = 0),
          (f.Ze = f.af = f.$e = null),
          (f.je = f.ie = null),
          (f.ug = null),
          (f.vg = null),
          (f.ic = 1e-6),
          (f.ag = 0),
          (a.gb = f),
          (f = 1));
        break;
      case rd:
      case sd:
        (a.gb = null),
          null == a.Za &&
            (we && x("lpf_create_it: warning: debug mode enabled"),
            (f = { valid: 0 }),
            (f.Nc = f.ef = 0),
            (f.ia = ve()),
            (f.g = 0),
            (f.Cf = null),
            (f.K = 50),
            (f.i = 0),
            (f.Md = f.Ld = null),
            (f.Od = f.Nd = null),
            (f.Qc = null),
            (f.Ge = f.Fe = null),
            (f.Ie = f.He = null),
            (f.Jd = 1e3),
            (f.od = 0),
            (f.Wb = null),
            (f.Xb = null),
            (f.jb = f.Gc = null),
            (a.Za = f),
            (f = 1));
    }
    null != a.gb ? (e = a.gb.ia) : null != a.Za && (e = a.Za.ia);
    f && (e.Va = a.Cd);
    e.cc = a.cc;
    e.xc = a.xc;
    e.gc = a.gc;
    e.Mb = a.Mb;
    e.sc = a.sc;
    null != a.gb && (f && (a.gb.$d = a.$c), (a.gb.ic = a.ic));
    null != a.Za && (f && (a.Za.K = a.vc), f && (a.Za.Jd = a.kd));
    if (null != a.gb) {
      a: {
        e = a.gb;
        1 > b && w("fhv_factorize: m = " + b + "; invalid parameter");
        1e8 < b && w("fhv_factorize: m = " + b + "; matrix too big");
        e.g = b;
        e.valid = 0;
        null == e.Ze && (e.Ze = new Int32Array(1 + e.$d));
        null == e.af && (e.af = new Int32Array(1 + e.$d));
        null == e.$e && (e.$e = new Int32Array(1 + e.$d));
        e.hb < b &&
          ((e.hb = b + 100),
          (e.je = new Int32Array(1 + e.hb)),
          (e.ie = new Int32Array(1 + e.hb)),
          (e.ug = new Int32Array(1 + e.hb)),
          (e.vg = new Float64Array(1 + e.hb)));
        switch (xe(e.ia, b, c, d)) {
          case ye:
            b = ze;
            break a;
          case Ae:
            b = Be;
            break a;
        }
        e.valid = 1;
        e.Pb = 0;
        ga(e.je, 1, e.ia.kb, 1, b);
        ga(e.ie, 1, e.ia.vb, 1, b);
        b = e.ag = 0;
      }
      switch (b) {
        case ze:
          return (a = pd);
        case Be:
          return (a = qd);
      }
    } else if (null != a.Za) {
      a: {
        e = a.Za;
        if (we) var g, h, k, l, p, m;
        1 > b && w("lpf_factorize: m = " + b + "; invalid parameter");
        1e8 < b && w("lpf_factorize: m = " + b + "; matrix too big");
        e.ef = e.g = b;
        e.valid = 0;
        null == e.Md && (e.Md = new Int32Array(1 + e.K));
        null == e.Ld && (e.Ld = new Int32Array(1 + e.K));
        null == e.Od && (e.Od = new Int32Array(1 + e.K));
        null == e.Nd && (e.Nd = new Int32Array(1 + e.K));
        null == e.Qc &&
          ((f = e.K),
          Ce && x("scf_create_it: warning: debug mode enabled"),
          (1 <= f && 32767 >= f) ||
            w("scf_create_it: n_max = " + f + "; invalid parameter"),
          (g = {}),
          (g.K = f),
          (g.i = 0),
          (g.Nb = new Float64Array(1 + f * f)),
          (g.v = new Float64Array(1 + (f * (f + 1)) / 2)),
          (g.s = new Int32Array(1 + f)),
          (g.fg = De),
          (g.pa = 0),
          (g.l = Ce ? new Float64Array(1 + f * f) : null),
          (g.ig = new Float64Array(1 + f)),
          (e.Qc = g));
        null == e.Wb && (e.Wb = new Int32Array(1 + e.Jd));
        null == e.Xb && (e.Xb = new Float64Array(1 + e.Jd));
        e.Nc < b &&
          ((e.Nc = b + 100),
          (e.Ge = new Int32Array(1 + e.Nc + e.K)),
          (e.Fe = new Int32Array(1 + e.Nc + e.K)),
          (e.Ie = new Int32Array(1 + e.Nc + e.K)),
          (e.He = new Int32Array(1 + e.Nc + e.K)),
          (e.jb = new Float64Array(1 + e.Nc + e.K)),
          (e.Gc = new Float64Array(1 + e.Nc + e.K)));
        switch (xe(e.ia, b, c, d)) {
          case ye:
            b = Ee;
            break a;
          case Ae:
            b = LPF_ECOND;
            break a;
        }
        e.valid = 1;
        if (we) {
          e.Cf = p = new Float64Array(1 + b * b);
          l = new Int32Array(1 + b);
          m = new Float64Array(1 + b);
          for (f = 1; f <= b * b; f++) p[f] = 0;
          for (h = 1; h <= b; h++)
            for (k = c(d, h, l, m), f = 1; f <= k; f++)
              (g = l[f]), (p[(g - 1) * b + h] = m[f]);
        }
        e.i = 0;
        c = e.Qc;
        c.i = c.pa = 0;
        for (f = 1; f <= b; f++)
          (e.Ge[f] = e.Fe[f] = f), (e.Ie[f] = e.He[f] = f);
        e.od = 1;
        b = 0;
      }
      switch (b) {
        case 0:
          switch (a.type) {
            case rd:
              a.Za.Qc.fg = De;
              break;
            case sd:
              a.Za.Qc.fg = Fe;
          }
          break;
        case Ee:
          return (a = pd);
        case LPF_ECOND:
          return (a = qd);
      }
    }
    a.valid = 1;
    return (a.hg = 0);
  }
  function wd(a, b) {
    if (null != a.gb) {
      var c = a.gb,
        d = c.ia.kb,
        e = c.ia.vb,
        f = c.je,
        g = c.ie;
      c.valid || w("fhv_ftran: the factorization is not valid");
      c.ia.kb = f;
      c.ia.vb = g;
      Ge(c.ia, 0, b);
      c.ia.kb = d;
      c.ia.vb = e;
      He(c, 0, b);
      Ie(c.ia, 0, b);
    } else if (null != a.Za) {
      var c = a.Za,
        d = c.ef,
        e = c.g,
        h = c.i,
        k = c.Fe,
        f = c.He,
        g = c.jb,
        l,
        p;
      if (we) var m;
      c.valid || w("lpf_ftran: the factorization is not valid");
      if (we) for (m = new Float64Array(1 + e), l = 1; l <= e; l++) m[l] = b[l];
      for (l = 1; l <= d + h; l++) g[l] = (p = k[l]) <= e ? b[p] : 0;
      Ge(c.ia, 0, g);
      Je(c, g, d, g);
      Ke(c.Qc, 0, g, d);
      h = c.i;
      k = c.Md;
      l = c.Ld;
      p = c.Wb;
      var q = c.Xb,
        r,
        n,
        t,
        y;
      for (r = 1; r <= h; r++)
        if (0 != g[r + d])
          for (y = -1 * g[r + d], n = k[r], t = n + l[r]; n < t; n++)
            g[p[n]] += y * q[n];
      Ie(c.ia, 0, g);
      for (l = 1; l <= e; l++) b[l] = g[f[l]];
      we && Le(c, 0, b, m);
    }
  }
  function yd(a, b) {
    if (null != a.gb) {
      var c = a.gb,
        d = c.ia.kb,
        e = c.ia.vb,
        f = c.je,
        g = c.ie;
      c.valid || w("fhv_btran: the factorization is not valid");
      Ie(c.ia, 1, b);
      He(c, 1, b);
      c.ia.kb = f;
      c.ia.vb = g;
      Ge(c.ia, 1, b);
      c.ia.kb = d;
      c.ia.vb = e;
    } else if (null != a.Za) {
      var c = a.Za,
        d = c.ef,
        e = c.g,
        h = c.i,
        f = c.Ge,
        k = c.Ie,
        g = c.jb,
        l,
        p;
      if (we) var m;
      c.valid || w("lpf_btran: the factorization is not valid");
      if (we) for (m = new Float64Array(1 + e), l = 1; l <= e; l++) m[l] = b[l];
      for (l = 1; l <= d + h; l++) g[l] = (p = k[l]) <= e ? b[p] : 0;
      Ie(c.ia, 1, g);
      Me(c, g, d, g);
      Ke(c.Qc, 1, g, d);
      h = c.i;
      k = c.Od;
      l = c.Nd;
      p = c.Wb;
      var q = c.Xb,
        r,
        n,
        t,
        y;
      for (r = 1; r <= h; r++)
        if (0 != g[r + d])
          for (y = -1 * g[r + d], n = k[r], t = n + l[r]; n < t; n++)
            g[p[n]] += y * q[n];
      Ge(c.ia, 1, g);
      for (l = 1; l <= e; l++) b[l] = g[f[l]];
      we && Le(c, 1, b, m);
    }
  }
  function Ne(a, b, c, d, e, f) {
    if (null != a.gb)
      switch (Oe(a.gb, b, c, d, e, f)) {
        case ze:
          return (a.valid = 0), (a = pd);
        case Pe:
          return (a.valid = 0), (a = se);
        case Qe:
          return (a.valid = 0), (a = te);
        case Re:
          return (a.valid = 0), (a = ue);
      }
    else if (null != a.Za) {
      a: {
        var g = a.Za,
          h = g.ef,
          k = g.g;
        if (we) var l = g.Cf;
        var p = g.i,
          m = g.Md,
          q = g.Ld,
          r = g.Od,
          n = g.Nd,
          t = g.Ge,
          y = g.Fe,
          E = g.Ie,
          C = g.He,
          D = g.od,
          H = g.Wb,
          R = g.Xb,
          V = g.Gc,
          O = g.jb,
          Q = g.Gc,
          F,
          W,
          X;
        g.valid || w("lpf_update_it: the factorization is not valid");
        (1 <= b && b <= k) ||
          w("lpf_update_it: j = " + b + "; column number out of range");
        if (p == g.K) (g.valid = 0), (b = LPF_ELIMIT);
        else {
          for (F = 1; F <= k; F++) V[F] = 0;
          for (X = 1; X <= c; X++)
            (F = d[e + X]),
              (1 <= F && F <= k) ||
                w(
                  "lpf_update_it: ind[" +
                    X +
                    "] = " +
                    F +
                    "; row number out of range"
                ),
              0 != V[F] &&
                w(
                  "lpf_update_it: ind[" +
                    X +
                    "] = " +
                    F +
                    "; duplicate row index not allowed"
                ),
              0 == f[X] &&
                w(
                  "lpf_update_it: val[" +
                    X +
                    "] = " +
                    f[X] +
                    "; zero element not allowed"
                ),
              (V[F] = f[X]);
          if (we) for (F = 1; F <= k; F++) l[(F - 1) * k + b] = V[F];
          for (F = 1; F <= h + p; F++) O[F] = (W = y[F]) <= k ? V[W] : 0;
          for (F = 1; F <= h + p; F++) Q[F] = 0;
          Q[C[b]] = 1;
          Ge(g.ia, 0, O);
          Ie(g.ia, 1, Q);
          if (g.Jd < D + h + h) {
            F = g.Jd;
            c = g.od - 1;
            d = g.Wb;
            for (e = g.Xb; F < D + h + h; ) F += F;
            g.Jd = F;
            g.Wb = new Int32Array(1 + F);
            g.Xb = new Float64Array(1 + F);
            ga(g.Wb, 1, d, 1, c);
            ga(g.Xb, 1, e, 1, c);
            H = g.Wb;
            R = g.Xb;
          }
          m[p + 1] = D;
          for (F = 1; F <= h; F++)
            0 != O[F] && ((H[D] = F), (R[D] = O[F]), D++);
          q[p + 1] = D - g.od;
          g.od = D;
          r[p + 1] = D;
          for (F = 1; F <= h; F++)
            0 != Q[F] && ((H[D] = F), (R[D] = Q[F]), D++);
          n[p + 1] = D - g.od;
          g.od = D;
          Je(g, O, 0, O);
          Me(g, Q, 0, Q);
          q = 0;
          for (F = 1; F <= h; F++) q -= Q[F] * O[F];
          m = g.Qc;
          D = q;
          F = m.K;
          q = m.i;
          c = m.Nb;
          d = m.v;
          e = m.s;
          if (Ce) var ca = m.l;
          n = m.ig;
          r = 0;
          if (q == F) r = Se;
          else {
            m.i = ++q;
            f = 1;
            for (k = (f - 1) * m.K + q; f < q; f++, k += F) c[k] = 0;
            k = 1;
            for (f = (q - 1) * m.K + k; k < q; k++, f++) c[f] = 0;
            for (f = c[(q - 1) * m.K + q] = 1; f < q; f++) {
              H = 0;
              k = 1;
              for (l = (f - 1) * m.K + 1; k < q; k++, l++) H += c[l] * O[k + h];
              d[Te(m, f, q)] = H;
            }
            for (k = 1; k < q; k++) n[k] = Q[e[k] + h];
            n[q] = D;
            e[q] = q;
            if (Ce) {
              f = 1;
              for (k = (f - 1) * m.K + q; f < q; f++, k += F) ca[k] = O[f + h];
              k = 1;
              for (f = (q - 1) * m.K + k; k < q; k++, f++) ca[f] = Q[k + h];
              ca[(q - 1) * m.K + q] = D;
            }
            for (O = 1; O < q && 0 == n[O]; O++);
            switch (m.fg) {
              case De:
                Q = m.i;
                ca = m.Nb;
                for (D = m.v; O < Q; O++) {
                  e = Te(m, O, O);
                  c = (O - 1) * m.K + 1;
                  f = (Q - 1) * m.K + 1;
                  if (Math.abs(D[e]) < Math.abs(n[O])) {
                    F = O;
                    for (d = e; F <= Q; F++, d++)
                      (l = D[d]), (D[d] = n[F]), (n[F] = l);
                    F = 1;
                    d = c;
                    for (k = f; F <= Q; F++, d++, k++)
                      (l = ca[d]), (ca[d] = ca[k]), (ca[k] = l);
                  }
                  Math.abs(D[e]) < Ue && (D[e] = n[O] = 0);
                  if (0 != n[O]) {
                    l = n[O] / D[e];
                    F = O + 1;
                    for (d = e + 1; F <= Q; F++, d++) n[F] -= l * D[d];
                    F = 1;
                    d = c;
                    for (k = f; F <= Q; F++, d++, k++) ca[k] -= l * ca[d];
                  }
                }
                Math.abs(n[Q]) < Ue && (n[Q] = 0);
                D[Te(m, Q, Q)] = n[Q];
                break;
              case Fe:
                Ve(m, O, n);
            }
            F = m.K;
            O = m.i;
            Q = m.v;
            D = 0;
            ca = 1;
            for (n = Te(m, ca, ca); ca <= O; ca++, n += F, F--)
              0 != Q[n] && D++;
            m.pa = D;
            m.pa != q && (r = We);
            Ce && Le(m, "scf_update_exp");
          }
          switch (r) {
            case We:
              g.valid = 0;
              b = Ee;
              break a;
          }
          t[h + p + 1] = y[h + p + 1] = h + p + 1;
          E[h + p + 1] = C[h + p + 1] = h + p + 1;
          F = C[b];
          W = C[h + p + 1];
          E[F] = h + p + 1;
          C[h + p + 1] = F;
          E[W] = b;
          C[b] = W;
          g.i++;
          b = 0;
        }
      }
      switch (b) {
        case Ee:
          return (a.valid = 0), (a = pd);
        case LPF_ELIMIT:
          return (a.valid = 0), (a = te);
      }
    }
    a.hg++;
    return 0;
  }
  var Xe = "!\"#$%&()/,.;?@_`'{}|~",
    Ye = (exports.glp_read_lp = function (a, b, c) {
      function d(a, b) {
        throw Error(a.count + ": " + b);
      }
      function e(a, b) {
        x(a.count + ": warning: " + b);
      }
      function f(a) {
        var b;
        "\n" == a.l && a.count++;
        b = a.Qg();
        0 > b
          ? "\n" == a.l
            ? (a.count--, (b = -1))
            : (e(a, "missing final end of line"), (b = "\n"))
          : "\n" != b &&
            (0 <= " \t\n\v\f\r".indexOf(b)
              ? (b = " ")
              : ta(b) && d(a, "invalid control character " + b.charCodeAt(0)));
        a.l = b;
      }
      function g(a) {
        a.h += a.l;
        f(a);
      }
      function h(a, b) {
        return a.toLowerCase() == b.toLowerCase() ? 1 : 0;
      }
      function k(a) {
        function b() {
          for (a.b = Q; va(a.l) || 0 <= Xe.indexOf(a.l); ) g(a);
          c &&
            (h(a.h, "minimize")
              ? (a.b = y)
              : h(a.h, "minimum")
              ? (a.b = y)
              : h(a.h, "min")
              ? (a.b = y)
              : h(a.h, "maximize")
              ? (a.b = E)
              : h(a.h, "maximum")
              ? (a.b = E)
              : h(a.h, "max")
              ? (a.b = E)
              : h(a.h, "subject")
              ? " " == a.l &&
                (f(a),
                "t" == a.l.toLowerCase() &&
                  ((a.b = C),
                  (a.h += " "),
                  g(a),
                  "o" != a.l.toLowerCase() &&
                    d(a, "keyword `subject to' incomplete"),
                  g(a),
                  ua(a.l) &&
                    d(a, "keyword `" + a.h + a.l + "...' not recognized")))
              : h(a.h, "such")
              ? " " == a.l &&
                (f(a),
                "t" == a.l.toLowerCase() &&
                  ((a.b = C),
                  (a.h += " "),
                  g(a),
                  "h" != a.l.toLowerCase() &&
                    d(a, "keyword `such that' incomplete"),
                  g(a),
                  "a" != a.l.toLowerCase() &&
                    d(a, "keyword `such that' incomplete"),
                  g(a),
                  "t" != a.l.toLowerCase() &&
                    d(a, "keyword `such that' incomplete"),
                  g(a),
                  ua(a.l) &&
                    d(a, "keyword `" + a.h + a.l + "...' not recognized")))
              : h(a.h, "st")
              ? (a.b = C)
              : h(a.h, "s.t.")
              ? (a.b = C)
              : h(a.h, "st.")
              ? (a.b = C)
              : h(a.h, "bounds")
              ? (a.b = D)
              : h(a.h, "bound")
              ? (a.b = D)
              : h(a.h, "general")
              ? (a.b = H)
              : h(a.h, "generals")
              ? (a.b = H)
              : h(a.h, "gen")
              ? (a.b = H)
              : h(a.h, "integer")
              ? (a.b = R)
              : h(a.h, "integers")
              ? (a.b = R)
              : h(a.h, "int")
              ? (a.b = R)
              : h(a.h, "binary")
              ? (a.b = V)
              : h(a.h, "binaries")
              ? (a.b = V)
              : h(a.h, "bin")
              ? (a.b = V)
              : h(a.h, "end") && (a.b = O));
        }
        var c;
        a.b = -1;
        a.h = "";
        for (a.value = 0; ; ) {
          for (c = 0; " " == a.l; ) f(a);
          if (-1 == a.l) a.b = t;
          else if ("\n" == a.l)
            if ((f(a), ua(a.l))) (c = 1), b();
            else continue;
          else if ("\\" == a.l) {
            for (; "\n" != a.l; ) f(a);
            continue;
          } else if (ua(a.l) || ("." != a.l && 0 <= Xe.indexOf(a.l))) b();
          else if (wa(a.l) || "." == a.l) {
            for (a.b = F; wa(a.l); ) g(a);
            if ("." == a.l)
              for (
                g(a),
                  1 != a.h.length ||
                    wa(a.l) ||
                    d(a, "invalid use of decimal point");
                wa(a.l);

              )
                g(a);
            if ("e" == a.l || "E" == a.l)
              for (
                g(a),
                  ("+" != a.l && "-" != a.l) || g(a),
                  wa(a.l) || d(a, "numeric constant `" + a.h + "' incomplete");
                wa(a.l);

              )
                g(a);
            a.value = Number(a.h);
            a.value == Number.NaN &&
              d(a, "numeric constant `" + a.h + "' out of range");
          } else
            "+" == a.l
              ? ((a.b = W), g(a))
              : "-" == a.l
              ? ((a.b = X), g(a))
              : ":" == a.l
              ? ((a.b = ca), g(a))
              : "<" == a.l
              ? ((a.b = ka), g(a), "=" == a.l && g(a))
              : ">" == a.l
              ? ((a.b = P), g(a), "=" == a.l && g(a))
              : "=" == a.l
              ? ((a.b = u),
                g(a),
                "<" == a.l
                  ? ((a.b = ka), g(a))
                  : ">" == a.l && ((a.b = P), g(a)))
              : d(a, "character `" + a.l + "' not recognized");
          break;
        }
        for (; " " == a.l; ) f(a);
      }
      function l(a, b) {
        var c = xb(a.Ka, b);
        if (0 == c) {
          c = Oa(a.Ka, 1);
          Qa(a.Ka, c, b);
          if (a.K < c) {
            var d = a.K,
              e = a.Z,
              f = a.j,
              g = a.ba,
              h = a.c,
              k = a.d;
            a.K += a.K;
            a.Z = new Int32Array(1 + a.K);
            ga(a.Z, 1, e, 1, d);
            a.j = new Float64Array(1 + a.K);
            ga(a.j, 1, f, 1, d);
            a.ba = new Int8Array(1 + a.K);
            ha(a.ba, 1, 0, a.K);
            ga(a.ba, 1, g, 1, d);
            a.c = new Float64Array(1 + a.K);
            ga(a.c, 1, h, 1, d);
            a.d = new Float64Array(1 + a.K);
            ga(a.d, 1, k, 1, d);
          }
          a.c[c] = +s;
          a.d[c] = -s;
        }
        return c;
      }
      function p(a) {
        for (var b, c = 0, e, f, g; ; )
          if (
            (a.b == W ? ((f = 1), k(a)) : a.b == X ? ((f = -1), k(a)) : (f = 1),
            a.b == F ? ((g = a.value), k(a)) : (g = 1),
            a.b != Q && d(a, "missing variable name"),
            (b = l(a, a.h)),
            a.ba[b] &&
              d(a, "multiple use of variable `" + a.h + "' not allowed"),
            c++,
            (a.Z[c] = b),
            (a.j[c] = f * g),
            (a.ba[b] = 1),
            k(a),
            a.b != W && a.b != X)
          ) {
            for (b = 1; b <= c; b++) a.ba[a.Z[b]] = 0;
            e = 0;
            for (b = 1; b <= c; b++)
              0 != a.j[b] && (e++, (a.Z[e] = a.Z[b]), (a.j[e] = a.j[b]));
            break;
          }
        return e;
      }
      function m(a, b, c) {
        a.c[b] != +s &&
          e(a, "lower bound of variable `" + nb(a.Ka, b) + "' redefined");
        a.c[b] = c;
      }
      function q(a, b, c) {
        a.d[b] != -s &&
          e(a, "upper bound of variable `" + nb(a.Ka, b) + "' redefined");
        a.d[b] = c;
      }
      function r(a) {
        var b, c, e, f;
        for (k(a); a.b == W || a.b == X || a.b == F || a.b == Q; )
          a.b == W || a.b == X
            ? ((c = 1),
              (f = a.b == W ? 1 : -1),
              k(a),
              a.b == F
                ? ((e = f * a.value), k(a))
                : h(a.h, "infinity") || h(a.h, "inf")
                ? (0 < f && d(a, "invalid use of `+inf' as lower bound"),
                  (e = -s),
                  k(a))
                : d(a, "missing lower bound"))
            : a.b == F
            ? ((c = 1), (e = a.value), k(a))
            : (c = 0),
            c &&
              (a.b != ka &&
                d(a, "missing `<', `<=', or `=<' after lower bound"),
              k(a)),
            a.b != Q && d(a, "missing variable name"),
            (b = l(a, a.h)),
            c && m(a, b, e),
            k(a),
            a.b == ka
              ? (k(a),
                a.b == W || a.b == X
                  ? ((f = a.b == W ? 1 : -1),
                    k(a),
                    a.b == F
                      ? (q(a, b, f * a.value), k(a))
                      : h(a.h, "infinity") || h(a.h, "inf")
                      ? (0 > f && d(a, "invalid use of `-inf' as upper bound"),
                        q(a, b, +s),
                        k(a))
                      : d(a, "missing upper bound"))
                  : a.b == F
                  ? (q(a, b, a.value), k(a))
                  : d(a, "missing upper bound"))
              : a.b == P
              ? (c && d(a, "invalid bound definition"),
                k(a),
                a.b == W || a.b == X
                  ? ((f = a.b == W ? 1 : -1),
                    k(a),
                    a.b == F
                      ? (m(a, b, f * a.value), k(a))
                      : h(a.h, "infinity") || 0 == h(a.h, "inf")
                      ? (0 < f && d(a, "invalid use of `+inf' as lower bound"),
                        m(a, b, -s),
                        k(a))
                      : d(a, "missing lower bound"))
                  : a.b == F
                  ? (m(a, b, a.value), k(a))
                  : d(a, "missing lower bound"))
              : a.b == u
              ? (c && d(a, "invalid bound definition"),
                k(a),
                a.b == W || a.b == X
                  ? ((f = a.b == W ? 1 : -1),
                    k(a),
                    a.b == F
                      ? (m(a, b, f * a.value), q(a, b, f * a.value), k(a))
                      : d(a, "missing fixed value"))
                  : a.b == F
                  ? (m(a, b, a.value), q(a, b, a.value), k(a))
                  : d(a, "missing fixed value"))
              : h(a.h, "free")
              ? (c && d(a, "invalid bound definition"),
                m(a, b, -s),
                q(a, b, +s),
                k(a))
              : c || d(a, "invalid bound definition");
      }
      function n(a) {
        var b, c;
        a.b == H
          ? ((c = 0), k(a))
          : a.b == R
          ? ((c = 0), k(a))
          : a.b == V && ((c = 1), k(a));
        for (; a.b == Q; )
          (b = l(a, a.h)), Hc(a.Ka, b, Fc), c && (m(a, b, 0), q(a, b, 1)), k(a);
      }
      var t = 0,
        y = 1,
        E = 2,
        C = 3,
        D = 4,
        H = 5,
        R = 6,
        V = 7,
        O = 8,
        Q = 9,
        F = 10,
        W = 11,
        X = 12,
        ca = 13,
        ka = 14,
        P = 15,
        u = 16,
        z = {};
      x("Reading problem data");
      null == b && (b = {});
      z.Ka = a;
      z.p = b;
      z.Qg = c;
      z.count = 0;
      z.l = "\n";
      z.b = t;
      z.h = "";
      z.value = 0;
      z.K = 100;
      z.Z = new Int32Array(1 + z.K);
      z.j = new Float64Array(1 + z.K);
      z.ba = new Int8Array(1 + z.K);
      ha(z.ba, 1, 0, z.K);
      z.c = new Float64Array(1 + z.K);
      z.d = new Float64Array(1 + z.K);
      db(a);
      vb(a);
      k(z);
      z.b != y && z.b != E && d(z, "`minimize' or `maximize' keyword missing");
      (function (a) {
        var b, c;
        a.b == y ? Fa(a.Ka, za) : a.b == E && Fa(a.Ka, Ea);
        k(a);
        a.b == Q && ":" == a.l ? (Da(a.Ka, a.h), k(a), k(a)) : Da(a.Ka, "obj");
        c = p(a);
        for (b = 1; b <= c; b++) Xa(a.Ka, a.Z[b], a.j[b]);
      })(z);
      z.b != C && d(z, "constraints section missing");
      (function (a) {
        var b, c, e;
        for (
          k(a);
          (b = La(a.Ka, 1)),
            a.b == Q && ":" == a.l
              ? (0 != wb(a.Ka, a.h) &&
                  d(a, "constraint `" + a.h + "' multiply defined"),
                Pa(a.Ka, b, a.h),
                k(a),
                k(a))
              : Pa(a.Ka, b, "r." + a.count),
            (c = p(a)),
            Ya(a.Ka, b, c, a.Z, a.j),
            a.b == ka
              ? ((e = Ta), k(a))
              : a.b == P
              ? ((e = Sa), k(a))
              : a.b == u
              ? ((e = B), k(a))
              : d(a, "missing constraint sense"),
            a.b == W ? ((c = 1), k(a)) : a.b == X ? ((c = -1), k(a)) : (c = 1),
            a.b != F && d(a, "missing right-hand side"),
            Va(a.Ka, b, e, c * a.value, c * a.value),
            "\n" != a.l &&
              -1 != a.l &&
              d(a, "invalid symbol(s) beyond right-hand side"),
            k(a),
            a.b == W || a.b == X || a.b == F || a.b == Q;

        );
      })(z);
      for (z.b == D && r(z); z.b == H || z.b == R || z.b == V; ) n(z);
      z.b == O
        ? k(z)
        : z.b == t
        ? e(z, "keyword `end' missing")
        : d(z, "symbol " + z.h + " in wrong position");
      z.b != t && d(z, "extra symbol(s) detected beyond `end'");
      var L, v;
      for (b = 1; b <= a.i; b++)
        (L = z.c[b]),
          (v = z.d[b]),
          L == +s && (L = 0),
          v == -s && (v = +s),
          (c =
            L == -s && v == +s
              ? Ka
              : v == +s
              ? Sa
              : L == -s
              ? Ta
              : L != v
              ? I
              : B),
          Wa(z.Ka, b, c, L, v);
      x(
        a.g +
          " row" +
          (1 == a.g ? "" : "s") +
          ", " +
          a.i +
          " column" +
          (1 == a.i ? "" : "s") +
          ", " +
          a.L +
          " non-zero" +
          (1 == a.L ? "" : "s")
      );
      0 < Jc(a) &&
        ((b = Jc(a)),
        (c = Kc(a)),
        1 == b
          ? 0 == c
            ? x("One variable is integer")
            : x("One variable is binary")
          : ((L = b + " integer variables, "),
            x(
              (0 == c
                ? L + "none"
                : 1 == c
                ? L + "one"
                : c == b
                ? L + "all"
                : L + c) +
                " of which " +
                (1 == c ? "is" : "are") +
                " binary"
            )));
      x(z.count + " lines were read");
      yb(a);
      $a(a);
      return 0;
    });
  exports.glp_write_lp = function (a, b, c) {
    function d(a) {
      if ("." == a[0] || wa(a[0])) return 1;
      for (var b = 0; b < a.length; b++)
        if (!va(a[b]) && 0 > Xe.indexOf(a[b])) return 1;
      return 0;
    }
    function e(a) {
      for (var b = 0; b < a.length; b++)
        " " == a[b]
          ? (a[b] = "_")
          : "-" == a[b]
          ? (a[b] = "~")
          : "[" == a[b]
          ? (a[b] = "(")
          : "]" == a[b] && (a[b] = ")");
    }
    function f(a, b) {
      var c;
      c = 0 == b ? ib(a.Ka) : mb(a.Ka, b);
      if (null == c) return 0 == b ? "obj" : "r_" + b;
      e(c);
      return d(c) ? (0 == b ? "obj" : "r_" + b) : c;
    }
    function g(a, b) {
      var c = nb(a.Ka, b);
      if (null == c) return "x_" + b;
      e(c);
      return d(c) ? "x_" + b : c;
    }
    function h() {
      c("End");
      r++;
      x(r + " lines were written");
      return 0;
    }
    var k = {},
      l,
      p,
      m,
      q,
      r,
      n;
    x("Writing problem data");
    null == b && (b = {});
    k.Ka = a;
    k.p = b;
    r = 0;
    c("\\* Problem: " + (null == a.name ? "Unknown" : a.name) + " *\\");
    r++;
    c("");
    r++;
    if (!(0 < a.g && 0 < a.i))
      return (
        x("Warning: problem has no rows/columns"),
        c("\\* WARNING: PROBLEM HAS NO ROWS/COLUMNS *\\"),
        r++,
        c(""),
        r++,
        h()
      );
    a.dir == za ? (c("Minimize"), r++) : a.dir == Ea && (c("Maximize"), r++);
    b = f(k, 0);
    n = " " + b + ":";
    p = 0;
    for (m = 1; m <= a.i; m++)
      if (((l = a.f[m]), 0 != l.u || null == l.k))
        p++,
          (b = g(k, m)),
          (q =
            0 == l.u
              ? " + 0 " + b
              : 1 == l.u
              ? " + " + b
              : -1 == l.u
              ? " - " + b
              : 0 < l.u
              ? " + " + l.u + " " + b
              : " - " + -l.u + " " + b),
          72 < n.length + q.length && (c(n), (n = ""), r++),
          (n += q);
    0 == p && ((q = " 0 " + g(k, 1)), (n += q));
    c(n);
    r++;
    0 != a.ha && (c("\\* constant term = " + a.ha + " *\\"), r++);
    c("");
    r++;
    c("Subject To");
    r++;
    for (m = 1; m <= a.g; m++)
      if (((l = a.n[m]), l.type != Ka)) {
        b = f(k, m);
        n = " " + b + ":";
        for (p = l.k; null != p; p = p.B)
          (b = g(k, p.f.C)),
            (q =
              1 == p.j
                ? " + " + b
                : -1 == p.j
                ? " - " + b
                : 0 < p.j
                ? " + " + p.j + " " + b
                : " - " + -p.j + " " + b),
            72 < n.length + q.length && (c(n), (n = ""), r++),
            (n += q);
        l.type == I
          ? ((q = " - ~r_" + m),
            72 < n.length + q.length && (c(n), (n = ""), r++),
            (n += q))
          : null == l.k && ((q = " 0 " + g(k, 1)), (n += q));
        if (l.type == Sa) q = " >= " + l.c;
        else if (l.type == Ta) q = " <= " + l.d;
        else if (l.type == I || l.type == B) q = " = " + l.c;
        72 < n.length + q.length && (c(n), (n = ""), r++);
        n += q;
        c(n);
        r++;
      }
    c("");
    r++;
    q = 0;
    for (m = 1; m <= a.g; m++)
      (l = a.n[m]),
        l.type == I &&
          (q || (c("Bounds"), (q = 1), r++),
          c(" 0 <= ~r_" + m + " <= " + (l.d - l.c)),
          r++);
    for (m = 1; m <= a.i; m++)
      if (((l = a.f[m]), l.type != Sa || 0 != l.c))
        q || (c("Bounds"), (q = 1), r++),
          (b = g(k, m)),
          l.type == Ka
            ? (c(" " + b + " free"), r++)
            : l.type == Sa
            ? (c(" " + b + " >= " + l.c), r++)
            : l.type == Ta
            ? (c(" -Inf <= " + b + " <= " + l.d), r++)
            : l.type == I
            ? (c(" " + l.c + " <= " + b + " <= " + l.d), r++)
            : l.type == B && (c(" " + b + " = " + l.c), r++);
    q && c("");
    r++;
    q = 0;
    for (m = 1; m <= a.i; m++)
      (l = a.f[m]),
        l.kind != Ma &&
          (q || (c("Generals"), (q = 1), r++), c(" " + g(k, m)), r++);
    q && (c(""), r++);
    return h();
  };
  exports.glp_read_lp_from_string = function (a, b, c) {
    var d = 0;
    return Ye(a, b, function () {
      return d < c.length ? c[d++] : -1;
    });
  };
  var ze = 1,
    Be = 2,
    Pe = 3,
    Qe = 4,
    Re = 5;
  function He(a, b, c) {
    var d = a.Pb,
      e = a.Ze,
      f = a.af,
      g = a.$e,
      h = a.ia.wb,
      k = a.ia.xb,
      l,
      p,
      m;
    a.valid || w("fhv_h_solve: the factorization is not valid");
    if (b)
      for (b = d; 1 <= b; b--) {
        if (((a = e[b]), (m = c[a]), 0 != m))
          for (l = f[b], p = l + g[b] - 1; l <= p; l++) c[h[l]] -= k[l] * m;
      }
    else
      for (b = 1; b <= d; b++) {
        a = e[b];
        m = c[a];
        l = f[b];
        for (p = l + g[b] - 1; l <= p; l++) m -= k[l] * c[h[l]];
        c[a] = m;
      }
  }
  function Oe(a, b, c, d, e, f) {
    var g = a.g,
      h = a.ia,
      k = h.Fc,
      l = h.Ec,
      p = h.pd,
      m = h.Af,
      q = h.Dc,
      r = h.Cc,
      n = h.Rc,
      t = h.kb,
      y = h.vb,
      E = h.sf,
      C = h.me,
      D = h.wb,
      H = h.xb,
      R = h.ze,
      V = h.Mb,
      O = a.Ze,
      Q = a.af,
      F = a.$e,
      W = a.je,
      X = a.ie,
      ca = a.ug,
      ka = a.vg,
      P = a.ic,
      u,
      z;
    a.valid || w("fhv_update_it: the factorization is not valid");
    (1 <= b && b <= g) ||
      w("fhv_update_it: j = " + b + "; column number out of range");
    if (a.Pb == a.$d) return (a.valid = 0), (a = Qe);
    for (u = 1; u <= g; u++) ka[u] = 0;
    for (z = 1; z <= c; z++)
      (u = d[e + z]),
        (1 <= u && u <= g) ||
          w(
            "fhv_update_it: ind[" + z + "] = " + u + "; row number out of range"
          ),
        0 != ka[u] &&
          w(
            "fhv_update_it: ind[" +
              z +
              "] = " +
              u +
              "; duplicate row index not allowed"
          ),
        0 == f[z] &&
          w(
            "fhv_update_it: val[" +
              z +
              "] = " +
              f[z] +
              "; zero element not allowed"
          ),
        (ka[u] = f[z]);
    a.ia.kb = W;
    a.ia.vb = X;
    Ge(a.ia, 0, ka);
    a.ia.kb = t;
    a.ia.vb = y;
    He(a, 0, ka);
    c = 0;
    for (u = 1; u <= g; u++)
      (e = ka[u]), 0 == e || Math.abs(e) < V || (c++, (ca[c] = u), (ka[c] = e));
    X = q[b];
    for (W = X + r[b] - 1; X <= W; X++) {
      u = D[X];
      f = k[u];
      for (z = f + l[u] - 1; D[f] != b; f++);
      D[f] = D[z];
      H[f] = H[z];
      l[u]--;
    }
    h.Qb -= r[b];
    r[b] = 0;
    e = E[b];
    d = 0;
    for (z = 1; z <= c; z++) {
      u = ca[z];
      if (l[u] + 1 > p[u] && Ze(h, u, l[u] + 10))
        return (a.valid = 0), (h.Va = h.Ga + h.Ga), (a = Re);
      f = k[u] + l[u];
      D[f] = b;
      H[f] = ka[z];
      l[u]++;
      d < y[u] && (d = y[u]);
    }
    if (n[b] < c && $e(h, b, c))
      return (a.valid = 0), (h.Va = h.Ga + h.Ga), (a = Re);
    X = q[b];
    ga(D, X, ca, 1, c);
    ga(H, X, ka, 1, c);
    r[b] = c;
    h.Qb += c;
    if (e > d) return (a.valid = 0), (a = ze);
    u = t[e];
    b = C[e];
    for (z = e; z < d; z++)
      (t[z] = t[z + 1]), (y[t[z]] = z), (C[z] = C[z + 1]), (E[C[z]] = z);
    t[d] = u;
    y[u] = d;
    C[d] = b;
    E[b] = d;
    for (b = 1; b <= g; b++) R[b] = 0;
    f = k[u];
    for (z = f + l[u] - 1; f <= z; f++) {
      b = D[f];
      R[b] = H[f];
      X = q[b];
      for (W = X + r[b] - 1; D[X] != u; X++);
      D[X] = D[W];
      H[X] = H[W];
      r[b]--;
    }
    h.Qb -= l[u];
    l[u] = 0;
    a.Pb++;
    O[a.Pb] = u;
    F[a.Pb] = 0;
    if (h.Ma - h.Fa < d - e && (af(h), h.Ma - h.Fa < d - e))
      return (a.valid = h.valid = 0), (h.Va = h.Ga + h.Ga), (a = Re);
    for (z = e; z < d; z++)
      if (((b = t[z]), (c = C[z]), 0 != R[c])) {
        y = R[c] / m[b];
        E = k[b];
        for (c = E + l[b] - 1; E <= c; E++) R[D[E]] -= y * H[E];
        h.Ma--;
        D[h.Ma] = b;
        H[h.Ma] = y;
        F[a.Pb]++;
      }
    0 == F[a.Pb] ? a.Pb-- : ((Q[a.Pb] = h.Ma), (a.ag += F[a.Pb]));
    m[u] = R[C[d]];
    c = 0;
    for (z = d + 1; z <= g; z++)
      if (((b = C[z]), (e = R[b]), !(Math.abs(e) < V))) {
        if (r[b] + 1 > n[b] && $e(h, b, r[b] + 10))
          return (a.valid = 0), (h.Va = h.Ga + h.Ga), (a = Re);
        X = q[b] + r[b];
        D[X] = u;
        H[X] = e;
        r[b]++;
        c++;
        ca[c] = b;
        ka[c] = e;
      }
    if (p[u] < c && Ze(h, u, c))
      return (a.valid = 0), (h.Va = h.Ga + h.Ga), (a = Re);
    f = k[u];
    ga(D, f, ca, 1, c);
    ga(H, f, ka, 1, c);
    l[u] = c;
    h.Qb += c;
    e = 0;
    u = t[d];
    f = k[u];
    for (z = f + l[u] - 1; f <= z; f++)
      e < Math.abs(H[f]) && (e = Math.abs(H[f]));
    b = C[d];
    X = q[b];
    for (W = X + r[b] - 1; X <= W; X++)
      e < Math.abs(H[X]) && (e = Math.abs(H[X]));
    return Math.abs(m[u]) < P * e ? ((a.valid = 0), (a = Pe)) : 0;
  }
  function ic(a) {
    function b(a, b, c, d, k, l) {
      var p,
        m,
        q,
        r,
        n,
        t,
        y,
        E,
        C,
        D,
        H,
        R,
        V,
        O,
        Q = 0;
      (0 < a && 0 < b) ||
        w("triang: m = " + a + "; n = " + b + "; invalid dimension");
      p = new Int32Array(1 + (a >= b ? a : b));
      m = new Int32Array(1 + a);
      q = new Int32Array(1 + b);
      r = new Int32Array(1 + a);
      n = new Int32Array(1 + a);
      y = new Int32Array(1 + b);
      E = new Int32Array(1 + b);
      for (D = 1; D <= b; D++) (H = d(c, -D, p)), (y[D] = m[H]), (m[H] = D);
      for (H = t = 0; H <= a; H++)
        for (D = m[H]; 0 != D; D = y[D]) (E[D] = t), (t = D);
      H = 0;
      for (D = t; 0 != D; D = E[D]) (y[D] = H), (H = D);
      for (C = 1; C <= a; C++)
        (m[C] = H = d(c, +C, p)),
          (r[C] = 0),
          (n[C] = q[H]),
          0 != n[C] && (r[n[C]] = C),
          (q[H] = C);
      for (C = 1; C <= a; C++) k[C] = 0;
      for (D = 1; D <= b; D++) l[D] = 0;
      R = 1;
      for (V = b; R <= V; ) {
        C = q[1];
        if (0 != C) {
          D = 0;
          for (O = d(c, +C, p); 1 <= O; O--) (H = p[O]), 0 == l[H] && (D = H);
          k[C] = l[D] = R;
          R++;
          Q++;
        } else (D = t), (l[D] = V), V--;
        0 == y[D] ? (t = E[D]) : (E[y[D]] = E[D]);
        0 != E[D] && (y[E[D]] = y[D]);
        for (O = d(c, -D, p); 1 <= O; O--)
          (C = p[O]),
            (H = m[C]),
            0 == r[C] ? (q[H] = n[C]) : (n[r[C]] = n[C]),
            0 != n[C] && (r[n[C]] = r[C]),
            (m[C] = --H),
            (r[C] = 0),
            (n[C] = q[H]),
            0 != n[C] && (r[n[C]] = C),
            (q[H] = C);
      }
      for (C = 1; C <= a; C++) 0 == k[C] && (k[C] = R++);
      for (D = 1; D <= b; D++);
      for (r = 1; r <= a; r++) m[r] = 0;
      for (C = 1; C <= a; C++) (r = k[C]), (m[r] = C);
      for (H = 1; H <= b; H++) q[H] = 0;
      for (D = 1; D <= b; D++) (H = l[D]), (q[H] = D);
      for (r = 1; r <= Q; r++) for (C = m[r], O = d(c, +C, p); 1 <= O; O--);
      return Q;
    }
    function c(a, b, c) {
      var d = kb(a);
      lb(a);
      var k,
        l,
        p,
        m = 0;
      if (0 < b) {
        k = +b;
        p = ub(a, k, c, null);
        for (b = 1; b <= p; b++)
          bf(a, c[b], function (a) {
            a != cf && (c[++m] = d + c[b]);
          });
        df(a, k, function (a) {
          a != cf && (c[++m] = k);
        });
      } else
        (l = -b),
          (p = function (b) {
            b != cf && (l <= d ? (c[++m] = l) : (m = gb(a, l - d, c, null)));
          }),
          l <= d ? df(a, l, p) : bf(a, l - d, p);
      return m;
    }
    function d(a) {
      var d = kb(a),
        g = lb(a),
        h,
        k,
        l,
        p,
        m,
        q,
        r,
        n = new Int32Array(1 + d + g);
      x("Constructing initial basis...");
      if (0 == d || 0 == g) Hb(a);
      else {
        m = new Int32Array(1 + d);
        k = new Int32Array(1 + d + g);
        p = b(d, d + g, a, c, m, k);
        3 <= ef(a) && x("Size of triangular part = " + p + "");
        q = new Int32Array(1 + d);
        r = new Int32Array(1 + d + g);
        for (h = 1; h <= d; h++) q[m[h]] = h;
        for (h = 1; h <= d + g; h++) r[k[h]] = h;
        for (l = 1; l <= d + g; l++) n[l] = -1;
        for (k = 1; k <= p; k++) (h = r[k]), (n[h] = ff);
        for (k = p + 1; k <= d; k++) (h = q[k]), (n[h] = ff);
        for (l = 1; l <= d + g; l++)
          n[l] != ff &&
            ((p = function (a, b, c) {
              switch (a) {
                case gf:
                  n[l] = hf;
                  break;
                case jf:
                  n[l] = kf;
                  break;
                case lf:
                  n[l] = mf;
                  break;
                case nf:
                  n[l] = Math.abs(b) <= Math.abs(c) ? kf : mf;
                  break;
                case cf:
                  n[l] = of;
              }
            }),
            l <= d ? df(a, l, p) : bf(a, l - d, p));
        for (l = 1; l <= d + g; l++)
          l <= d ? Fb(a, l, n[l] - ff + A) : Gb(a, l - d, n[l] - ff + A);
      }
    }
    0 == a.g || 0 == a.i ? Hb(a) : d(a);
  }
  function Mc(a, b) {
    var c, d;
    if (0 == a.Sc)
      for (
        c = a.fe,
          d = a.ya,
          a.fe = 0 == c ? 20 : c + c,
          a.ya = Array(1 + a.fe),
          ia(a.ya, 0, 1 + a.fe),
          null != d && ga(a.ya, 1, d, 1, c),
          d = a.fe;
        d > c;
        d--
      )
        (a.ya[d].rb = null), (a.ya[d].e = a.Sc), (a.Sc = d);
    d = a.Sc;
    a.Sc = a.ya[d].e;
    a.ya[d].e = 0;
    c = d;
    d = {};
    a.ya[c].rb = d;
    d.s = c;
    d.R = b;
    d.La = null == b ? 0 : b.La + 1;
    d.count = 0;
    d.Na = null;
    d.zc = null;
    d.ec = null;
    d.eg = 0;
    d.rc = null == b ? (a.A.dir == za ? -s : +s) : b.rc;
    d.bound = null == b ? (a.A.dir == za ? -s : +s) : b.bound;
    d.Tc = 0;
    d.tg = 0;
    d.Ag = 0;
    d.Zc = 0;
    d.Rd = 0;
    d.data = 0 == a.p.Me ? null : {};
    d.ja = null;
    d.ca = a.Xa;
    d.e = null;
    null == a.head ? (a.head = d) : (a.Xa.e = d);
    a.Xa = d;
    a.Pd++;
    a.Zf++;
    a.Jg++;
    null != b && b.count++;
    return d;
  }
  function pf(a, b) {
    var c = a.A,
      d,
      e,
      f,
      g;
    d = a.ya[b].rb;
    a.N = d;
    e = a.ya[1].rb;
    if (d != e) {
      for (d.ja = null; null != d; d = d.R) null != d.R && (d.R.ja = d);
      for (d = e; null != d; d = d.ja) {
        var h = c.g;
        e = c.i;
        if (null == d.ja) {
          a.Fg = h;
          a.Gg < h + e &&
            ((f = h + e + 100),
            (a.Gg = f),
            (a.qf = new Int8Array(1 + f)),
            (a.of = new Float64Array(1 + f)),
            (a.rf = new Float64Array(1 + f)),
            (a.pf = new Int8Array(1 + f)));
          for (f = 1; f <= h; f++)
            (g = c.n[f]),
              (a.qf[f] = g.type),
              (a.of[f] = g.c),
              (a.rf[f] = g.d),
              (a.pf[f] = g.m);
          for (f = 1; f <= e; f++)
            (g = c.f[f]),
              (a.qf[c.g + f] = g.type),
              (a.of[c.g + f] = g.c),
              (a.rf[c.g + f] = g.d),
              (a.pf[c.g + f] = g.m);
        }
        for (f = d.Na; null != f; f = f.e)
          f.pc <= h
            ? Va(c, f.pc, f.type, f.c, f.d)
            : Wa(c, f.pc - h, f.type, f.c, f.d);
        for (f = d.zc; null != f; f = f.e)
          f.pc <= h ? Fb(c, f.pc, f.m) : Gb(c, f.pc - h, f.m);
        if (null != d.ec) {
          var k,
            l,
            h = new Int32Array(1 + e);
          l = new Float64Array(1 + e);
          for (e = d.ec; null != e; e = e.e) {
            f = La(c, 1);
            Pa(c, f, e.name);
            c.n[f].La = d.La;
            c.n[f].origin = e.origin;
            c.n[f].qc = e.qc;
            Va(c, f, e.type, e.c, e.d);
            k = 0;
            for (g = e.k; null != g; g = g.e) k++, (h[k] = g.C), (l[k] = g.j);
            Ya(c, f, k, h, l);
            zb(c, f, e.ma);
            Fb(c, f, e.m);
          }
        }
      }
      for (d = a.N; null != d.Na; ) (f = d.Na), (d.Na = f.e);
      for (; null != d.zc; ) (f = d.zc), (d.zc = f.e);
      for (; null != d.ec; )
        for (e = d.ec, d.ec = e.e; null != e.k; ) (g = e.k), (e.k = g.e);
    }
  }
  function qf(a) {
    var b = a.A,
      c = b.g,
      d = b.i,
      e = a.N,
      f,
      g,
      h;
    if (null == e.R)
      for (
        a.Ig = c,
          a.se = new Int8Array(1 + c + d),
          a.qe = new Float64Array(1 + c + d),
          a.te = new Float64Array(1 + c + d),
          a.re = new Int8Array(1 + c + d),
          f = 1;
        f <= c + d;
        f++
      )
        (h = f <= c ? b.n[f] : b.f[f - c]),
          (a.se[f] = h.type),
          (a.qe[f] = h.c),
          (a.te[f] = h.d),
          (a.re[f] = h.m);
    else {
      var k = a.Ig,
        l = a.Fg;
      for (f = 1; f <= l + d; f++) {
        var p, m, q, r, n, t;
        p = a.qf[f];
        q = a.of[f];
        r = a.rf[f];
        g = a.pf[f];
        h = f <= l ? b.n[f] : b.f[f - l];
        m = h.type;
        n = h.c;
        t = h.d;
        h = h.m;
        if (p != m || q != n || r != t)
          (p = {}),
            (p.pc = f),
            (p.type = m),
            (p.c = n),
            (p.d = t),
            (p.e = e.Na),
            (e.Na = p);
        g != h && ((g = {}), (g.pc = f), (g.m = h), (g.e = e.zc), (e.zc = g));
      }
      if (l < c)
        for (
          m = new Int32Array(1 + d), n = new Float64Array(1 + d), g = c;
          g > l;
          g--
        ) {
          h = b.n[g];
          t = {};
          f = mb(b, g);
          t.name = null == f ? null : f;
          t.type = h.type;
          t.c = h.c;
          t.d = h.d;
          t.k = null;
          p = ub(b, g, m, n);
          for (f = 1; f <= p; f++)
            (q = {}), (q.C = m[f]), (q.j = n[f]), (q.e = t.k), (t.k = q);
          t.ma = h.ma;
          t.m = h.m;
          t.e = e.ec;
          e.ec = t;
        }
      if (c != k) {
        c -= k;
        e = new Int32Array(1 + c);
        for (g = 1; g <= c; g++) e[g] = k + g;
        ab(b, c, e);
      }
      c = b.g;
      for (g = 1; g <= c; g++)
        Va(b, g, a.se[g], a.qe[g], a.te[g]), Fb(b, g, a.re[g]);
      for (k = 1; k <= d; k++)
        Wa(b, k, a.se[c + k], a.qe[c + k], a.te[c + k]), Gb(b, k, a.re[c + k]);
    }
    a.N = null;
  }
  function rf(a, b, c) {
    var d;
    b = a.ya[b].rb;
    null == b.ca ? (a.head = b.e) : (b.ca.e = b.e);
    null == b.e ? (a.Xa = b.ca) : (b.e.ca = b.ca);
    b.ca = b.e = null;
    a.Pd--;
    for (d = 1; 2 >= d; d++) c[d] = Mc(a, b).s;
  }
  function sf(a, b) {
    var c;
    c = a.ya[b].rb;
    null == c.ca ? (a.head = c.e) : (c.ca.e = c.e);
    null == c.e ? (a.Xa = c.ca) : (c.e.ca = c.ca);
    c.ca = c.e = null;
    for (a.Pd--; ; ) {
      for (var d; null != c.Na; ) (d = c.Na), (c.Na = d.e);
      for (; null != c.zc; ) (d = c.zc), (c.zc = d.e);
      for (; null != c.ec; ) {
        d = c.ec;
        for (d.name = null; null != d.k; ) d.k = d.k.e;
        c.ec = d.e;
      }
      b = c.s;
      a.ya[b].rb = null;
      a.ya[b].e = a.Sc;
      a.Sc = b;
      c = c.R;
      a.Zf--;
      if (null != c && (c.count--, 0 == c.count)) continue;
      break;
    }
  }
  function tf(a, b, c) {
    var d = a.A,
      e = d.g,
      f,
      g,
      h,
      k,
      l = a.Bg,
      p = a.Tg,
      m,
      q;
    xc(d);
    Ib(d);
    a = d.f[b].r;
    b = Bd(d, e + b, l, p);
    for (f = -1; 1 >= f; f += 2)
      if (
        ((h = l),
        (g = Fd(d, b, h, p, f, 1e-9)),
        (g = 0 == g ? 0 : h[g]),
        0 == g)
      )
        d.dir == za
          ? 0 > f
            ? (m = +s)
            : (q = +s)
          : d.dir == Ea && (0 > f ? (m = -s) : (q = -s));
      else {
        for (h = 1; h <= b && l[h] != g; h++);
        h = p[h];
        g <= e
          ? ((k = d.n[g].m), (g = d.n[g].J))
          : ((k = d.f[g - e].m), (g = d.f[g - e].J));
        if (d.dir == za) {
          if ((k == G && 0 > g) || (k == Ua && 0 < g) || k == Ra) g = 0;
        } else
          d.dir == Ea &&
            ((k == G && 0 < g) || (k == Ua && 0 > g) || k == Ra) &&
            (g = 0);
        k = (0 > f ? Math.floor(a) : Math.ceil(a)) - a;
        k /= h;
        h = g * k;
        0 > f ? (m = d.aa + h) : (q = d.aa + h);
      }
    c(m, q);
  }
  function uf(a, b) {
    var c = a.A,
      d = c.i,
      e,
      f,
      g,
      h = a.Bg,
      k;
    g = 0;
    k = c.ha;
    e = 0;
    for (f = 1; f <= d; f++) {
      var l = c.f[f];
      if (0 != l.u)
        if (l.type == B) k += l.u * l.r;
        else {
          if (l.kind != Fc || l.u != Math.floor(l.u)) return b;
          2147483647 >= Math.abs(l.u) ? (h[++g] = Math.abs(l.u) | 0) : (e = 1);
        }
    }
    if (0 == e) {
      if (0 == g) return b;
      d = 0;
      for (e = 1; e <= g; e++) {
        if (1 == e) d = h[1];
        else for (f = h[e], l = void 0; 0 < f; ) (l = d % f), (d = f), (f = l);
        if (1 == d) break;
      }
      e = d;
    }
    c.dir == za
      ? b != +s &&
        ((c = (b - k) / e),
        c >= Math.floor(c) + 0.001 && ((c = Math.ceil(c)), (b = e * c + k)))
      : c.dir == Ea &&
        b != -s &&
        ((c = (b - k) / e),
        c <= Math.ceil(c) - 0.001 && ((c = Math.floor(c)), (b = e * c + k)));
    return b;
  }
  function vf(a, b) {
    var c = a.A,
      d = 1,
      e;
    if (c.za == dc)
      switch (((e = a.p.we * (1 + Math.abs(c.ta))), c.dir)) {
        case za:
          b >= c.ta - e && (d = 0);
          break;
        case Ea:
          b <= c.ta + e && (d = 0);
      }
    else
      switch (c.dir) {
        case za:
          b == +s && (d = 0);
          break;
        case Ea:
          b == -s && (d = 0);
      }
    return d;
  }
  function wf(a) {
    var b = null;
    switch (a.A.dir) {
      case za:
        for (a = a.head; null != a; a = a.e)
          if (null == b || b.bound > a.bound) b = a;
        break;
      case Ea:
        for (a = a.head; null != a; a = a.e)
          if (null == b || b.bound < a.bound) b = a;
    }
    return null == b ? 0 : b.s;
  }
  var xf = (exports.glp_ios_relative_gap = function (a) {
    var b = a.A,
      c;
    b.za == dc
      ? ((b = b.ta),
        (c = wf(a)),
        0 == c
          ? (a = 0)
          : ((a = a.ya[c].rb.bound),
            (a = Math.abs(b - a) / (Math.abs(b) + 2.220446049250313e-16))))
      : (a = s);
    return a;
  });
  function yf(a) {
    var b = a.A,
      c = new kc();
    switch (a.p.o) {
      case lc:
        c.o = lc;
        break;
      case Lb:
        c.o = Lb;
        break;
      case fc:
      case Wb:
        c.o = fc;
        break;
      case mc:
        c.o = Wb;
    }
    c.cb = Qb;
    c.fb = a.p.o < mc ? a.p.fb : 0;
    if (b.za == dc)
      switch (a.A.dir) {
        case za:
          c.jf = b.ta;
          break;
        case Ea:
          c.hf = b.ta;
      }
    b = sc(b, c);
    a.N.eg++;
    return b;
  }
  function Oc(a) {
    for (; null != a.head; ) {
      var b = a.head;
      for (a.head = b.e; null != b.k; ) b.k = b.k.e;
    }
    a.size = 0;
    a.head = a.Xa = null;
    a.fh = 0;
    a.N = null;
  }
  function zf(a, b) {
    function c(a, b, c, d, e) {
      var f, g, h;
      g = h = 0;
      for (f = 1; f <= a; f++)
        if (0 < b[f])
          if (c[f] == -s)
            if (0 == g) g = f;
            else {
              h = -s;
              g = 0;
              break;
            }
          else h += b[f] * c[f];
        else if (0 > b[f])
          if (d[f] == +s)
            if (0 == g) g = f;
            else {
              h = -s;
              g = 0;
              break;
            }
          else h += b[f] * d[f];
      e.Wd = h;
      e.Uf = g;
      g = h = 0;
      for (f = 1; f <= a; f++)
        if (0 < b[f])
          if (d[f] == +s)
            if (0 == g) g = f;
            else {
              h = +s;
              g = 0;
              break;
            }
          else h += b[f] * d[f];
        else if (0 > b[f])
          if (c[f] == -s)
            if (0 == g) g = f;
            else {
              h = +s;
              g = 0;
              break;
            }
          else h += b[f] * c[f];
      e.Vd = h;
      e.Tf = g;
    }
    function d(a, b) {
      b(0 == a.Uf ? a.Wd : -s, 0 == a.Tf ? a.Vd : +s);
    }
    function e(a, b, c, d, e, f, g, h) {
      var k, l, m, n;
      c == -s || a.Vd == +s
        ? (k = -s)
        : 0 == a.Tf
        ? 0 < b[g]
          ? (k = c - (a.Vd - b[g] * f[g]))
          : 0 > b[g] && (k = c - (a.Vd - b[g] * e[g]))
        : (k = a.Tf == g ? c - a.Vd : -s);
      d == +s || a.Wd == -s
        ? (l = +s)
        : 0 == a.Uf
        ? 0 < b[g]
          ? (l = d - (a.Wd - b[g] * e[g]))
          : 0 > b[g] && (l = d - (a.Wd - b[g] * f[g]))
        : (l = a.Uf == g ? d - a.Wd : +s);
      1e-6 > Math.abs(b[g])
        ? ((m = -s), (n = +s))
        : 0 < b[g]
        ? ((m = k == -s ? -s : k / b[g]), (n = l == +s ? +s : l / b[g]))
        : 0 > b[g] &&
          ((m = l == +s ? -s : l / b[g]), (n = k == -s ? +s : k / b[g]));
      h(m, n);
    }
    function f(a, b, c, e, f) {
      var g = 0,
        h = b[c],
        k = e[f],
        l = null,
        m = null;
      d(a, function (a, b) {
        l = a;
        m = b;
      });
      if (
        (h != -s && ((a = 0.001 * (1 + Math.abs(h))), m < h - a)) ||
        (k != +s && ((a = 0.001 * (1 + Math.abs(k))), l > k + a))
      )
        return 1;
      h != -s && ((a = 1e-12 * (1 + Math.abs(h))), l > h - a && (b[c] = -s));
      k != +s && ((a = 1e-12 * (1 + Math.abs(k))), m < k + a && (e[f] = +s));
      return g;
    }
    function g(a, b, c, d, f, g, h, k, l) {
      var m = 0,
        n,
        q,
        p = null,
        r = null;
      n = f[k];
      q = g[k];
      e(a, b, c, d, f, g, k, function (a, b) {
        p = a;
        r = b;
      });
      h &&
        (p != -s &&
          (p = 0.001 > p - Math.floor(p) ? Math.floor(p) : Math.ceil(p)),
        r != +s &&
          (r = 0.001 > Math.ceil(r) - r ? Math.ceil(r) : Math.floor(r)));
      if (
        (n != -s && ((a = 0.001 * (1 + Math.abs(n))), r < n - a)) ||
        (q != +s && ((a = 0.001 * (1 + Math.abs(q))), p > q + a))
      )
        return 1;
      p != -s && ((a = 0.001 * (1 + Math.abs(p))), n < p - a && (n = p));
      r != +s && ((a = 0.001 * (1 + Math.abs(r))), q > r + a && (q = r));
      n != -s &&
        q != +s &&
        ((a = Math.abs(n)),
        (b = Math.abs(q)),
        n > q - 1e-10 * (1 + (a <= b ? a : b)) &&
          (n == f[k]
            ? (q = n)
            : q == g[k]
            ? (n = q)
            : a <= b
            ? (q = n)
            : (n = q)));
      l(n, q);
      return m;
    }
    function h(a, b, c, d, e) {
      var f,
        g = 0;
      b < d &&
        (a || b == -s
          ? g++
          : ((f = c == +s ? 1 + Math.abs(b) : 1 + (c - b)),
            d - b >= 0.25 * f && g++));
      c > e &&
        (a || c == +s
          ? g++
          : ((f = b == -s ? 1 + Math.abs(c) : 1 + (c - b)),
            c - e >= 0.25 * f && g++));
      return g;
    }
    var k = a.A,
      l = k.g,
      p = k.i,
      m,
      q,
      r,
      n = 0,
      t,
      y,
      E,
      C;
    t = new Float64Array(1 + l);
    y = new Float64Array(1 + l);
    switch (k.za) {
      case Aa:
        t[0] = -s;
        y[0] = +s;
        break;
      case dc:
        switch (k.dir) {
          case za:
            t[0] = -s;
            y[0] = k.ta - k.ha;
            break;
          case Ea:
            (t[0] = k.ta - k.ha), (y[0] = +s);
        }
    }
    for (m = 1; m <= l; m++) (t[m] = pb(k, m)), (y[m] = qb(k, m));
    E = new Float64Array(1 + p);
    C = new Float64Array(1 + p);
    for (m = 1; m <= p; m++) (E[m] = sb(k, m)), (C[m] = tb(k, m));
    q = l + 1;
    r = new Int32Array(1 + q);
    for (m = 1; m <= q; m++) r[m] = m - 1;
    if (
      (function (a, b, d, e, k, l, m, n) {
        var q = a.g,
          p = a.i,
          r = {},
          t,
          u,
          E = 0,
          C,
          v,
          y,
          M,
          ba,
          J,
          ea,
          fa;
        C = new Int32Array(1 + p);
        v = new Int32Array(1 + q + 1);
        y = new Int32Array(1 + q + 1);
        M = new Int32Array(1 + q + 1);
        ba = new Float64Array(1 + p);
        J = new Float64Array(1 + p);
        ea = new Float64Array(1 + p);
        u = 0;
        for (t = 1; t <= l; t++) (q = m[t]), (v[++u] = q), (y[q] = 1);
        for (; 0 < u; )
          if (((q = v[u--]), (y[q] = 0), M[q]++, b[q] != -s || d[q] != +s)) {
            l = 0;
            if (0 == q)
              for (m = 1; m <= p; m++)
                (fa = a.f[m]), 0 != fa.u && (l++, (C[l] = m), (ba[l] = fa.u));
            else
              for (m = a.n[q].k; null != m; m = m.B)
                l++, (C[l] = m.f.C), (ba[l] = m.j);
            for (t = 1; t <= l; t++) (m = C[t]), (J[t] = e[m]), (ea[t] = k[m]);
            c(l, ba, J, ea, r);
            if (f(r, b, q, d, q)) {
              E = 1;
              break;
            }
            if (b[q] != -s || d[q] != +s)
              for (t = 1; t <= l; t++) {
                var sa,
                  K = null,
                  oa = null;
                m = C[t];
                fa = a.f[m];
                sa = fa.kind != Ma;
                if (
                  g(r, ba, b[q], d[q], J, ea, sa, t, function (a, b) {
                    K = a;
                    oa = b;
                  })
                )
                  return (E = 1);
                sa = h(sa, e[m], k[m], K, oa);
                e[m] = K;
                k[m] = oa;
                if (0 < sa)
                  for (m = fa.k; null != m; m = m.I)
                    (fa = m.n.ea),
                      M[fa] >= n ||
                        (b[fa] == -s && d[fa] == +s) ||
                        0 != y[fa] ||
                        ((v[++u] = fa), (y[fa] = 1));
              }
          }
        return E;
      })(k, t, y, E, C, q, r, b)
    )
      return 1;
    for (m = 1; m <= l; m++)
      zc(k, m) == A &&
        (t[m] == -s && y[m] == +s
          ? Va(k, m, Ka, 0, 0)
          : y[m] == +s
          ? Va(k, m, Sa, t[m], 0)
          : t[m] == -s && Va(k, m, Ta, 0, y[m]));
    for (m = 1; m <= p; m++)
      Wa(
        k,
        m,
        E[m] == -s && C[m] == +s
          ? Ka
          : C[m] == +s
          ? Sa
          : E[m] == -s
          ? Ta
          : E[m] != C[m]
          ? I
          : B,
        E[m],
        C[m]
      );
    return n;
  }
  function Nc(a) {
    function b(a, b) {
      var c, d, e, f;
      d = a.A.za == dc ? String(a.A.ta) : "not found yet";
      c = wf(a);
      0 == c
        ? (e = "tree is empty")
        : ((c = a.ya[c].rb.bound),
          (e = c == -s ? "-inf" : c == +s ? "+inf" : c));
      a.A.dir == za ? (f = ">=") : a.A.dir == Ea && (f = "<=");
      c = xf(a);
      x(
        "+" +
          a.A.$ +
          ": " +
          (b ? ">>>>>" : "mip =") +
          " " +
          d +
          " " +
          f +
          " " +
          e +
          " " +
          (0 == c
            ? "  0.0%"
            : 0.001 > c
            ? " < 0.1%"
            : 9.999 >= c
            ? "  " + Number(100 * c).toFixed(1) + "%"
            : "") +
          " (" +
          a.Pd +
          "; " +
          (a.Jg - a.Zf) +
          ")"
      );
      a.Lg = ja();
    }
    function c(a, b) {
      return vf(a, a.ya[b].rb.bound);
    }
    function d(a) {
      var b = a.A,
        c,
        d,
        e = 0,
        f,
        g,
        h,
        k,
        l,
        m = 0;
      for (c = 1; c <= b.i; c++)
        if (((h = b.f[c]), (a.ad[c] = 0), h.kind == Fc && h.m == A)) {
          d = h.type;
          f = h.c;
          g = h.d;
          h = h.r;
          if (d == Sa || d == I || d == B) {
            k = f - a.p.Ub;
            l = f + a.p.Ub;
            if (k <= h && h <= l) continue;
            if (h < f) continue;
          }
          if (d == Ta || d == I || d == B) {
            k = g - a.p.Ub;
            l = g + a.p.Ub;
            if (k <= h && h <= l) continue;
            if (h > g) continue;
          }
          k = Math.floor(h + 0.5) - a.p.Ub;
          l = Math.floor(h + 0.5) + a.p.Ub;
          (k <= h && h <= l) ||
            ((a.ad[c] = 1),
            e++,
            (k = h - Math.floor(h)),
            (l = Math.ceil(h) - h),
            (m += k <= l ? k : l));
        }
      a.N.Ag = e;
      a.N.Zc = m;
      a.p.o >= mc &&
        (0 == e
          ? x("There are no fractional columns")
          : 1 == e
          ? x(
              "There is one fractional column, integer infeasibility is " +
                m +
                ""
            )
          : x(
              "There are " +
                e +
                " fractional columns, integer infeasibility is " +
                m +
                ""
            ));
    }
    function e(a) {
      var b = a.A,
        c;
      b.za = dc;
      b.ta = b.aa;
      for (c = 1; c <= b.g; c++) {
        var d = b.n[c];
        d.Sa = d.r;
      }
      for (c = 1; c <= b.i; c++)
        (d = b.f[c]),
          d.kind == Ma
            ? (d.Sa = d.r)
            : d.kind == Fc && (d.Sa = Math.floor(d.r + 0.5));
      a.qh++;
    }
    function f(a, b, c) {
      var d = a.A,
        e,
        f = d.g,
        g,
        h,
        k,
        l,
        m,
        n = Array(3),
        q,
        p,
        r,
        t,
        y = null,
        v = null,
        S;
      g = d.f[b].type;
      q = d.f[b].c;
      p = d.f[b].d;
      e = d.f[b].r;
      r = Math.floor(e);
      t = Math.ceil(e);
      switch (g) {
        case Ka:
          h = Ta;
          k = Sa;
          break;
        case Sa:
          h = q == r ? B : I;
          k = Sa;
          break;
        case Ta:
          h = Ta;
          k = t == p ? B : I;
          break;
        case I:
          (h = q == r ? B : I), (k = t == p ? B : I);
      }
      tf(a, b, function (a, b) {
        y = a;
        v = b;
      });
      g = uf(a, y);
      S = uf(a, v);
      l = !vf(a, g);
      m = !vf(a, S);
      if (l && m)
        return a.p.o >= mc && x("Both down- and up-branches are hopeless"), 2;
      if (m)
        return (
          a.p.o >= mc && x("Up-branch is hopeless"),
          Wa(d, b, h, q, r),
          (a.N.rc = y),
          d.dir == za
            ? a.N.bound < g && (a.N.bound = g)
            : d.dir == Ea && a.N.bound > g && (a.N.bound = g),
          1
        );
      if (l)
        return (
          a.p.o >= mc && x("Down-branch is hopeless"),
          Wa(d, b, k, t, p),
          (a.N.rc = v),
          d.dir == za
            ? a.N.bound < S && (a.N.bound = S)
            : d.dir == Ea && a.N.bound > S && (a.N.bound = S),
          1
        );
      a.p.o >= mc &&
        x("Branching on column " + b + ", primal value is " + e + "");
      l = a.N.s;
      a.N.Tc = b;
      a.N.tg = e;
      qf(a);
      rf(a, l, n);
      a.p.o >= mc &&
        x(
          "Node " +
            n[1] +
            " begins down branch, node " +
            n[2] +
            " begins up branch "
        );
      e = a.ya[n[1]].rb;
      e.Na = {};
      e.Na.pc = f + b;
      e.Na.type = h;
      e.Na.c = q;
      e.Na.d = r;
      e.Na.e = null;
      e.rc = y;
      d.dir == za
        ? e.bound < g && (e.bound = g)
        : d.dir == Ea && e.bound > g && (e.bound = g);
      e = a.ya[n[2]].rb;
      e.Na = {};
      e.Na.pc = f + b;
      e.Na.type = k;
      e.Na.c = t;
      e.Na.d = p;
      e.Na.e = null;
      e.rc = v;
      d.dir == za
        ? e.bound < S && (e.bound = S)
        : d.dir == Ea && e.bound > S && (e.bound = S);
      c == Af ? (a.ud = 0) : c == Bf ? (a.ud = n[1]) : c == Cf && (a.ud = n[2]);
      return 0;
    }
    function g(a) {
      var b = a.A,
        c,
        d,
        e = 0,
        f,
        g,
        h,
        k;
      f = b.aa;
      for (c = 1; c <= b.i; c++)
        if (((k = b.f[c]), k.kind == Fc))
          switch (((g = k.c), (h = k.d), (d = k.m), (k = k.J), b.dir)) {
            case za:
              d == G
                ? (0 > k && (k = 0), f + k >= b.ta && (Wa(b, c, B, g, g), e++))
                : d == Ua &&
                  (0 < k && (k = 0), f - k >= b.ta && (Wa(b, c, B, h, h), e++));
              break;
            case Ea:
              d == G
                ? (0 < k && (k = 0), f + k <= b.ta && (Wa(b, c, B, g, g), e++))
                : d == Ua &&
                  (0 > k && (k = 0), f - k <= b.ta && (Wa(b, c, B, h, h), e++));
          }
      a.p.o >= mc &&
        0 != e &&
        (1 == e
          ? x("One column has been fixed by reduced cost")
          : x(e + " columns have been fixed by reduced costs"));
    }
    function h(a) {
      var b,
        c = 0,
        d = null;
      for (b = a.wc + 1; b <= a.A.g; b++)
        a.A.n[b].origin == Ja &&
          a.A.n[b].La == a.N.La &&
          a.A.n[b].m == A &&
          (null == d && (d = new Int32Array(1 + a.A.g)), (d[++c] = b));
      0 < c && (ab(a.A, c, d), Jb(a.A));
    }
    function k(a) {
      var b = a.A,
        c,
        d = 0,
        e = 0,
        f = 0,
        g = 0,
        h = 0;
      for (c = b.g; 0 < c; c--) {
        var k = b.n[c];
        k.origin == Ja &&
          (k.qc == Df
            ? d++
            : k.qc == Ef
            ? e++
            : k.qc == Ff
            ? f++
            : k.qc == Gf
            ? g++
            : h++);
      }
      0 < d + e + f + g + h &&
        (x("Cuts on level " + a.N.La + ":"),
        0 < d && x(" gmi = " + d + ";"),
        0 < e && x(" mir = " + e + ";"),
        0 < f && x(" cov = " + f + ";"),
        0 < g && x(" clq = " + g + ";"),
        0 < h && x(" app = " + h + ";"),
        x(""));
    }
    function l(a) {
      if (a.p.Ed == bb || a.p.Ad == bb || a.p.xd == bb || a.p.vd == bb) {
        var b, c, d;
        c = a.i;
        1e3 > c && (c = 1e3);
        d = 0;
        for (b = a.wc + 1; b <= a.A.g; b++) a.A.n[b].origin == Ja && d++;
        if (!(d >= c)) {
          a.p.Ad == bb && 5 > a.N.Rd && Hf(a);
          a.p.Ed == bb && If(a, a.Xf);
          if (a.p.xd == bb) {
            b = a.A;
            c = kb(b);
            var e = lb(b),
              f,
              g,
              h,
              k,
              l;
            xc(b);
            d = new Int32Array(1 + e);
            k = new Float64Array(1 + e);
            l = new Float64Array(1 + e);
            for (e = 1; e <= c; e++)
              for (h = 1; 2 >= h; h++) {
                g = ob(b, e) - Ka + gf;
                if (1 == h) {
                  if (g != lf && g != nf) continue;
                  g = ub(b, e, d, k);
                  k[0] = Jf(b, e);
                } else {
                  if (g != jf && g != nf) continue;
                  g = ub(b, e, d, k);
                  for (f = 1; f <= g; f++) k[f] = -k[f];
                  k[0] = -Kf(b, e);
                }
                a: {
                  var m = b;
                  f = d;
                  for (
                    var n = k,
                      q = l,
                      p = null,
                      r = null,
                      t = Array(5),
                      y = void 0,
                      v = void 0,
                      S = void 0,
                      M = (S = void 0),
                      ba = void 0,
                      J = (M = M = void 0),
                      S = 0,
                      v = 1;
                    v <= g;
                    v++
                  )
                    (y = f[v]),
                      rb(m, y) - Ka + gf == cf
                        ? (n[0] -= n[v] * Lf(m, y))
                        : (S++, (f[S] = f[v]), (n[S] = n[v]));
                  g = S;
                  S = 0;
                  for (v = 1; v <= g; v++)
                    (y = f[v]),
                      (Ic(m, y) == Ma ? Mf : Nf) == Nf &&
                        rb(m, y) - Ka + gf == nf &&
                        0 == Lf(m, y) &&
                        1 == Of(m, y) &&
                        (S++,
                        (ba = f[S]),
                        (M = n[S]),
                        (f[S] = f[v]),
                        (n[S] = n[v]),
                        (f[v] = ba),
                        (n[v] = M));
                  if (2 > S) g = 0;
                  else {
                    ba = M = 0;
                    for (v = S + 1; v <= g; v++) {
                      y = f[v];
                      if (rb(m, y) - Ka + gf != nf) {
                        g = 0;
                        break a;
                      }
                      0 < n[v]
                        ? ((ba += n[v] * Lf(m, y)), (M += n[v] * Of(m, y)))
                        : ((ba += n[v] * Of(m, y)), (M += n[v] * Lf(m, y)));
                    }
                    M -= ba;
                    J = 0;
                    for (v = S + 1; v <= g; v++)
                      (y = f[v]), (J += n[v] * Dc(m, y));
                    J -= ba;
                    0 > J && (J = 0);
                    J > M && (J = M);
                    n[0] -= ba;
                    for (v = 1; v <= S; v++)
                      (y = f[v]),
                        (q[v] = Dc(m, y)),
                        0 > q[v] && (q[v] = 0),
                        1 < q[v] && (q[v] = 1);
                    for (v = 1; v <= S; v++)
                      0 > n[v] &&
                        ((f[v] = -f[v]),
                        (n[v] = -n[v]),
                        (n[0] += n[v]),
                        (q[v] = 1 - q[v]));
                    m = S;
                    v = n[0];
                    y = J;
                    J = void 0;
                    for (J = 1; J <= m; J++);
                    for (J = 1; J <= m; J++);
                    J = void 0;
                    b: {
                      for (
                        var ea = (J = void 0),
                          fa = 0,
                          sa = 0,
                          K = void 0,
                          oa = void 0,
                          Bb = 0.001,
                          K = 0.001 * (1 + Math.abs(v)),
                          J = 1;
                        J <= m;
                        J++
                      )
                        for (ea = J + 1; ea <= m; ea++) {
                          fa++;
                          if (1e3 < fa) {
                            J = sa;
                            break b;
                          }
                          n[J] + n[ea] + y > v + K &&
                            ((oa = n[J] + n[ea] - v),
                            (p = 1 / (oa + M)),
                            (r = 2 - p * oa),
                            (oa = q[J] + q[ea] + p * y - r),
                            Bb < oa &&
                              ((Bb = oa), (t[1] = J), (t[2] = ea), (sa = 1)));
                        }
                      J = sa;
                    }
                    ea = void 0;
                    if (J) ea = 2;
                    else {
                      J = J = void 0;
                      b: {
                        for (
                          var fa = (ea = J = void 0),
                            K = (sa = 0),
                            Bb = (oa = void 0),
                            ec = 0.001,
                            oa = 0.001 * (1 + Math.abs(v)),
                            J = 1;
                          J <= m;
                          J++
                        )
                          for (ea = J + 1; ea <= m; ea++)
                            for (fa = ea + 1; fa <= m; fa++) {
                              sa++;
                              if (1e3 < sa) {
                                J = K;
                                break b;
                              }
                              n[J] + n[ea] + n[fa] + y > v + oa &&
                                ((Bb = n[J] + n[ea] + n[fa] - v),
                                (p = 1 / (Bb + M)),
                                (r = 3 - p * Bb),
                                (Bb = q[J] + q[ea] + q[fa] + p * y - r),
                                ec < Bb &&
                                  ((ec = Bb),
                                  (t[1] = J),
                                  (t[2] = ea),
                                  (t[3] = fa),
                                  (K = 1)));
                            }
                        J = K;
                      }
                      if (J) J = 3;
                      else {
                        J = void 0;
                        b: {
                          for (
                            var sa = (fa = ea = J = void 0),
                              oa = (K = 0),
                              ec = (Bb = void 0),
                              Wj = 0.001,
                              Bb = 0.001 * (1 + Math.abs(v)),
                              J = 1;
                            J <= m;
                            J++
                          )
                            for (ea = J + 1; ea <= m; ea++)
                              for (fa = ea + 1; fa <= m; fa++)
                                for (sa = fa + 1; sa <= m; sa++) {
                                  K++;
                                  if (1e3 < K) {
                                    J = oa;
                                    break b;
                                  }
                                  n[J] + n[ea] + n[fa] + n[sa] + y > v + Bb &&
                                    ((ec = n[J] + n[ea] + n[fa] + n[sa] - v),
                                    (p = 1 / (ec + M)),
                                    (r = 4 - p * ec),
                                    (ec =
                                      q[J] + q[ea] + q[fa] + q[sa] + p * y - r),
                                    Wj < ec &&
                                      ((Wj = ec),
                                      (t[1] = J),
                                      (t[2] = ea),
                                      (t[3] = fa),
                                      (t[4] = sa),
                                      (oa = 1)));
                                }
                          J = oa;
                        }
                        J = J ? 4 : 0;
                      }
                      ea = J;
                    }
                    M = ea;
                    if (0 == M) g = 0;
                    else {
                      f[0] = 0;
                      n[0] = r;
                      for (y = 1; y <= M; y++) t[y] = f[t[y]];
                      for (v = 1; v <= M; v++)
                        0 < t[v]
                          ? ((f[v] = +t[v]), (n[v] = 1))
                          : ((f[v] = -t[v]), (n[v] = -1), (n[0] -= 1));
                      for (v = S + 1; v <= g; v++)
                        M++, (f[M] = f[v]), (n[M] = p * n[v]);
                      n[0] += p * ba;
                      g = M;
                    }
                  }
                }
                if (0 != g) {
                  f = b;
                  n = g;
                  q = d;
                  p = k;
                  r = lb(f);
                  S = t = void 0;
                  ba = 0;
                  0 > n &&
                    w("lpx_eval_row: len = " + n + "; invalid row length");
                  for (S = 1; S <= n; S++)
                    (t = q[S]),
                      (1 <= t && t <= r) ||
                        w(
                          "lpx_eval_row: j = " +
                            t +
                            "; column number out of range"
                        ),
                      (ba += p[S] * Dc(f, t));
                  f = ba - k[0];
                  0.001 > f || Id(a, Ff, g, d, k, Ta, k[0]);
                }
              }
          }
          a.p.vd == bb &&
            null != a.Ne &&
            ((0 == a.N.La && 50 > a.N.Rd) || (0 < a.N.La && 5 > a.N.Rd)) &&
            ((c = a.Ne),
            (d = lb(a.A)),
            (b = new Int32Array(1 + d)),
            (d = new Float64Array(1 + d)),
            (c = Pf(a.A, c, b, d)),
            0 < c && Id(a, Gf, c, b, d, Ta, d[0]));
        }
      }
    }
    function p(a) {
      var b,
        d,
        e = 0;
      for (b = a.head; null != b; b = d)
        (d = b.e), c(a, b.s) || (sf(a, b.s), e++);
      a.p.o >= mc &&
        (1 == e
          ? x("One hopeless branch has been pruned")
          : 1 < e && x(e + " hopeless branches have been pruned"));
    }
    var m,
      q,
      r,
      n,
      t = 0,
      y = a.hc;
    for (q = 0; ; ) {
      r = null;
      switch (q) {
        case 0:
          if (null == a.head) {
            a.p.o >= mc && x("Active list is empty!");
            n = 0;
            r = 3;
            break;
          }
          if (
            null != a.p.ob &&
            ((a.reason = Qf), a.p.ob(a, a.p.Uc), (a.reason = 0), a.stop)
          ) {
            n = Rc;
            r = 3;
            break;
          }
          0 == a.gf && (a.gf = 1 == a.Pd ? a.head.s : 0 != a.ud ? a.ud : Rf(a));
          pf(a, a.gf);
          a.gf = a.ud = 0;
          null != a.N.R && a.N.R.s != t && (t = 0);
          m = a.N.s;
          a.p.o >= mc &&
            (x(
              "------------------------------------------------------------------------"
            ),
            x("Processing node " + m + " at level " + a.N.La + ""));
          1 == m &&
            (a.p.Ad == bb && a.p.o >= Wb && x("Gomory's cuts enabled"),
            a.p.Ed == bb &&
              (a.p.o >= Wb && x("MIR cuts enabled"), (a.Xf = Sf(a))),
            a.p.xd == bb && a.p.o >= Wb && x("Cover cuts enabled"),
            a.p.vd == bb &&
              (a.p.o >= Wb && x("Clique cuts enabled"), (a.Ne = Tf(a.A))));
        case 1:
          (a.p.o >= mc || (a.p.o >= fc && a.p.bc - 1 <= 1e3 * la(a.Lg))) &&
            b(a, 0);
          a.p.o >= Wb &&
            60 <= la(y) &&
            (x("Time used: " + la(a.hc) + " secs"), (y = ja()));
          if (0 < a.p.ce && xf(a) <= a.p.ce) {
            a.p.o >= mc &&
              x("Relative gap tolerance reached; search terminated ");
            n = Pc;
            r = 3;
            break;
          }
          if (2147483647 > a.p.sb && a.p.sb - 1 <= 1e3 * la(a.hc)) {
            a.p.o >= mc && x("Time limit exhausted; search terminated");
            n = Qc;
            r = 3;
            break;
          }
          if (
            null != a.p.ob &&
            ((a.reason = Uf), a.p.ob(a, a.p.Uc), (a.reason = 0), a.stop)
          ) {
            n = Rc;
            r = 3;
            break;
          }
          if (a.p.ed != hd)
            if (a.p.ed == id) {
              if (0 == a.N.La && zf(a, 100)) {
                r = 2;
                break;
              }
            } else if (a.p.ed == jd && zf(a, 0 == a.N.La ? 100 : 10)) {
              r = 2;
              break;
            }
          if (!c(a, m)) {
            x("*** not tested yet ***");
            r = 2;
            break;
          }
          a.p.o >= mc && x("Solving LP relaxation...");
          n = yf(a);
          if (0 != n && n != Vf && n != Wf) {
            a.p.o >= Lb &&
              x(
                "ios_driver: unable to solve current LP relaxation; glp_simplex returned " +
                  n +
                  ""
              );
            n = Sb;
            r = 3;
            break;
          }
          q = a.A.na;
          r = a.A.sa;
          if (q == dc && r == dc)
            a.p.o >= mc && x("Found optimal solution to LP relaxation");
          else if (r == jc) {
            a.p.o >= Lb &&
              x(
                "ios_driver: current LP relaxation has no dual feasible solution"
              );
            n = Sb;
            r = 3;
            break;
          } else if (q == Ad && r == dc) {
            a.p.o >= mc &&
              x(
                "LP relaxation has no solution better than incumbent objective value"
              );
            r = 2;
            break;
          } else if (q == jc) {
            a.p.o >= mc && x("LP relaxation has no feasible solution");
            r = 2;
            break;
          }
          q = a.N.rc = a.A.aa;
          q = uf(a, q);
          a.A.dir == za
            ? a.N.bound < q && (a.N.bound = q)
            : a.A.dir == Ea && a.N.bound > q && (a.N.bound = q);
          a.p.o >= mc && x("Local bound is " + q + "");
          if (!c(a, m)) {
            a.p.o >= mc && x("Current branch is hopeless and can be pruned");
            r = 2;
            break;
          }
          if (null != a.p.ob) {
            a.reason = Ga;
            a.p.ob(a, a.p.Uc);
            a.reason = 0;
            if (a.stop) {
              n = Rc;
              r = 3;
              break;
            }
            if (a.pe) {
              a.pe = a.tf = 0;
              r = 1;
              break;
            }
            a.tf && ((a.tf = 0), Jb(a.A));
          }
          d(a);
          if (0 == a.N.Ag) {
            a.p.o >= mc && x("New integer feasible solution found");
            a.p.o >= Wb && k(a);
            e(a);
            a.p.o >= fc && b(a, 1);
            if (
              null != a.p.ob &&
              ((a.reason = Xf), a.p.ob(a, a.p.Uc), (a.reason = 0), a.stop)
            ) {
              n = Rc;
              r = 3;
              break;
            }
            r = 2;
            break;
          }
          a.A.za == dc && g(a);
          if (null != a.p.ob) {
            a.reason = Yf;
            a.p.ob(a, a.p.Uc);
            a.reason = 0;
            if (a.stop) {
              n = Rc;
              r = 3;
              break;
            }
            if (!c(a, m)) {
              a.p.o >= mc &&
                x("Current branch became hopeless and can be pruned");
              r = 2;
              break;
            }
          }
          if (a.p.Xe && ((a.reason = Yf), Zf(a), (a.reason = 0), !c(a, m))) {
            a.p.o >= mc &&
              x("Current branch became hopeless and can be pruned");
            r = 2;
            break;
          }
          if (
            null != a.p.ob &&
            ((a.reason = Ia), a.p.ob(a, a.p.Uc), (a.reason = 0), a.stop)
          ) {
            n = Rc;
            r = 3;
            break;
          }
          if (0 == a.N.La || 0 == t) (a.reason = Ia), l(a), (a.reason = 0);
          0 < a.Bd.size && ((a.reason = Ia), $f(a), (a.reason = 0));
          Oc(a.Bd);
          if (a.pe) {
            a.pe = 0;
            a.N.Rd++;
            r = 1;
            break;
          }
          h(a);
          a.p.o >= Wb && 0 == a.N.La && k(a);
          null != a.Fd && ag(a);
          if (
            null != a.p.ob &&
            ((a.reason = bg), a.p.ob(a, a.p.Uc), (a.reason = 0), a.stop)
          ) {
            n = Rc;
            r = 3;
            break;
          }
          0 == a.Tc &&
            (a.Tc = cg(a, function (b) {
              a.Jf = b;
            }));
          q = a.N.s;
          n = f(a, a.Tc, a.Jf);
          a.Tc = a.Jf = 0;
          if (0 == n) {
            t = q;
            r = 0;
            break;
          } else if (1 == n) {
            a.N.eg = a.N.Rd = 0;
            r = 1;
            break;
          } else if (2 == n) {
            r = 2;
            break;
          }
        case 2:
          a.p.o >= mc && x("Node " + m + " fathomed");
          qf(a);
          sf(a, m);
          a.A.za == dc && p(a);
          r = t = 0;
          break;
        case 3:
          return a.p.o >= fc && b(a, 0), (a.Xf = null), (a.Ne = null), n;
      }
      if (null == r) break;
      q = r;
    }
  }
  function dg(a) {
    var b;
    b = {};
    b.i = a;
    b.L = 0;
    b.Ja = new Int32Array(1 + a);
    b.Z = new Int32Array(1 + a);
    b.j = new Float64Array(1 + a);
    return b;
  }
  function eg(a, b, c) {
    var d = a.Ja[b];
    0 == c
      ? 0 != d &&
        ((a.Ja[b] = 0),
        d < a.L &&
          ((a.Ja[a.Z[a.L]] = d), (a.Z[d] = a.Z[a.L]), (a.j[d] = a.j[a.L])),
        a.L--)
      : (0 == d && ((d = ++a.L), (a.Ja[b] = d), (a.Z[d] = b)), (a.j[d] = c));
  }
  function fg(a) {
    for (var b = 1; b <= a.L; b++) a.Ja[a.Z[b]] = 0;
    a.L = 0;
  }
  function gg(a, b) {
    for (var c = 0, d = 1; d <= a.L; d++)
      0 == Math.abs(a.j[d]) || Math.abs(a.j[d]) < b
        ? (a.Ja[a.Z[d]] = 0)
        : (c++, (a.Ja[a.Z[d]] = c), (a.Z[c] = a.Z[d]), (a.j[c] = a.j[d]));
    a.L = c;
  }
  function hg(a, b) {
    fg(a);
    a.L = b.L;
    ga(a.Z, 1, b.Z, 1, a.L);
    ga(a.j, 1, b.j, 1, a.L);
    for (var c = 1; c <= a.L; c++) a.Ja[a.Z[c]] = c;
  }
  function Hf(a) {
    function b(a) {
      return a - Math.floor(a);
    }
    function c(a, c, d) {
      var e = a.A,
        f = e.g,
        g = e.i,
        h = c.Z,
        k = c.j;
      c = c.kh;
      var l, D, H, R, V, O, Q, F, W, X, ca;
      D = Bd(e, f + d, h, k);
      F = e.f[d].r;
      for (l = 1; l <= f + g; l++) c[l] = 0;
      ca = b(F);
      for (d = 1; d <= D; d++) {
        l = h[d];
        l <= f ? ((R = e.n[l]), (H = Ma)) : ((R = e.f[l - f]), (H = R.kind));
        V = R.c;
        O = R.d;
        R = R.m;
        W = k[d];
        if (1e5 < Math.abs(W)) return;
        if (!(1e-10 > Math.abs(W))) {
          switch (R) {
            case Ra:
              return;
            case G:
              Q = -W;
              break;
            case Ua:
              Q = +W;
              break;
            case Na:
              continue;
          }
          switch (H) {
            case Fc:
              if (1e-10 > Math.abs(Q - Math.floor(Q + 0.5))) continue;
              else X = b(Q) <= b(F) ? b(Q) : (b(F) / (1 - b(F))) * (1 - b(Q));
              break;
            case Ma:
              X = 0 <= Q ? +Q : (b(F) / (1 - b(F))) * -Q;
          }
          switch (R) {
            case G:
              c[l] = +X;
              ca += X * V;
              break;
            case Ua:
              (c[l] = -X), (ca -= X * O);
          }
        }
      }
      for (d = 1; d <= f; d++)
        if (!(1e-10 > Math.abs(c[d])))
          for (R = e.n[d], Q = R.k; null != Q; Q = Q.B)
            c[f + Q.f.C] += c[d] * Q.j;
      D = 0;
      for (d = 1; d <= g; d++)
        1e-10 > Math.abs(c[f + d]) ||
          ((R = e.f[d]),
          R.type == B
            ? (ca -= c[f + d] * R.c)
            : (D++, (h[D] = d), (k[D] = c[f + d])));
      1e-12 > Math.abs(ca) && (ca = 0);
      for (l = 1; l <= D; l++)
        if (0.001 > Math.abs(k[l]) || 1e3 < Math.abs(k[l])) return;
      Id(a, Df, D, h, k, Sa, ca);
    }
    var d = a.A,
      e = d.g,
      f = d.i,
      g,
      h,
      k = {};
    g = Array(1 + f);
    k.Z = new Int32Array(1 + f);
    k.j = new Float64Array(1 + f);
    k.kh = new Float64Array(1 + e + f);
    e = 0;
    for (h = 1; h <= f; h++) {
      var l = d.f[h];
      l.kind == Fc &&
        l.type != B &&
        l.m == A &&
        ((l = b(l.r)),
        0.05 <= l && 0.95 >= l && (e++, (g[e].C = h), (g[e].Nb = l)));
    }
    ma(g, e, function (a, b) {
      return a.Nb > b.Nb ? -1 : a.Nb < b.Nb ? 1 : 0;
    });
    f = Hd(a);
    for (d = 1; d <= e && !(50 <= Hd(a) - f); d++) c(a, k, g[d].C);
  }
  var ig = 0,
    jg = 5,
    kg = 0,
    lg = 1,
    mg = 2;
  function Sf(a) {
    var b = a.A,
      c = b.g,
      b = b.i,
      d;
    ig && x("ios_mir_init: warning: debug mode enabled");
    d = {};
    d.g = c;
    d.i = b;
    d.Eb = new Int8Array(1 + c);
    d.ab = new Int8Array(1 + c + b);
    d.c = new Float64Array(1 + c + b);
    d.Yb = new Int32Array(1 + c + b);
    d.d = new Float64Array(1 + c + b);
    d.zb = new Int32Array(1 + c + b);
    d.x = new Float64Array(1 + c + b);
    d.Ef = new Int32Array(1 + jg);
    d.Ya = dg(c + b);
    d.lb = new Int8Array(1 + c + b);
    d.oa = dg(c + b);
    d.G = dg(c + b);
    (function (a, b) {
      var c = a.A,
        d = b.g,
        k;
      for (k = 1; k <= d; k++) {
        var l = c.n[k];
        b.Eb[k] = 0;
        b.ab[k] = 0;
        switch (l.type) {
          case Ka:
            b.c[k] = -s;
            b.d[k] = +s;
            break;
          case Sa:
            b.c[k] = l.c;
            b.d[k] = +s;
            break;
          case Ta:
            b.c[k] = -s;
            b.d[k] = l.d;
            break;
          case I:
            b.c[k] = l.c;
            b.d[k] = l.d;
            break;
          case B:
            b.c[k] = b.d[k] = l.c;
        }
        b.Yb[k] = b.zb[k] = 0;
      }
    })(a, d);
    (function (a, b) {
      var c = a.A,
        d = b.g,
        k = b.i,
        l;
      for (l = d + 1; l <= d + k; l++) {
        var p = c.f[l - d];
        switch (p.kind) {
          case Ma:
            b.ab[l] = 0;
            break;
          case Fc:
            b.ab[l] = 1;
        }
        switch (p.type) {
          case Ka:
            b.c[l] = -s;
            b.d[l] = +s;
            break;
          case Sa:
            b.c[l] = p.c;
            b.d[l] = +s;
            break;
          case Ta:
            b.c[l] = -s;
            b.d[l] = p.d;
            break;
          case I:
            b.c[l] = p.c;
            b.d[l] = p.d;
            break;
          case B:
            b.c[l] = b.d[l] = p.c;
        }
        b.Yb[l] = b.zb[l] = 0;
      }
    })(a, d);
    (function (a, b) {
      var c = a.A,
        d = b.g,
        k,
        l,
        p,
        m,
        q,
        r;
      for (l = 1; l <= d; l++)
        if ((0 == b.c[l] && b.d[l] == +s) || (b.c[l] == -s && 0 == b.d[l]))
          if (
            ((k = c.n[l].k),
            null != k &&
              ((p = d + k.f.C),
              (q = k.j),
              (k = k.B),
              null != k && ((m = d + k.f.C), (r = k.j), null == k.B)))
          ) {
            if (b.ab[p] || !b.ab[m])
              if (b.ab[p] && !b.ab[m])
                (m = p), (r = q), (p = d + k.f.C), (q = k.j);
              else continue;
            b.c[m] != -s &&
              b.d[m] != +s &&
              b.c[m] != b.d[m] &&
              (0 == b.d[l] && ((q = -q), (r = -r)),
              0 < q
                ? 0 == b.Yb[p] &&
                  ((b.c[p] = -r / q), (b.Yb[p] = m), (b.Eb[l] = 1))
                : 0 == b.zb[p] &&
                  ((b.d[p] = -r / q), (b.zb[p] = m), (b.Eb[l] = 1)));
          }
    })(a, d);
    (function (a, b) {
      var c = a.A,
        d = b.g,
        k,
        l,
        p,
        m;
      for (l = 1; l <= d; l++)
        if (b.c[l] == -s && b.d[l] == +s) b.Eb[l] = 1;
        else {
          m = 0;
          for (k = c.n[l].k; null != k; k = k.B) {
            p = d + k.f.C;
            if (b.c[p] == -s && b.d[p] == +s) {
              b.Eb[l] = 1;
              break;
            }
            if ((b.ab[p] && b.c[p] == -s) || (b.ab[p] && b.d[p] == +s)) {
              b.Eb[l] = 1;
              break;
            }
            (0 == b.Yb[p] && 0 == b.zb[p] && b.c[p] == b.d[p]) || m++;
          }
          0 == m && (b.Eb[l] = 1);
        }
    })(a, d);
    return d;
  }
  function If(a, b) {
    function c(a, b, c, d, e, f, g) {
      function k(a, b, c, d, e, f, g) {
        var h;
        h = c;
        for (c = 1; c <= a; c++)
          (g[c] = b[c] / f), e[c] && (g[c] = -g[c]), (h -= b[c] * d[c]);
        b = h / f;
        var l;
        if (0.01 > Math.abs(b - Math.floor(b + 0.5))) b = 1;
        else {
          h = b - Math.floor(b);
          for (c = 1; c <= a; c++)
            (l = g[c] - Math.floor(g[c]) - h),
              (g[c] =
                0 >= l ? Math.floor(g[c]) : Math.floor(g[c]) + l / (1 - h));
          q = Math.floor(b);
          r = 1 / (1 - h);
          b = 0;
        }
        if (b) return 1;
        for (c = 1; c <= a; c++) e[c] && ((g[c] = -g[c]), (q += g[c] * d[c]));
        r /= f;
        return 0;
      }
      var h, l, m, n, p;
      m = Array(4);
      var t, v, y, D;
      y = new Int8Array(1 + a);
      D = Array(1 + a);
      for (l = 1; l <= a; l++) y[l] = e[l] >= 0.5 * d[l];
      v = p = 0;
      for (l = 1; l <= a; l++)
        if (
          ((h = 1e-9 * (1 + Math.abs(d[l]))),
          !(
            e[l] < h ||
            e[l] > d[l] - h ||
            ((h = k(a, b, c, d, y, Math.abs(b[l]), g)), h)
          ))
        ) {
          t = -q - r * f;
          for (h = 1; h <= a; h++) t += g[h] * e[h];
          v < t && ((v = t), (p = Math.abs(b[l])));
        }
      0.001 > v && (v = 0);
      if (0 == v) return v;
      m[1] = p / 2;
      m[2] = p / 4;
      m[3] = p / 8;
      for (l = 1; 3 >= l; l++)
        if (((h = k(a, b, c, d, y, m[l], g)), !h)) {
          t = -q - r * f;
          for (h = 1; h <= a; h++) t += g[h] * e[h];
          v < t && ((v = t), (p = m[l]));
        }
      m = 0;
      for (l = 1; l <= a; l++)
        (h = 1e-9 * (1 + Math.abs(d[l]))),
          e[l] < h ||
            e[l] > d[l] - h ||
            (m++, (D[m].C = l), (D[m].xf = Math.abs(e[l] - 0.5 * d[l])));
      ma(D, m, function (a, b) {
        return a.xf < b.xf ? -1 : a.xf > b.xf ? 1 : 0;
      });
      for (n = 1; n <= m; n++)
        if (
          ((l = D[n].C),
          (y[l] = !y[l]),
          (h = k(a, b, c, d, y, p, g)),
          (y[l] = !y[l]),
          !h)
        ) {
          t = -q - r * f;
          for (h = 1; h <= a; h++) t += g[h] * e[h];
          v < t && ((v = t), (y[l] = !y[l]));
        }
      h = k(a, b, c, d, y, p, g);
      return v;
    }
    function d(a, b, c) {
      var d = a.A;
      a = b.g;
      b.Eb[c] = 2;
      b.Je = 1;
      b.Ef[1] = c;
      fg(b.Ya);
      eg(b.Ya, c, 1);
      for (c = d.n[c].k; null != c; c = c.B) eg(b.Ya, a + c.f.C, -c.j);
      b.Df = 0;
    }
    function e(a) {
      var b, c;
      for (b = 1; b <= a.Ya.L; b++)
        (c = a.Ya.Z[b]),
          0 == a.Yb[c] &&
            0 == a.zb[c] &&
            a.c[c] == a.d[c] &&
            ((a.Df -= a.Ya.j[b] * a.c[c]), (a.Ya.j[b] = 0));
      gg(a.Ya, 2.220446049250313e-16);
    }
    function f(a) {
      var b, c, d, e;
      for (b = 1; b <= a.Ya.L; b++)
        (c = a.Ya.Z[b]),
          a.ab[c] ||
            ((d = a.Yb[c]),
            (e =
              0 == d
                ? a.c[c] == -s
                  ? s
                  : a.x[c] - a.c[c]
                : a.x[c] - a.c[c] * a.x[d]),
            (d = a.zb[c]),
            (d =
              0 == d
                ? a.zb[c] == +s
                  ? s
                  : a.d[c] - a.x[c]
                : a.d[c] * a.x[d] - a.x[c]),
            (a.lb[c] = e <= d ? lg : mg));
    }
    function g(a) {
      var b, c, d, e;
      hg(a.oa, a.Ya);
      a.tc = a.Df;
      for (b = a.oa.L; 1 <= b; b--)
        (d = a.oa.Z[b]),
          a.ab[d] ||
            (a.lb[d] == lg
              ? ((e = a.Yb[d]),
                0 == e
                  ? (a.tc -= a.oa.j[b] * a.c[d])
                  : ((c = a.oa.Ja[e]),
                    0 == c &&
                      (eg(a.oa, e, 1), (c = a.oa.Ja[e]), (a.oa.j[c] = 0)),
                    (a.oa.j[c] += a.oa.j[b] * a.c[d])))
              : a.lb[d] == mg &&
                ((e = a.zb[d]),
                0 == e
                  ? (a.tc -= a.oa.j[b] * a.d[d])
                  : ((c = a.oa.Ja[e]),
                    0 == c &&
                      (eg(a.oa, e, 1), (c = a.oa.Ja[e]), (a.oa.j[c] = 0)),
                    (a.oa.j[c] += a.oa.j[b] * a.d[d])),
                (a.oa.j[b] = -a.oa.j[b])));
      for (b = 1; b <= a.oa.L; b++)
        (d = a.oa.Z[b]),
          a.ab[d] &&
            (Math.abs(a.c[d]) <= Math.abs(a.d[d])
              ? ((a.lb[d] = lg), (a.tc -= a.oa.j[b] * a.c[d]))
              : ((a.lb[d] = mg),
                (a.tc -= a.oa.j[b] * a.d[d]),
                (a.oa.j[b] = -a.oa.j[b])));
    }
    function h(a) {
      var b = a.g,
        d = a.i,
        e,
        f,
        g,
        h,
        k,
        l,
        m;
      k = 0;
      hg(a.G, a.oa);
      a.Lb = a.tc;
      gg(a.G, 2.220446049250313e-16);
      for (e = 1; e <= a.G.L; e++)
        (f = a.G.Z[e]), !a.ab[f] && 0 < a.G.j[e] && (a.G.j[e] = 0);
      gg(a.G, 0);
      h = 0;
      for (e = 1; e <= a.G.L; e++)
        (f = a.G.Z[e]),
          a.ab[f] &&
            (h++,
            (g = a.G.Z[h]),
            (a.G.Ja[f] = h),
            (a.G.Ja[g] = e),
            (a.G.Z[h] = f),
            (a.G.Z[e] = g),
            (f = a.G.j[h]),
            (a.G.j[h] = a.G.j[e]),
            (a.G.j[e] = f));
      if (0 == h) return k;
      l = new Float64Array(1 + h);
      g = new Float64Array(1 + h);
      m = new Float64Array(1 + h);
      for (e = 1; e <= h; e++)
        (f = a.G.Z[e]),
          (l[e] = a.d[f] - a.c[f]),
          a.lb[f] == lg
            ? (g[e] = a.x[f] - a.c[f])
            : a.lb[f] == mg && (g[e] = a.d[f] - a.x[f]),
          0 > g[e] && (g[e] = 0);
      k = 0;
      for (e = h + 1; e <= a.G.L; e++)
        (f = a.G.Z[e]),
          a.lb[f] == lg
            ? ((g = a.Yb[f]),
              (g = 0 == g ? a.x[f] - a.c[f] : a.x[f] - a.c[f] * a.x[g]))
            : a.lb[f] == mg &&
              ((g = a.zb[f]),
              (g = 0 == g ? a.d[f] - a.x[f] : a.d[f] * a.x[g] - a.x[f])),
          0 > g && (g = 0),
          (k -= a.G.j[e] * g);
      k = c(h, a.G.j, a.Lb, l, g, k, m);
      if (0 == k) return k;
      for (e = 1; e <= h; e++) a.G.j[e] = m[e];
      for (e = h + 1; e <= a.G.L; e++)
        (f = a.G.Z[e]), f <= b + d && (a.G.j[e] *= 0);
      a.Lb = null;
      return k;
    }
    function k(a) {
      var b, c, d, e;
      for (b = 1; b <= a.G.L; b++)
        (d = a.G.Z[b]),
          a.ab[d] &&
            (a.lb[d] == lg
              ? (a.Lb += a.G.j[b] * a.c[d])
              : a.lb[d] == mg &&
                ((a.Lb -= a.G.j[b] * a.d[d]), (a.G.j[b] = -a.G.j[b])));
      for (b = 1; b <= a.G.L; b++)
        (d = a.G.Z[b]),
          a.ab[d] ||
            (a.lb[d] == lg
              ? ((e = a.Yb[d]),
                0 == e
                  ? (a.Lb += a.G.j[b] * a.c[d])
                  : ((c = a.G.Ja[e]),
                    0 == c && (eg(a.G, e, 1), (c = a.G.Ja[e]), (a.G.j[c] = 0)),
                    (a.G.j[c] -= a.G.j[b] * a.c[d])))
              : a.lb[d] == mg &&
                ((e = a.zb[d]),
                0 == e
                  ? (a.Lb -= a.G.j[b] * a.d[d])
                  : ((c = a.G.Ja[e]),
                    0 == c && (eg(a.G, e, 1), (c = a.G.Ja[e]), (a.G.j[c] = 0)),
                    (a.G.j[c] += a.G.j[b] * a.d[d])),
                (a.G.j[b] = -a.G.j[b])));
    }
    function l(a, b) {
      var c = a.A,
        d = b.g,
        e,
        f,
        g,
        h;
      for (f = b.G.L; 1 <= f; f--)
        if (((e = b.G.Z[f]), !(e > d))) {
          for (e = c.n[e].k; null != e; e = e.B)
            (g = d + e.f.C),
              (h = b.G.Ja[g]),
              0 == h && (eg(b.G, g, 1), (h = b.G.Ja[g]), (b.G.j[h] = 0)),
              (b.G.j[h] += b.G.j[f] * e.j);
          b.G.j[f] = 0;
        }
      gg(b.G, 0);
    }
    function p(a, b) {
      var c = b.g,
        d = b.i,
        e,
        f,
        g = new Int32Array(1 + d),
        h = new Float64Array(1 + d);
      f = 0;
      for (d = b.G.L; 1 <= d; d--)
        (e = b.G.Z[d]), f++, (g[f] = e - c), (h[f] = b.G.j[d]);
      Id(a, Ef, f, g, h, Ta, b.Lb);
    }
    function m(a, b) {
      var c = a.A,
        d = b.g,
        e = b.i,
        f,
        g,
        h,
        k = 0,
        l = 0,
        m,
        n = 0;
      for (f = 1; f <= b.Ya.L; f++)
        (g = b.Ya.Z[f]),
          g <= d ||
            b.ab[g] ||
            0.001 > Math.abs(b.Ya.j[f]) ||
            ((h = b.Yb[g]),
            (m =
              0 == h
                ? b.c[g] == -s
                  ? s
                  : b.x[g] - b.c[g]
                : b.x[g] - b.c[g] * b.x[h]),
            (h = b.zb[g]),
            (h =
              0 == h
                ? b.zb[g] == +s
                  ? s
                  : b.d[g] - b.x[g]
                : b.d[g] * b.x[h] - b.x[g]),
            (m = m <= h ? m : h),
            !(0.001 > m) && n < m && ((n = m), (k = g)));
      if (0 == k) return 1;
      for (g = 1; g <= d; g++)
        if (!b.Eb[g]) {
          for (f = c.n[g].k; null != f && f.f.C != k - d; f = f.B);
          if (null != f && 0.001 <= Math.abs(f.j)) break;
        }
      if (g > d) return 2;
      b.Je++;
      b.Ef[b.Je] = g;
      b.Eb[g] = 2;
      e = dg(d + e);
      eg(e, g, 1);
      for (f = c.n[g].k; null != f; f = f.B) eg(e, d + f.f.C, -f.j);
      f = b.Ya.Ja[k];
      c = b.Ya;
      d = -b.Ya.j[f] / e.j[e.Ja[k]];
      for (g = 1; g <= e.L; g++)
        (f = e.Z[g]),
          (n = void 0),
          (n = c.Ja[f]),
          (n = 0 == n ? 0 : c.j[n]),
          (m = e.j[g]),
          eg(c, f, n + d * m);
      eg(b.Ya, k, 0);
      return l;
    }
    var q,
      r,
      n = b.g,
      t = b.i,
      y,
      E,
      C;
    (function (a, b) {
      var c = a.A,
        d = b.g,
        e = b.i,
        f;
      for (f = 1; f <= d; f++) b.x[f] = c.n[f].r;
      for (f = d + 1; f <= d + e; f++) b.x[f] = c.f[f - d].r;
    })(a, b);
    ha(b.lb, 1, kg, n + t);
    for (y = 1; y <= n; y++)
      if (!b.Eb[y]) {
        for (d(a, b, y); ; ) {
          e(b);
          if (ig) for (E = 1; E <= n + t; E++);
          f(b);
          g(b);
          C = h(b);
          0 < C && (k(b), l(a, b), p(a, b));
          for (var D = 1; D <= b.oa.L; D++) (E = b.oa.Z[D]), (b.lb[E] = kg);
          if (!(0 == C && b.Je < jg && 0 == m(a, b))) break;
        }
        for (E = 1; E <= b.Je; E++) (C = b.Ef[E]), (b.Eb[C] = 0);
      }
  }
  function Tf(a) {
    function b(a, b) {
      var c;
      switch (ob(a, b) - Ka + gf) {
        case gf:
        case lf:
          c = -s;
          break;
        case jf:
        case nf:
        case cf:
          c = Kf(a, b);
      }
      return c;
    }
    function c(a, b) {
      var c;
      switch (ob(a, b) - Ka + gf) {
        case gf:
        case jf:
          c = +s;
          break;
        case lf:
        case nf:
        case cf:
          c = Jf(a, b);
      }
      return c;
    }
    function d(a, b) {
      var c;
      switch (rb(a, b) - Ka + gf) {
        case gf:
        case lf:
          c = -s;
          break;
        case jf:
        case nf:
        case cf:
          c = Lf(a, b);
      }
      return c;
    }
    function e(a, b) {
      var c;
      switch (rb(a, b) - Ka + gf) {
        case gf:
        case jf:
          c = +s;
          break;
        case lf:
        case nf:
        case cf:
          c = Of(a, b);
      }
      return c;
    }
    function f(a, b) {
      return (
        (Ic(a, b) == Ma ? Mf : Nf) == Nf &&
        rb(a, b) - Ka + gf == nf &&
        0 == Lf(a, b) &&
        1 == Of(a, b)
      );
    }
    function g(a, b, c, f) {
      var g, h, k;
      k = 0;
      for (h = 1; h <= b; h++)
        if (((g = c[h]), 0 < f[h])) {
          g = d(a, g);
          if (g == -s) {
            k = -s;
            break;
          }
          k += f[h] * g;
        } else if (0 > f[h]) {
          g = e(a, g);
          if (g == +s) {
            k = -s;
            break;
          }
          k += f[h] * g;
        }
      return k;
    }
    function h(a, b, c, f) {
      var g, h, k;
      k = 0;
      for (h = 1; h <= b; h++)
        if (((g = c[h]), 0 < f[h])) {
          g = e(a, g);
          if (g == +s) {
            k = +s;
            break;
          }
          k += f[h] * g;
        } else if (0 > f[h]) {
          g = d(a, g);
          if (g == -s) {
            k = +s;
            break;
          }
          k += f[h] * g;
        }
      return k;
    }
    function k(a, b, c, d, e, f, g, h) {
      b != -s && g && (b -= a[f]);
      c != +s && g && (c -= a[f]);
      d != -s && (0 > a[f] && (d -= a[f]), 0 > a[h] && (d -= a[h]));
      e != +s && (0 < a[f] && (e -= a[f]), 0 < a[h] && (e -= a[h]));
      f =
        0 < a[h]
          ? b == -s || e == +s
            ? -s
            : (b - e) / a[h]
          : c == +s || d == -s
          ? -s
          : (c - d) / a[h];
      if (0.001 < f) return 2;
      f =
        0 < a[h]
          ? c == +s || d == -s
            ? +s
            : (c - d) / a[h]
          : b == -s || e == +s
          ? +s
          : (b - e) / a[h];
      return 0.999 > f ? 1 : 0;
    }
    var l = null,
      p,
      m,
      q,
      r,
      n,
      t,
      y,
      E,
      C,
      D,
      H,
      R,
      V,
      O,
      Q,
      F;
    x("Creating the conflict graph...");
    p = kb(a);
    m = lb(a);
    q = 0;
    D = new Int32Array(1 + m);
    H = new Int32Array(1 + m);
    C = new Int32Array(1 + m);
    F = new Float64Array(1 + m);
    for (r = 1; r <= p; r++)
      if (((R = b(a, r)), (V = c(a, r)), R != -s || V != +s))
        if (((E = ub(a, r, C, F)), !(500 < E)))
          for (O = g(a, E, C, F), Q = h(a, E, C, F), t = 1; t <= E; t++)
            if (f(a, C[t]))
              for (y = t + 1; y <= E; y++)
                f(a, C[y]) &&
                  (k(F, R, V, O, Q, t, 0, y) || k(F, R, V, O, Q, t, 1, y)) &&
                  ((n = C[t]),
                  0 == D[n] && (q++, (D[n] = q), (H[q] = n)),
                  (n = C[y]),
                  0 == D[n] && (q++, (D[n] = q), (H[q] = n)));
    if (0 == q || 4e3 < q)
      return x("The conflict graph is either empty or too big"), l;
    l = {};
    l.i = m;
    l.Cb = q;
    l.Dg = 0;
    l.yf = D;
    l.ge = H;
    E = q + q;
    E = ((E * (E - 1)) / 2 + 0) / 1;
    l.Jc = Array(E);
    for (n = 1; n <= q; n++) ng(l, +H[n], -H[n]);
    for (r = 1; r <= p; r++)
      if (((R = b(a, r)), (V = c(a, r)), R != -s || V != +s))
        if (((E = ub(a, r, C, F)), !(500 < E)))
          for (O = g(a, E, C, F), Q = h(a, E, C, F), t = 1; t <= E; t++)
            if (f(a, C[t]))
              for (y = t + 1; y <= E; y++)
                if (f(a, C[y])) {
                  switch (k(F, R, V, O, Q, t, 0, y)) {
                    case 1:
                      ng(l, -C[t], +C[y]);
                      break;
                    case 2:
                      ng(l, -C[t], -C[y]);
                  }
                  switch (k(F, R, V, O, Q, t, 1, y)) {
                    case 1:
                      ng(l, +C[t], +C[y]);
                      break;
                    case 2:
                      ng(l, +C[t], -C[y]);
                  }
                }
    x("The conflict graph has 2*" + l.Cb + " vertices and " + l.Dg + " edges");
    return l;
  }
  function ng(a, b, c) {
    var d;
    0 < b ? (b = a.yf[b]) : ((b = a.yf[-b]), (b += a.Cb));
    0 < c ? (c = a.yf[c]) : ((c = a.yf[-c]), (c += a.Cb));
    b < c && ((d = b), (b = c), (c = d));
    d = ((b - 1) * (b - 2)) / 2 + (c - 1);
    a.Jc[d / 1] |= 1 << (0 - (d % 1));
    a.Dg++;
  }
  function Pf(a, b, c, d) {
    function e(a, b, c) {
      return b == c
        ? 0
        : b > c
        ? f(a, (b * (b - 1)) / 2 + c)
        : f(a, (c * (c - 1)) / 2 + b);
    }
    function f(a, b) {
      return a.Jc[b / 1] & (1 << (0 - (b % 1)));
    }
    function g(a, b, c, d, f, h) {
      var k, l, m, n, q, p, r, ca;
      ca = new Int32Array(a.i);
      if (0 >= b) {
        if ((0 == b && ((a.set[d++] = c[0]), (f += h)), f > a.gd))
          for (a.gd = f, a.bg = d, k = 0; k < d; k++) a.oh[k + 1] = a.set[k];
      } else
        for (k = b; 0 <= k && !(0 == d && k < b); k--) {
          m = c[k];
          if (0 < d && a.wg[m] <= a.gd - f) break;
          a.set[d] = m;
          n = f + a.Ic[m + 1];
          h -= a.Ic[m + 1];
          if (h <= a.gd - n) break;
          for (q = r = p = 0; r < c + k; )
            (l = c[r]),
              r++,
              e(a, l, m) && ((ca[p] = l), p++, (q += a.Ic[l + 1]));
          q <= a.gd - n || g(a, p - 1, ca, d + 1, n, q);
        }
    }
    var h = lb(a),
      k,
      l,
      p,
      m = 0,
      q,
      r,
      n;
    p = new Int32Array(1 + 2 * b.Cb);
    q = new Int32Array(1 + 2 * b.Cb);
    n = new Float64Array(1 + h);
    for (l = 1; l <= b.Cb; l++)
      (k = b.ge[l]),
        (k = Dc(a, k)),
        (k = (100 * k + 0.5) | 0),
        0 > k && (k = 0),
        100 < k && (k = 100),
        (p[l] = k),
        (p[b.Cb + l] = 100 - k);
    p = (function (a, b, c, d) {
      var f = {},
        h,
        k,
        l,
        m,
        n,
        q;
      f.i = a;
      f.Ic = b;
      f.Jc = c;
      f.gd = 0;
      f.bg = 0;
      f.oh = d;
      f.wg = new Int32Array(f.i);
      f.set = new Int32Array(f.i);
      n = new Int32Array(f.i);
      q = new Int32Array(f.i);
      b = new Int32Array(f.i);
      c = ja();
      for (a = 0; a < f.i; a++)
        for (h = q[a] = 0; h < f.i; h++) e(f, a, h) && (q[a] += f.Ic[h + 1]);
      for (a = 0; a < f.i; a++) n[a] = 0;
      for (a = f.i - 1; 0 <= a; a--) {
        m = l = -1;
        for (h = 0; h < f.i; h++)
          !n[h] &&
            (f.Ic[h + 1] > l || (f.Ic[h + 1] == l && q[h] > m)) &&
            ((l = f.Ic[h + 1]), (m = q[h]), (k = h));
        b[a] = k;
        n[k] = 1;
        for (h = 0; h < f.i; h++)
          !n[h] && h != k && e(f, k, h) && (q[h] -= f.Ic[k + 1]);
      }
      for (a = k = 0; a < f.i; a++)
        (k += f.Ic[b[a] + 1]),
          g(f, a, b, 0, 0, k),
          (f.wg[b[a]] = f.gd),
          4.999 <= la(c) &&
            (x("level = " + a + 1 + " (" + f.i + "); best = " + f.gd + ""),
            (c = ja()));
      for (a = 1; a <= f.bg; a++) d[a]++;
      return f.bg;
    })(2 * b.Cb, p, b.Jc, q);
    r = 0;
    for (l = 1; l <= p; l++)
      (k = q[l]),
        k <= b.Cb
          ? ((k = b.ge[k]), (k = Dc(a, k)), (r += k))
          : ((k = b.ge[k - b.Cb]), (k = Dc(a, k)), (r += 1 - k));
    if (1.01 <= r) {
      for (l = a = 1; l <= p; l++)
        (k = q[l]),
          k <= b.Cb
            ? ((k = b.ge[k]), (n[k] += 1))
            : ((k = b.ge[k - b.Cb]), (n[k] -= 1), (a -= 1));
      for (k = 1; k <= h; k++) 0 != n[k] && (m++, (c[m] = k), (d[m] = n[k]));
      c[0] = 0;
      d[0] = a;
    }
    return m;
  }
  function cg(a, b) {
    var c;
    if (a.p.Jb == Zc) {
      var d, e;
      for (d = 1; d <= a.i && !a.ad[d]; d++);
      e = Dc(a.A, d);
      b(e - Math.floor(e) < Math.ceil(e) - e ? Bf : Cf);
      c = d;
    } else if (a.p.Jb == $c) {
      for (d = a.i; 1 <= d && !a.ad[d]; d--);
      e = Dc(a.A, d);
      b(e - Math.floor(e) < Math.ceil(e) - e ? Bf : Cf);
      c = d;
    } else if (a.p.Jb == ad) c = og(a, b);
    else if (a.p.Jb == bd) {
      c = a.A;
      var f = c.g,
        g = c.i,
        h = a.ad,
        k,
        l,
        p,
        m,
        q,
        r,
        n,
        t,
        y,
        E,
        C,
        D,
        H,
        R;
      xc(c);
      n = new Int32Array(1 + g);
      R = new Float64Array(1 + g);
      l = 0;
      H = -1;
      for (k = 1; k <= g; k++)
        if (h[k]) {
          t = Dc(c, k);
          r = Bd(c, f + k, n, R);
          for (q = -1; 1 >= q; q += 2) {
            p = Fd(c, r, n, R, q, 1e-9);
            0 != p && (p = n[p]);
            if (0 == p) p = a.A.dir == za ? +s : -s;
            else {
              for (m = 1; m <= r && n[m] != p; m++);
              m = R[m];
              y = (0 > q ? Math.floor(t) : Math.ceil(t)) - t;
              y /= m;
              p > f &&
                Ic(c, p - f) != Ma &&
                0.001 < Math.abs(y - Math.floor(y + 0.5)) &&
                (y = 0 < y ? Math.ceil(y) : Math.floor(y));
              p <= f
                ? ((m = zc(c, p)), (p = Bc(c, p)))
                : ((m = Cc(c, p - f)), (p = Ec(c, p - f)));
              switch (a.A.dir) {
                case za:
                  if ((m == G && 0 > p) || (m == Ua && 0 < p) || m == Ra) p = 0;
                  break;
                case Ea:
                  if ((m == G && 0 < p) || (m == Ua && 0 > p) || m == Ra) p = 0;
              }
              p *= y;
            }
            0 > q ? (e = p) : (E = p);
          }
          if (H < Math.abs(e) || H < Math.abs(E))
            if (
              ((l = k),
              Math.abs(e) < Math.abs(E)
                ? ((d = Bf), (H = Math.abs(E)))
                : ((d = Cf), (H = Math.abs(e))),
              (C = e),
              (D = E),
              H == s)
            )
              break;
        }
      H < 1e-6 * (1 + 0.001 * Math.abs(c.aa))
        ? (l = og(a, b))
        : (a.p.o >= mc &&
            (x("branch_drtom: column " + l + " chosen to branch on"),
            Math.abs(C) == s
              ? x("branch_drtom: down-branch is infeasible")
              : x("branch_drtom: down-branch bound is " + (yc(c) + C) + ""),
            Math.abs(D) == s
              ? x("branch_drtom: up-branch   is infeasible")
              : x("branch_drtom: up-branch   bound is " + (yc(c) + D) + "")),
          b(d));
      c = l;
    } else a.p.Jb == cd && (c = pg(a, b));
    return c;
  }
  function og(a, b) {
    var c, d, e, f, g, h;
    d = 0;
    g = s;
    for (c = 1; c <= a.i; c++)
      a.ad[c] &&
        ((f = Dc(a.A, c)),
        (h = Math.floor(f) + 0.5),
        g > Math.abs(f - h) &&
          ((d = c), (g = Math.abs(f - h)), (e = f < h ? Bf : Cf)));
    b(e);
    return d;
  }
  function qg(a) {
    a = a.i;
    var b,
      c = {};
    c.yd = new Int32Array(1 + a);
    c.Ud = new Float64Array(1 + a);
    c.Id = new Int32Array(1 + a);
    c.ye = new Float64Array(1 + a);
    for (b = 1; b <= a; b++) (c.yd[b] = c.Id[b] = 0), (c.Ud[b] = c.ye[b] = 0);
    return c;
  }
  function ag(a) {
    var b,
      c,
      d = a.Fd;
    null != a.N.R &&
      ((b = a.N.R.Tc),
      (c = a.A.f[b].r - a.N.R.tg),
      (a = a.A.aa - a.N.R.rc),
      (a = Math.abs(a / c)),
      0 > c ? (d.yd[b]++, (d.Ud[b] += a)) : (d.Id[b]++, (d.ye[b] += a)));
  }
  function pg(a, b) {
    function c(a, b, c) {
      var d, e;
      xc(a);
      d = Ba();
      hb(d, a, 0);
      Wa(d, b, B, c, c);
      b = new kc();
      b.o = lc;
      b.cb = Tb;
      b.oc = 30;
      b.fb = 1e3;
      b.cb = Tb;
      b = sc(d, b);
      0 == b || b == rg
        ? tc(d) == jc
          ? (e = s)
          : uc(d) == dc
          ? (a.dir == za ? (e = d.aa - a.aa) : a.dir == Ea && (e = a.aa - d.aa),
            e < 1e-6 * (1 + 0.001 * Math.abs(a.aa)) && (e = 0))
          : (e = 0)
        : (e = 0);
      return e;
    }
    function d(a, b, d) {
      var e = a.Fd,
        f;
      if (d == Bf) {
        if (0 == e.yd[b]) {
          d = a.A.f[b].r;
          a = c(a.A, b, Math.floor(d));
          if (a == s) return (f = s);
          e.yd[b] = 1;
          e.Ud[b] = a / (d - Math.floor(d));
        }
        f = e.Ud[b] / e.yd[b];
      } else if (d == Cf) {
        if (0 == e.Id[b]) {
          d = a.A.f[b].r;
          a = c(a.A, b, Math.ceil(d));
          if (a == s) return (f = s);
          e.Id[b] = 1;
          e.ye[b] = a / (Math.ceil(d) - d);
        }
        f = e.ye[b] / e.Id[b];
      }
      return f;
    }
    function e(a) {
      var b = a.Fd,
        c,
        d = 0,
        e = 0;
      for (c = 1; c <= a.i; c++)
        Jd(a, c) && (d++, 0 < b.yd[c] && 0 < b.Id[c] && e++);
      x("Pseudocosts initialized for " + e + " of " + d + " variables");
    }
    var f = ja(),
      g,
      h,
      k,
      l,
      p,
      m,
      q;
    null == a.Fd && (a.Fd = qg(a));
    h = 0;
    q = -1;
    for (g = 1; g <= a.i; g++)
      if (Jd(a, g)) {
        l = a.A.f[g].r;
        p = d(a, g, Bf);
        if (p == s) return (h = g), (k = Bf), b(k), h;
        m = p * (l - Math.floor(l));
        p = d(a, g, Cf);
        if (p == s) return (h = g), (k = Cf), b(k), h;
        l = p * (Math.ceil(l) - l);
        p = m > l ? m : l;
        q < p && ((q = p), (h = g), (k = m <= l ? Bf : Cf));
        a.p.o >= bb && 10 <= la(f) && (e(a), (f = ja()));
      }
    if (0 == q) return (h = og(a, b));
    b(k);
    return h;
  }
  function Zf(a) {
    var b = a.A,
      c = b.i,
      d = null,
      e = null,
      f = null,
      g,
      h,
      k,
      l,
      p,
      m,
      q,
      r;
    for (l = 0; ; ) {
      var n = null;
      switch (l) {
        case 0:
          xc(b);
          if (0 != a.N.La || 1 != a.N.eg) {
            n = 5;
            break;
          }
          q = 0;
          for (k = 1; k <= c; k++)
            if (((g = b.f[k]), g.kind != Ma && g.type != B))
              if (g.type == I && 0 == g.c && 1 == g.d) q++;
              else {
                a.p.o >= Wb &&
                  x(
                    "FPUMP heuristic cannot be applied due to general integer variables"
                  );
                n = 5;
                break;
              }
          if (null != n) break;
          if (0 == q) {
            n = 5;
            break;
          }
          a.p.o >= Wb && x("Applying FPUMP heuristic...");
          e = Array(1 + q);
          ia(e, 1, q);
          l = 0;
          for (k = 1; k <= c; k++)
            (g = b.f[k]), g.kind == Fc && g.type == I && (e[++l].C = k);
          d = Ba();
        case 1:
          hb(d, b, cb);
          if (b.za == dc) {
            La(d, 1);
            p = new Int32Array(1 + c);
            m = new Float64Array(1 + c);
            for (k = 1; k <= c; k++) (p[k] = k), (m[k] = b.f[k].u);
            Ya(d, d.g, c, p, m);
            p = 0.1 * b.aa + 0.9 * b.ta;
            b.dir == za
              ? Va(d, d.g, Ta, 0, p - b.ha)
              : b.dir == Ea && Va(d, d.g, Sa, p - b.ha, 0);
          }
          m = 0;
          for (l = 1; l <= q; l++) e[l].x = -1;
        case 2:
          if (
            (m++, a.p.o >= Wb && x("Pass " + m + ""), (r = s), (p = 0), 1 < m)
          ) {
            null == f && (f = sg());
            for (l = 1; l <= q; l++)
              (k = e[l].C),
                (g = d.f[k]),
                (n = tg(f)),
                0 > n && (n = 0),
                (g = Math.abs(e[l].x - g.r)),
                0.5 < g + n && (e[l].x = 1 - e[l].x);
            n = 4;
            break;
          }
        case 3:
          for (l = h = 1; l <= q; l++)
            (g = d.f[e[l].C]),
              (g = 0.5 > g.r ? 0 : 1),
              e[l].x != g && ((h = 0), (e[l].x = g));
          if (h) {
            for (l = 1; l <= q; l++)
              (g = d.f[e[l].C]), (e[l].Td = Math.abs(g.r - e[l].x));
            ma(e, q, function (a, b) {
              return a.Td > b.Td ? -1 : a.Td < b.Td ? 1 : 0;
            });
            for (l = 1; l <= q && !((5 <= l && 0.35 > e[l].Td) || 10 <= l); l++)
              e[l].x = 1 - e[l].x;
          }
        case 4:
          if (2147483647 > a.p.sb && a.p.sb - 1 <= 1e3 * la(a.hc)) {
            n = 5;
            break;
          }
          d.dir = za;
          d.ha = 0;
          for (k = 1; k <= c; k++) d.f[k].u = 0;
          for (l = 1; l <= q; l++)
            (k = e[l].C),
              0 == e[l].x ? (d.f[k].u = 1) : ((d.f[k].u = -1), (d.ha += 1));
          h = new kc();
          a.p.o <= Lb
            ? (h.o = a.p.o)
            : a.p.o <= Wb && ((h.o = fc), (h.fb = 1e4));
          l = sc(d, h);
          if (0 != l) {
            a.p.o >= Lb && x("Warning: glp_simplex returned " + l + "");
            n = 5;
            break;
          }
          l = xc(d);
          if (l != vc) {
            a.p.o >= Lb && x("Warning: glp_get_status returned " + l + "");
            n = 5;
            break;
          }
          a.p.o >= mc && x("delta = " + d.aa + "");
          k = 0.3 * a.p.Ub;
          for (
            l = 1;
            l <= q && !((g = d.f[e[l].C]), k < g.r && g.r < 1 - k);
            l++
          );
          if (l > q) {
            g = new Float64Array(1 + c);
            for (k = 1; k <= c; k++)
              (g[k] = d.f[k].r),
                b.f[k].kind == Fc && (g[k] = Math.floor(g[k] + 0.5));
            d.ha = b.ha;
            d.dir = b.dir;
            for (l = 1; l <= q; l++)
              (d.f[e[l].C].c = g[e[l].C]),
                (d.f[e[l].C].d = g[e[l].C]),
                (d.f[e[l].C].type = B);
            for (k = 1; k <= c; k++) d.f[k].u = b.f[k].u;
            l = sc(d, h);
            if (0 != l) {
              a.p.o >= Lb && x("Warning: glp_simplex returned " + l + "");
              n = 5;
              break;
            }
            l = xc(d);
            if (l != vc) {
              a.p.o >= Lb && x("Warning: glp_get_status returned " + l + "");
              n = 5;
              break;
            }
            for (k = 1; k <= c; k++) b.f[k].kind != Fc && (g[k] = d.f[k].r);
            l = Kd(a, g);
            if (0 == l) {
              n = vf(a, a.N.bound) ? 1 : 5;
              break;
            }
          }
          r == s || d.aa <= r - 1e-6 * (1 + r) ? ((p = 0), (r = d.aa)) : p++;
          if (3 > p) {
            n = 3;
            break;
          }
          5 > m && (n = 2);
      }
      if (null == n) break;
      l = n;
    }
  }
  function $f(a) {
    function b(a, b, c) {
      var d,
        e = 0,
        f = 0,
        g = 0;
      for (d = a.k; null != d; d = d.e) (c[d.C] = d.j), (f += d.j * d.j);
      for (d = b.k; null != d; d = d.e) (e += c[d.C] * d.j), (g += d.j * d.j);
      for (d = a.k; null != d; d = d.e) c[d.C] = 0;
      a = Math.sqrt(f) * Math.sqrt(g);
      4.930380657631324e-32 > a && (a = 2.220446049250313e-16);
      return e / a;
    }
    var c, d, e, f, g, h, k, l, p, m;
    c = a.Bd;
    f = Array(1 + c.size);
    l = new Int32Array(1 + a.i);
    p = new Float64Array(1 + a.i);
    m = new Float64Array(1 + a.i);
    g = 0;
    for (d = c.head; null != d; d = d.e) g++, (f[g].Re = d), (f[g].ba = 0);
    for (g = 1; g <= c.size; g++) {
      var q = null,
        r = null;
      d = f[g].Re;
      h = k = 0;
      for (e = d.k; null != e; e = e.e)
        k++, (l[k] = e.C), (p[k] = e.j), (h += e.j * e.j);
      4.930380657631324e-32 > h && (h = 2.220446049250313e-16);
      k = Dd(a.A, k, l, p);
      d = Gd(a.A, k, l, p, d.type, d.cg, function (a, b, c, d, e, f) {
        q = e;
        r = f;
      });
      0 == d
        ? ((f[g].Xc = Math.abs(q) / Math.sqrt(h)),
          a.A.dir == za
            ? (0 > r && (r = 0), (f[g].Ab = +r))
            : (0 < r && (r = 0), (f[g].Ab = -r)))
        : 1 == d
        ? (f[g].Xc = f[g].Ab = 0)
        : 2 == d && ((f[g].Xc = 1), (f[g].Ab = s));
      0.01 > f[g].Ab && (f[g].Ab = 0);
    }
    ma(f, c.size, function (a, b) {
      if (0 == a.Ab && 0 == b.Ab) {
        if (a.Xc > b.Xc) return -1;
        if (a.Xc < b.Xc) return 1;
      } else {
        if (a.Ab > b.Ab) return -1;
        if (a.Ab < b.Ab) return 1;
      }
      return 0;
    });
    h = 0 == a.N.La ? 90 : 10;
    h > c.size && (h = c.size);
    for (g = 1; g <= h; g++)
      if (!(0.01 > f[g].Ab && 0.01 > f[g].Xc)) {
        for (c = 1; c < g && !(f[c].ba && 0.9 < b(f[g].Re, f[c].Re, m)); c++);
        if (!(c < g)) {
          d = f[g].Re;
          f[g].ba = 1;
          c = La(a.A, 1);
          null != d.name && Pa(a.A, c, d.name);
          a.A.n[c].qc = d.qc;
          k = 0;
          for (e = d.k; null != e; e = e.e) k++, (l[k] = e.C), (p[k] = e.j);
          Ya(a.A, c, k, l, p);
          Va(a.A, c, d.type, d.cg, d.cg);
        }
      }
  }
  function Rf(a) {
    function b(a) {
      var b, c;
      b = 0;
      c = s;
      for (a = a.head; null != a; a = a.e)
        c > a.R.Zc && ((b = a.s), (c = a.R.Zc));
      return b;
    }
    function c(a) {
      var b, c, d, e, p;
      b = a.ya[1].rb;
      e = (a.A.ta - b.bound) / b.Zc;
      c = 0;
      d = s;
      for (b = a.head; null != b; b = b.e)
        (p = b.R.bound + e * b.R.Zc),
          a.A.dir == Ea && (p = -p),
          d > p && ((c = b.s), (d = p));
      return c;
    }
    function d(a) {
      var b,
        c = null,
        d,
        e;
      switch (a.A.dir) {
        case za:
          d = +s;
          for (b = a.head; null != b; b = b.e) d > b.bound && (d = b.bound);
          e = 0.001 * (1 + Math.abs(d));
          for (b = a.head; null != b; b = b.e)
            b.bound <= d + e && (null == c || c.R.Zc > b.R.Zc) && (c = b);
          break;
        case Ea:
          d = -s;
          for (b = a.head; null != b; b = b.e) d < b.bound && (d = b.bound);
          e = 0.001 * (1 + Math.abs(d));
          for (b = a.head; null != b; b = b.e)
            b.bound >= d - e && (null == c || c.rc < b.rc) && (c = b);
      }
      return c.s;
    }
    var e;
    a.p.kc == dd
      ? (e = a.Xa.s)
      : a.p.kc == ed
      ? (e = a.head.s)
      : a.p.kc == fd
      ? (e = d(a))
      : a.p.kc == gd && (e = a.A.za == Aa ? b(a) : c(a));
    return e;
  }
  var pa = (exports.GLP_MAJOR_VERSION = 4),
    qa = (exports.GLP_MINOR_VERSION = 49),
    za = (exports.GLP_MIN = 1),
    Ea = (exports.GLP_MAX = 2),
    Ma = (exports.GLP_CV = 1),
    Fc = (exports.GLP_IV = 2),
    Gc = (exports.GLP_BV = 3),
    Ka = (exports.GLP_FR = 1),
    Sa = (exports.GLP_LO = 2),
    Ta = (exports.GLP_UP = 3),
    I = (exports.GLP_DB = 4),
    B = (exports.GLP_FX = 5),
    A = (exports.GLP_BS = 1),
    G = (exports.GLP_NL = 2),
    Ua = (exports.GLP_NU = 3),
    Ra = (exports.GLP_NF = 4),
    Na = (exports.GLP_NS = 5),
    Uc = (exports.GLP_SF_GM = 1),
    Vc = (exports.GLP_SF_EQ = 16),
    Wc = (exports.GLP_SF_2N = 32),
    Xc = (exports.GLP_SF_SKIP = 64),
    hc = (exports.GLP_SF_AUTO = 128),
    Zb = (exports.GLP_SOL = 1),
    le = (exports.GLP_IPT = 2),
    Sc = (exports.GLP_MIP = 3),
    Aa = (exports.GLP_UNDEF = 1),
    dc = (exports.GLP_FEAS = 2),
    Ad = (exports.GLP_INFEAS = 3),
    jc = (exports.GLP_NOFEAS = 4),
    vc = (exports.GLP_OPT = 5),
    wc = (exports.GLP_UNBND = 6),
    md = (exports.GLP_BF_FT = 1),
    rd = (exports.GLP_BF_BG = 2),
    sd = (exports.GLP_BF_GR = 3),
    lc = (exports.GLP_MSG_OFF = 0),
    Lb = (exports.GLP_MSG_ERR = 1),
    fc = (exports.GLP_MSG_ON = 2),
    Wb = (exports.GLP_MSG_ALL = 3),
    mc = (exports.GLP_MSG_DBG = 4),
    Ob = (exports.GLP_PRIMAL = 1),
    Qb = (exports.GLP_DUALP = 2),
    Tb = (exports.GLP_DUAL = 3),
    nc = (exports.GLP_PT_STD = 17),
    oc = (exports.GLP_PT_PSE = 34),
    pc = (exports.GLP_RT_STD = 17),
    qc = (exports.GLP_RT_HAR = 34);
  exports.GLP_ORD_NONE = 0;
  exports.GLP_ORD_QMD = 1;
  exports.GLP_ORD_AMD = 2;
  exports.GLP_ORD_SYMAMD = 3;
  var Zc = (exports.GLP_BR_FFV = 1),
    $c = (exports.GLP_BR_LFV = 2),
    ad = (exports.GLP_BR_MFV = 3),
    bd = (exports.GLP_BR_DTH = 4),
    cd = (exports.GLP_BR_PCH = 5),
    dd = (exports.GLP_BT_DFS = 1),
    ed = (exports.GLP_BT_BFS = 2),
    fd = (exports.GLP_BT_BLB = 3),
    gd = (exports.GLP_BT_BPH = 4),
    hd = (exports.GLP_PP_NONE = 0),
    id = (exports.GLP_PP_ROOT = 1),
    jd = (exports.GLP_PP_ALL = 2);
  exports.GLP_RF_REG = 0;
  var Ha = (exports.GLP_RF_LAZY = 1),
    Ja = (exports.GLP_RF_CUT = 2),
    Df = (exports.GLP_RF_GMI = 1),
    Ef = (exports.GLP_RF_MIR = 2),
    Ff = (exports.GLP_RF_COV = 3),
    Gf = (exports.GLP_RF_CLQ = 4),
    bb = (exports.GLP_ON = 1),
    cb = (exports.GLP_OFF = 0),
    Ga = (exports.GLP_IROWGEN = 1),
    Xf = (exports.GLP_IBINGO = 2),
    Yf = (exports.GLP_IHEUR = 3),
    Ia = (exports.GLP_ICUTGEN = 4),
    bg = (exports.GLP_IBRANCH = 5),
    Qf = (exports.GLP_ISELECT = 6),
    Uf = (exports.GLP_IPREPRO = 7),
    Af = (exports.GLP_NO_BRNCH = 0),
    Bf = (exports.GLP_DN_BRNCH = 1),
    Cf = (exports.GLP_UP_BRNCH = 2),
    Kb = (exports.GLP_EBADB = 1),
    Mb = (exports.GLP_ESING = 2),
    Nb = (exports.GLP_ECOND = 3),
    rc = (exports.GLP_EBOUND = 4),
    Sb = (exports.GLP_EFAIL = 5),
    Vf = (exports.GLP_EOBJLL = 6),
    Wf = (exports.GLP_EOBJUL = 7),
    rg = (exports.GLP_EITLIM = 8),
    Qc = (exports.GLP_ETMLIM = 9),
    ac = (exports.GLP_ENOPFS = 10),
    bc = (exports.GLP_ENODFS = 11),
    Lc = (exports.GLP_EROOT = 12),
    Rc = (exports.GLP_ESTOP = 13),
    Pc = (exports.GLP_EMIPGAP = 14);
  exports.GLP_ENOFEAS = 15;
  exports.GLP_ENOCVG = 16;
  exports.GLP_EINSTAB = 17;
  exports.GLP_EDATA = 18;
  exports.GLP_ERANGE = 19;
  exports.GLP_KKT_PE = 1;
  exports.GLP_KKT_PB = 2;
  exports.GLP_KKT_DE = 3;
  exports.GLP_KKT_DB = 4;
  exports.GLP_KKT_CS = 5;
  exports.GLP_MPS_DECK = 1;
  exports.GLP_MPS_FILE = 2;
  exports.GLP_ASN_MIN = 1;
  exports.GLP_ASN_MAX = 2;
  exports.GLP_ASN_MMP = 3;
  function ug(a) {
    var b = Math.floor(Math.log(a) / Math.log(2)) + 1;
    return Math.pow(2, 0.75 >= a / Math.pow(2, b) ? b - 1 : b);
  }
  function vg(a, b) {
    var c = Number(a);
    if (isNaN(c)) return 2;
    switch (c) {
      case Number.POSITIVE_INFINITY:
      case Number.NEGATIVE_INFINITY:
        return 1;
      default:
        return b(c), 0;
    }
  }
  function wg(a, b) {
    var c = Number(a);
    if (isNaN(c)) return 2;
    switch (c) {
      case Number.POSITIVE_INFINITY:
      case Number.NEGATIVE_INFINITY:
        return 1;
      default:
        return 0 == c % 1 ? (b(c), 0) : 2;
    }
  }
  function xg(a, b, c) {
    var d, e;
    if (!(1 <= a && 31 >= a && 1 <= b && 12 >= b && 1 <= c && 4e3 >= c))
      return -1;
    3 <= b ? (b -= 3) : ((b += 9), c--);
    d = (c / 100) | 0;
    c = (((146097 * d) / 4) | 0) + (((1461 * (c - 100 * d)) / 4) | 0);
    c += ((153 * b + 2) / 5) | 0;
    c += a + 1721119;
    yg(c, function (a) {
      e = a;
    });
    a != e && (c = -1);
    return c;
  }
  function yg(a, b) {
    var c, d, e;
    1721426 <= a &&
      3182395 >= a &&
      ((a -= 1721119),
      (e = ((4 * a - 1) / 146097) | 0),
      (c = (((4 * a - 1) % 146097) / 4) | 0),
      (a = ((4 * c + 3) / 1461) | 0),
      (c = ((((4 * c + 3) % 1461) + 4) / 4) | 0),
      (d = ((5 * c - 3) / 153) | 0),
      (c = (5 * c - 3) % 153),
      (c = ((c + 5) / 5) | 0),
      (e = 100 * e + a),
      9 >= d ? (d += 3) : ((d -= 9), e++),
      b(c, d, e));
  }
  var Ee = 1;
  LPF_ECOND = 2;
  LPF_ELIMIT = 3;
  var we = 0;
  function Me(a, b, c, d) {
    var e = a.i,
      f = a.Md,
      g = a.Ld,
      h = a.Wb;
    a = a.Xb;
    var k, l, p, m;
    for (k = 1; k <= e; k++) {
      m = 0;
      l = f[k];
      for (p = l + g[k]; l < p; l++) m += a[l] * d[h[l]];
      b[k + c] += -1 * m;
    }
  }
  function Je(a, b, c, d) {
    var e = a.i,
      f = a.Od,
      g = a.Nd,
      h = a.Wb;
    a = a.Xb;
    var k, l, p, m;
    for (k = 1; k <= e; k++) {
      m = 0;
      l = f[k];
      for (p = l + g[k]; l < p; l++) m += a[l] * d[h[l]];
      b[k + c] += -1 * m;
    }
  }
  if (we)
    var Le = function (a, b, c, d) {
      var e = a.g;
      a = a.Cf;
      var f,
        g,
        h = 0,
        k,
        l,
        p;
      for (f = 1; f <= e; f++) {
        k = 0;
        for (g = p = 1; g <= e; g++)
          (l = b ? a[e * (g - 1) + f] * c[g] : a[e * (f - 1) + g] * c[g]),
            p < Math.abs(l) && (p = Math.abs(l)),
            (k += l);
        g = Math.abs(k - d[f]) / p;
        h < g && (h = g);
      }
      1e-8 < h &&
        x(
          (b ? "lpf_btran" : "lpf_ftran") +
            ": dmax = " +
            h +
            "; relative error too large"
        );
    };
  exports.LPX_LP = 100;
  exports.LPX_MIP = 101;
  var gf = (exports.LPX_FR = 110),
    jf = (exports.LPX_LO = 111),
    lf = (exports.LPX_UP = 112),
    nf = (exports.LPX_DB = 113),
    cf = (exports.LPX_FX = 114);
  exports.LPX_MIN = 120;
  exports.LPX_MAX = 121;
  exports.LPX_P_UNDEF = 132;
  exports.LPX_P_FEAS = 133;
  exports.LPX_P_INFEAS = 134;
  exports.LPX_P_NOFEAS = 135;
  exports.LPX_D_UNDEF = 136;
  exports.LPX_D_FEAS = 137;
  exports.LPX_D_INFEAS = 138;
  exports.LPX_D_NOFEAS = 139;
  var ff = (exports.LPX_BS = 140),
    kf = (exports.LPX_NL = 141),
    mf = (exports.LPX_NU = 142),
    hf = (exports.LPX_NF = 143),
    of = (exports.LPX_NS = 144);
  exports.LPX_T_UNDEF = 150;
  exports.LPX_T_OPT = 151;
  var Mf = (exports.LPX_CV = 160),
    Nf = (exports.LPX_IV = 161);
  exports.LPX_I_UNDEF = 170;
  exports.LPX_I_OPT = 171;
  exports.LPX_I_FEAS = 172;
  exports.LPX_I_NOFEAS = 173;
  exports.LPX_OPT = 180;
  exports.LPX_FEAS = 181;
  exports.LPX_INFEAS = 182;
  exports.LPX_NOFEAS = 183;
  exports.LPX_UNBND = 184;
  exports.LPX_UNDEF = 185;
  exports.LPX_E_OK = 200;
  exports.LPX_E_EMPTY = 201;
  exports.LPX_E_BADB = 202;
  exports.LPX_E_INFEAS = 203;
  exports.LPX_E_FAULT = 204;
  exports.LPX_E_OBJLL = 205;
  exports.LPX_E_OBJUL = 206;
  exports.LPX_E_ITLIM = 207;
  exports.LPX_E_TMLIM = 208;
  exports.LPX_E_NOFEAS = 209;
  exports.LPX_E_INSTAB = 210;
  exports.LPX_E_SING = 211;
  exports.LPX_E_NOCONV = 212;
  exports.LPX_E_NOPFS = 213;
  exports.LPX_E_NODFS = 214;
  exports.LPX_E_MIPGAP = 215;
  var zg = (exports.LPX_K_MSGLEV = 300),
    Ag = (exports.LPX_K_SCALE = 301),
    Bg = (exports.LPX_K_DUAL = 302),
    Cg = (exports.LPX_K_PRICE = 303);
  exports.LPX_K_RELAX = 304;
  exports.LPX_K_TOLBND = 305;
  exports.LPX_K_TOLDJ = 306;
  exports.LPX_K_TOLPIV = 307;
  var Dg = (exports.LPX_K_ROUND = 308);
  exports.LPX_K_OBJLL = 309;
  exports.LPX_K_OBJUL = 310;
  var Eg = (exports.LPX_K_ITLIM = 311),
    Fg = (exports.LPX_K_ITCNT = 312);
  exports.LPX_K_TMLIM = 313;
  var Gg = (exports.LPX_K_OUTFRQ = 314);
  exports.LPX_K_OUTDLY = 315;
  var Hg = (exports.LPX_K_BRANCH = 316),
    Ig = (exports.LPX_K_BTRACK = 317);
  exports.LPX_K_TOLINT = 318;
  exports.LPX_K_TOLOBJ = 319;
  var Jg = (exports.LPX_K_MPSINFO = 320),
    Kg = (exports.LPX_K_MPSOBJ = 321),
    Lg = (exports.LPX_K_MPSORIG = 322),
    Mg = (exports.LPX_K_MPSWIDE = 323),
    Ng = (exports.LPX_K_MPSFREE = 324),
    Og = (exports.LPX_K_MPSSKIP = 325),
    Pg = (exports.LPX_K_LPTORIG = 326),
    Qg = (exports.LPX_K_PRESOL = 327),
    Rg = (exports.LPX_K_BINARIZE = 328),
    Sg = (exports.LPX_K_USECUTS = 329),
    Tg = (exports.LPX_K_BFTYPE = 330);
  exports.LPX_K_MIPGAP = 331;
  exports.LPX_C_COVER = 1;
  exports.LPX_C_CLIQUE = 2;
  exports.LPX_C_GOMORY = 4;
  exports.LPX_C_MIR = 8;
  exports.LPX_C_ALL = 255;
  function Kf(a, b) {
    var c = pb(a, b);
    c == -s && (c = 0);
    return c;
  }
  function Jf(a, b) {
    var c = qb(a, b);
    c == +s && (c = 0);
    return c;
  }
  function df(a, b, c) {
    c(ob(a, b) - Ka + gf, Kf(a, b), Jf(a, b));
  }
  function Lf(a, b) {
    var c = sb(a, b);
    c == -s && (c = 0);
    return c;
  }
  function Of(a, b) {
    var c = tb(a, b);
    c == +s && (c = 0);
    return c;
  }
  function bf(a, b, c) {
    c(rb(a, b) - Ka + gf, Lf(a, b), Of(a, b));
  }
  function ef(a) {
    var b = zg,
      c;
    null == a.ke &&
      ((a.ke = {}),
      (c = a.ke),
      (c.o = 3),
      (c.scale = 1),
      (c.J = 0),
      (c.mh = 1),
      (c.yh = 0.07),
      (c.Gb = 1e-7),
      (c.tb = 1e-7),
      (c.xe = 1e-9),
      (c.round = 0),
      (c.hf = -s),
      (c.jf = +s),
      (c.oc = -1),
      (c.sb = -1),
      (c.bc = 200),
      (c.fb = 0),
      (c.Og = 2),
      (c.Pg = 3),
      (c.Ub = 1e-5),
      (c.we = 1e-7),
      (c.Xg = 1),
      (c.Yg = 2),
      (c.Zg = 0),
      (c.ah = 1),
      (c.Wg = 0),
      (c.$g = 0),
      (c.Vg = 0),
      (c.lh = 0),
      (c.sd = 0),
      (c.wh = 0),
      (c.ce = 0));
    c = a.ke;
    var d = 0;
    switch (b) {
      case zg:
        d = c.o;
        break;
      case Ag:
        d = c.scale;
        break;
      case Bg:
        d = c.J;
        break;
      case Cg:
        d = c.mh;
        break;
      case Dg:
        d = c.round;
        break;
      case Eg:
        d = c.oc;
        break;
      case Fg:
        d = a.$;
        break;
      case Gg:
        d = c.bc;
        break;
      case Hg:
        d = c.Og;
        break;
      case Ig:
        d = c.Pg;
        break;
      case Jg:
        d = c.Xg;
        break;
      case Kg:
        d = c.Yg;
        break;
      case Lg:
        d = c.Zg;
        break;
      case Mg:
        d = c.ah;
        break;
      case Ng:
        d = c.Wg;
        break;
      case Og:
        d = c.$g;
        break;
      case Pg:
        d = c.Vg;
        break;
      case Qg:
        d = c.lh;
        break;
      case Rg:
        d = c.sd;
        break;
      case Sg:
        d = c.wh;
        break;
      case Tg:
        b = {};
        eb(a, b);
        switch (b.type) {
          case md:
            d = 1;
            break;
          case rd:
            d = 2;
            break;
          case sd:
            d = 3;
        }
        break;
      default:
        w("lpx_get_int_parm: parm = " + b + "; invalid parameter");
    }
    return d;
  }
  var ye = 1,
    Ae = 2;
  function ve() {
    var a = {};
    a.K = a.i = 0;
    a.valid = 0;
    a.Rf = a.Ye = null;
    a.Yd = a.Xd = null;
    a.Fc = a.Ec = a.pd = null;
    a.Af = null;
    a.Dc = a.Cc = a.Rc = null;
    a.kb = a.vb = null;
    a.sf = a.me = null;
    a.Ga = 0;
    a.Fa = a.Ma = 0;
    a.wb = null;
    a.xb = null;
    a.ld = a.Sb = 0;
    a.Hd = a.md = null;
    a.zf = null;
    a.vf = a.dg = a.wf = null;
    a.Oe = a.Qe = a.Pe = null;
    a.ba = null;
    a.ze = null;
    a.Va = 0;
    a.cc = 0.1;
    a.xc = 4;
    a.gc = 1;
    a.Mb = 1e-15;
    a.sc = 1e10;
    a.dh = a.$f = a.Qb = 0;
    a.Cg = a.rd = 0;
    a.pa = 0;
    return a;
  }
  function af(a) {
    var b = a.i,
      c = a.Fc,
      d = a.Ec,
      e = a.pd,
      f = a.Dc,
      g = a.Cc,
      h = a.Rc,
      k = a.wb,
      l = a.xb,
      p = a.md,
      m = 1,
      q,
      r;
    for (r = a.ld; 0 != r; r = p[r])
      if (r <= b) {
        q = r;
        if (c[q] != m) break;
        e[q] = d[q];
        m += e[q];
      } else {
        q = r - b;
        if (f[q] != m) break;
        h[q] = g[q];
        m += h[q];
      }
    for (; 0 != r; r = p[r])
      r <= b
        ? ((q = r),
          ga(k, m, k, c[q], d[q]),
          ga(l, m, l, c[q], d[q]),
          (c[q] = m),
          (e[q] = d[q]),
          (m += e[q]))
        : ((q = r - b),
          ga(k, m, k, f[q], g[q]),
          ga(l, m, l, f[q], g[q]),
          (f[q] = m),
          (h[q] = g[q]),
          (m += h[q]));
    a.Fa = m;
  }
  function Ze(a, b, c) {
    var d = a.i,
      e = a.Fc,
      f = a.Ec,
      g = a.pd,
      h = a.Rc,
      k = a.wb,
      l = a.xb,
      p = a.Hd,
      m = a.md,
      q = 0,
      r;
    if (a.Ma - a.Fa < c && (af(a), a.Ma - a.Fa < c)) return 1;
    r = g[b];
    ga(k, a.Fa, k, e[b], f[b]);
    ga(l, a.Fa, l, e[b], f[b]);
    e[b] = a.Fa;
    g[b] = c;
    a.Fa += c;
    0 == p[b]
      ? (a.ld = m[b])
      : ((c = p[b]), c <= d ? (g[c] += r) : (h[c - d] += r), (m[p[b]] = m[b]));
    0 == m[b] ? (a.Sb = p[b]) : (p[m[b]] = p[b]);
    p[b] = a.Sb;
    m[b] = 0;
    0 == p[b] ? (a.ld = b) : (m[p[b]] = b);
    a.Sb = b;
    return q;
  }
  function $e(a, b, c) {
    var d = a.i,
      e = a.pd,
      f = a.Dc,
      g = a.Cc,
      h = a.Rc,
      k = a.wb,
      l = a.xb,
      p = a.Hd,
      m = a.md,
      q = 0,
      r;
    if (a.Ma - a.Fa < c && (af(a), a.Ma - a.Fa < c)) return 1;
    r = h[b];
    ga(k, a.Fa, k, f[b], g[b]);
    ga(l, a.Fa, l, f[b], g[b]);
    f[b] = a.Fa;
    h[b] = c;
    a.Fa += c;
    b = d + b;
    0 == p[b]
      ? (a.ld = m[b])
      : ((c = p[b]), c <= d ? (e[c] += r) : (h[c - d] += r), (m[p[b]] = m[b]));
    0 == m[b] ? (a.Sb = p[b]) : (p[m[b]] = p[b]);
    p[b] = a.Sb;
    m[b] = 0;
    0 == p[b] ? (a.ld = b) : (m[p[b]] = b);
    a.Sb = b;
    return q;
  }
  function Ug(a, b) {
    var c = a.K;
    a.i = b;
    b <= c ||
      ((a.K = c = b + 100),
      (a.Rf = new Int32Array(1 + c)),
      (a.Ye = new Int32Array(1 + c)),
      (a.Yd = new Int32Array(1 + c)),
      (a.Xd = new Int32Array(1 + c)),
      (a.Fc = new Int32Array(1 + c)),
      (a.Ec = new Int32Array(1 + c)),
      (a.pd = new Int32Array(1 + c)),
      (a.Af = new Float64Array(1 + c)),
      (a.Dc = new Int32Array(1 + c)),
      (a.Cc = new Int32Array(1 + c)),
      (a.Rc = new Int32Array(1 + c)),
      (a.kb = new Int32Array(1 + c)),
      (a.vb = new Int32Array(1 + c)),
      (a.sf = new Int32Array(1 + c)),
      (a.me = new Int32Array(1 + c)),
      (a.Hd = new Int32Array(1 + c + c)),
      (a.md = new Int32Array(1 + c + c)),
      (a.zf = new Float64Array(1 + c)),
      (a.vf = new Int32Array(1 + c)),
      (a.dg = new Int32Array(1 + c)),
      (a.wf = new Int32Array(1 + c)),
      (a.Oe = new Int32Array(1 + c)),
      (a.Qe = new Int32Array(1 + c)),
      (a.Pe = new Int32Array(1 + c)),
      (a.ba = new Int32Array(1 + c)),
      (a.ze = new Float64Array(1 + c)));
  }
  function Vg(a, b, c) {
    var d = a.i,
      e = a.Yd,
      f = a.Xd,
      g = a.Fc,
      h = a.Ec,
      k = a.pd,
      l = a.Dc,
      p = a.Cc,
      m = a.Rc,
      q = a.kb,
      r = a.vb,
      n = a.sf,
      t = a.me,
      y = a.wb,
      E = a.xb,
      C = a.Hd,
      D = a.md,
      H = a.zf,
      R = a.vf,
      V = a.dg,
      O = a.wf,
      Q = a.Oe,
      F = a.Qe,
      W = a.Pe,
      X = a.ba,
      ca = a.ze,
      ka = 0,
      P,
      u,
      z,
      L,
      v,
      S,
      M;
    L = 1;
    v = a.Ga + 1;
    for (u = 1; u <= d; u++) (e[u] = v), (f[u] = 0);
    for (P = 1; P <= d; P++) (h[P] = k[P] = 0), (X[P] = 0);
    f = e = 0;
    for (u = 1; u <= d; u++) {
      var ba = q,
        J = ca;
      z = b(c, u, ba, J);
      (0 <= z && z <= d) ||
        w(
          "luf_factorize: j = " + u + "; len = " + z + "; invalid column length"
        );
      if (v - L < z) return (ka = 1);
      l[u] = L;
      p[u] = m[u] = z;
      e += z;
      for (S = 1; S <= z; S++)
        (P = ba[S]),
          (M = J[S]),
          (1 <= P && P <= d) ||
            w("luf_factorize: i = " + P + "; j = " + u + "; invalid row index"),
          X[P] &&
            w(
              "luf_factorize: i = " +
                P +
                "; j = " +
                u +
                "; duplicate element not allowed"
            ),
          0 == M &&
            w(
              "luf_factorize: i = " +
                P +
                "; j = " +
                u +
                "; zero element not allowed"
            ),
          (y[L] = P),
          (E[L] = M),
          L++,
          0 > M && (M = -M),
          f < M && (f = M),
          (X[P] = 1),
          k[P]++;
      for (S = 1; S <= z; S++) X[ba[S]] = 0;
    }
    for (P = 1; P <= d; P++) {
      z = k[P];
      if (v - L < z) return (ka = 1);
      g[P] = L;
      L += z;
    }
    for (u = 1; u <= d; u++)
      for (P = l[u], b = P + p[u] - 1, k = P; k <= b; k++)
        (P = y[k]),
          (M = E[k]),
          (c = g[P] + h[P]),
          (y[c] = u),
          (E[c] = M),
          h[P]++;
    for (k = 1; k <= d; k++) q[k] = r[k] = n[k] = t[k] = k;
    a.Fa = L;
    a.Ma = v;
    a.ld = d + 1;
    a.Sb = d;
    for (P = 1; P <= d; P++) (C[P] = P - 1), (D[P] = P + 1);
    C[1] = d + d;
    D[d] = 0;
    for (u = 1; u <= d; u++) (C[d + u] = d + u - 1), (D[d + u] = d + u + 1);
    C[d + 1] = 0;
    for (k = D[d + d] = 1; k <= d; k++) (X[k] = 0), (ca[k] = 0);
    a.dh = e;
    a.$f = 0;
    a.Qb = e;
    a.Cg = f;
    a.rd = f;
    a.pa = -1;
    for (P = 1; P <= d; P++) H[P] = -1;
    for (z = 0; z <= d; z++) R[z] = 0;
    for (P = 1; P <= d; P++)
      (z = h[P]),
        (V[P] = 0),
        (O[P] = R[z]),
        0 != O[P] && (V[O[P]] = P),
        (R[z] = P);
    for (z = 0; z <= d; z++) Q[z] = 0;
    for (u = 1; u <= d; u++)
      (z = p[u]),
        (F[u] = 0),
        (W[u] = Q[z]),
        0 != W[u] && (F[W[u]] = u),
        (Q[z] = u);
    return ka;
  }
  function Wg(a, b) {
    function c() {
      b(D, H);
      return 0 == D;
    }
    var d = a.i,
      e = a.Fc,
      f = a.Ec,
      g = a.Dc,
      h = a.Cc,
      k = a.wb,
      l = a.xb,
      p = a.zf,
      m = a.vf,
      q = a.wf,
      r = a.Oe,
      n = a.Qe,
      t = a.Pe,
      y = a.cc,
      E = a.xc,
      C = a.gc,
      D,
      H,
      R,
      V,
      O,
      Q,
      F,
      W,
      X,
      ca,
      ka,
      P,
      u,
      z,
      L,
      v,
      S,
      M;
    D = H = 0;
    v = s;
    ka = 0;
    W = r[1];
    if (0 != W) return (D = k[g[W]]), (H = W), c();
    V = m[1];
    if (0 != V) return (D = V), (H = k[e[V]]), c();
    for (R = 2; R <= d; R++) {
      for (W = r[R]; 0 != W; W = P) {
        V = g[W];
        X = V + h[W] - 1;
        P = t[W];
        u = z = 0;
        L = 2147483647;
        for (ca = V; ca <= X; ca++)
          if (((V = k[ca]), (O = e[V]), (Q = O + f[V] - 1), !(f[V] >= L))) {
            S = p[V];
            if (0 > S) {
              for (F = O; F <= Q; F++)
                (M = l[F]), 0 > M && (M = -M), S < M && (S = M);
              p[V] = S;
            }
            for (F = e[V]; k[F] != W; F++);
            M = l[F];
            0 > M && (M = -M);
            if (!(M < y * S) && ((u = V), (z = W), (L = f[V]), L <= R))
              return (D = u), (H = z), c();
          }
        if (0 != u) {
          if (
            (ka++,
            (W = (L - 1) * (R - 1)),
            W < v && ((D = u), (H = z), (v = W)),
            ka == E)
          )
            return c();
        } else
          C &&
            (0 == n[W] ? (r[R] = t[W]) : (t[n[W]] = t[W]),
            0 != t[W] && (n[t[W]] = n[W]),
            (n[W] = t[W] = W));
      }
      for (V = m[R]; 0 != V; V = q[V]) {
        O = e[V];
        Q = O + f[V] - 1;
        S = p[V];
        if (0 > S) {
          for (F = O; F <= Q; F++)
            (M = l[F]), 0 > M && (M = -M), S < M && (S = M);
          p[V] = S;
        }
        u = z = 0;
        L = 2147483647;
        for (F = O; F <= Q; F++)
          if (
            ((W = k[F]),
            !(h[W] >= L) &&
              ((M = l[F]),
              0 > M && (M = -M),
              !(M < y * S) && ((u = V), (z = W), (L = h[W]), L <= R)))
          )
            return (D = u), (H = z), c();
        if (
          0 != u &&
          (ka++,
          (W = (R - 1) * (L - 1)),
          W < v && ((D = u), (H = z), (v = W)),
          ka == E)
        )
          return c();
      }
    }
    return c();
  }
  function Xg(a, b, c) {
    var d = a.i,
      e = a.Yd,
      f = a.Xd,
      g = a.Fc,
      h = a.Ec,
      k = a.pd,
      l = a.Af,
      p = a.Dc,
      m = a.Cc,
      q = a.Rc,
      r = a.wb,
      n = a.xb,
      t = a.Hd,
      y = a.md,
      E = a.zf,
      C = a.vf,
      D = a.dg,
      H = a.wf,
      R = a.Oe,
      V = a.Qe,
      O = a.Pe,
      Q = a.ba,
      F = a.ze,
      W = a.Mb,
      X = a.Ye,
      ca = 0,
      ka,
      P,
      u,
      z,
      L,
      v,
      S,
      M,
      ba,
      J,
      ea,
      fa;
    0 == D[b] ? (C[h[b]] = H[b]) : (H[D[b]] = H[b]);
    0 != H[b] && (D[H[b]] = D[b]);
    0 == V[c] ? (R[m[c]] = O[c]) : (O[V[c]] = O[c]);
    0 != O[c] && (V[O[c]] = V[c]);
    M = g[b];
    ba = M + h[b] - 1;
    for (P = M; r[P] != c; P++);
    fa = l[b] = n[P];
    r[P] = r[ba];
    n[P] = n[ba];
    h[b]--;
    ba--;
    l = p[c];
    J = l + m[c] - 1;
    for (u = l; r[u] != b; u++);
    r[u] = r[J];
    m[c]--;
    J--;
    for (P = M; P <= ba; P++) {
      z = r[P];
      Q[z] = 1;
      F[z] = n[P];
      0 == V[z] ? (R[m[z]] = O[z]) : (O[V[z]] = O[z]);
      0 != O[z] && (V[O[z]] = V[z]);
      v = p[z];
      for (S = v + m[z] - 1; r[v] != b; v++);
      r[v] = r[S];
      m[z]--;
    }
    for (; l <= J; ) {
      u = r[l];
      0 == D[u] ? (C[h[u]] = H[u]) : (H[D[u]] = H[u]);
      0 != H[u] && (D[H[u]] = D[u]);
      z = g[u];
      ka = z + h[u] - 1;
      for (L = z; r[L] != c; L++);
      ea = n[L] / fa;
      r[L] = r[ka];
      n[L] = n[ka];
      h[u]--;
      ka--;
      r[l] = r[J];
      m[c]--;
      J--;
      P = h[b];
      for (L = z; L <= ka; L++)
        if (((z = r[L]), Q[z]))
          if (
            ((v = n[L] -= ea * F[z]),
            0 > v && (v = -v),
            (Q[z] = 0),
            P--,
            0 == v || v < W)
          ) {
            r[L] = r[ka];
            n[L] = n[ka];
            h[u]--;
            L--;
            ka--;
            v = p[z];
            for (S = v + m[z] - 1; r[v] != u; v++);
            r[v] = r[S];
            m[z]--;
          } else a.rd < v && (a.rd = v);
      if (h[u] + P > k[u]) {
        if (Ze(a, u, h[u] + P)) return (ca = 1);
        M = g[b];
        ba = M + h[b] - 1;
        l = p[c];
        J = l + m[c] - 1;
      }
      ka = 0;
      for (P = M; P <= ba; P++)
        (z = r[P]),
          Q[z]
            ? ((v = S = -ea * F[z]),
              0 > v && (v = -v),
              0 == v ||
                v < W ||
                ((L = g[u] + h[u]),
                (r[L] = z),
                (n[L] = S),
                h[u]++,
                (X[++ka] = z),
                a.rd < v && (a.rd = v)))
            : (Q[z] = 1);
      for (L = 1; L <= ka; L++) {
        z = X[L];
        if (m[z] + 1 > q[z]) {
          if ($e(a, z, m[z] + 10)) return (ca = 1);
          M = g[b];
          ba = M + h[b] - 1;
          l = p[c];
          J = l + m[c] - 1;
        }
        v = p[z] + m[z];
        r[v] = u;
        m[z]++;
      }
      D[u] = 0;
      H[u] = C[h[u]];
      0 != H[u] && (D[H[u]] = u);
      C[h[u]] = u;
      E[u] = -1;
      if (1 > a.Ma - a.Fa) {
        af(a);
        if (1 > a.Ma - a.Fa) return (ca = 1);
        M = g[b];
        ba = M + h[b] - 1;
        l = p[c];
        J = l + m[c] - 1;
      }
      a.Ma--;
      r[a.Ma] = u;
      n[a.Ma] = ea;
      f[b]++;
    }
    q[c] = 0;
    L = d + c;
    0 == t[L] ? (a.ld = y[L]) : (y[t[L]] = y[L]);
    0 == y[L] ? (a.Sb = t[L]) : (t[y[L]] = t[L]);
    e[b] = a.Ma;
    for (P = M; P <= ba; P++)
      if (
        ((z = r[P]),
        (Q[z] = 0),
        (F[z] = 0),
        1 == m[z] || V[z] != z || O[z] != z)
      )
        (V[z] = 0), (O[z] = R[m[z]]), 0 != O[z] && (V[O[z]] = z), (R[m[z]] = z);
    return ca;
  }
  function Yg(a) {
    var b = a.i,
      c = a.Fc,
      d = a.Ec,
      e = a.Dc,
      f = a.Cc,
      g = a.Rc,
      h = a.wb,
      k = a.xb,
      l = a.Hd,
      p = a.md,
      m,
      q,
      r,
      n;
    n = 0;
    for (m = 1; m <= b; m++) {
      q = c[m];
      for (r = q + d[m] - 1; q <= r; q++) g[h[q]]++;
      n += d[m];
    }
    a.Qb = n;
    if (a.Ma - a.Fa < n) return 1;
    for (n = 1; n <= b; n++) (e[n] = a.Fa), (a.Fa += g[n]);
    for (m = 1; m <= b; m++)
      for (q = c[m], r = q + d[m] - 1; q <= r; q++)
        (n = h[q]), (g = e[n] + f[n]), (h[g] = m), (k[g] = k[q]), f[n]++;
    for (c = b + 1; c <= b + b; c++) (l[c] = c - 1), (p[c] = c + 1);
    l[b + 1] = a.Sb;
    p[a.Sb] = b + 1;
    p[b + b] = 0;
    a.Sb = b + b;
    return 0;
  }
  function Zg(a) {
    var b = a.i,
      c = a.Rf,
      d = a.Ye,
      e = a.Yd,
      f = a.Xd,
      g = a.wb,
      h = a.xb,
      k,
      l,
      p,
      m;
    for (k = 1; k <= b; k++) d[k] = 0;
    k = 0;
    for (l = 1; l <= b; l++) {
      p = e[l];
      for (m = p + f[l] - 1; p <= m; p++) d[g[p]]++;
      k += f[l];
    }
    a.$f = k;
    if (a.Ma - a.Fa < k) return 1;
    for (k = 1; k <= b; k++) (c[k] = a.Ma), (a.Ma -= d[k]);
    for (l = 1; l <= b; l++)
      for (p = e[l], m = p + f[l] - 1; p <= m; p++)
        (k = g[p]), (a = --c[k]), (g[a] = l), (h[a] = h[p]);
    return 0;
  }
  function xe(a, b, c, d) {
    function e() {
      0 < a.Va &&
        ((a.Ga = a.Va),
        (a.wb = new Int32Array(1 + a.Ga)),
        (a.xb = new Float64Array(1 + a.Ga)),
        (a.Va = 0));
      if (Vg(a, c, d)) return (a.Va = a.Ga + a.Ga), !0;
      for (q = 1; q <= b; q++) {
        if (
          Wg(a, function (a, b) {
            r = a;
            n = b;
          })
        )
          return (a.pa = q - 1), (y = ye), !1;
        p = g[r];
        m = h[n];
        t = f[q];
        f[p] = t;
        g[t] = p;
        f[q] = r;
        g[r] = q;
        t = k[q];
        k[m] = t;
        h[t] = m;
        k[q] = n;
        h[n] = q;
        if (Xg(a, r, n)) return (a.Va = a.Ga + a.Ga), !0;
        if (a.rd > l * a.Cg) return (a.pa = q - 1), (y = Ae), !1;
      }
      af(a);
      return Yg(a) || Zg(a) ? ((a.Va = a.Ga + a.Ga), !0) : !1;
    }
    var f,
      g,
      h,
      k,
      l = a.sc,
      p,
      m,
      q,
      r,
      n,
      t,
      y = null;
    1 > b && w("luf_factorize: n = " + b + "; invalid parameter");
    1e8 < b && w("luf_factorize: n = " + b + "; matrix too big");
    a.valid = 0;
    Ug(a, b);
    f = a.kb;
    g = a.vb;
    h = a.sf;
    k = a.me;
    0 == a.Ga && 0 == a.Va && (a.Va = 5 * (b + 10));
    for (; e(); );
    if (null != y) return y;
    a.valid = 1;
    a.pa = b;
    y = 0;
    t = 3 * (b + a.Qb) + 2 * a.$f;
    if (a.Ga < t) for (a.Va = a.Ga; a.Va < t; ) (q = a.Va), (a.Va = q + q);
    return y;
  }
  function Ge(a, b, c) {
    var d = a.i,
      e = a.Rf,
      f = a.Ye,
      g = a.Yd,
      h = a.Xd,
      k = a.kb,
      l = a.wb,
      p = a.xb,
      m;
    a.valid || w("luf_f_solve: LU-factorization is not valid");
    if (b)
      for (; 1 <= d; d--) {
        if (((m = k[d]), (a = c[m]), 0 != a))
          for (b = e[m], m = b + f[m] - 1; b <= m; b++) c[l[b]] -= p[b] * a;
      }
    else
      for (e = 1; e <= d; e++)
        if (((m = k[e]), (a = c[m]), 0 != a))
          for (b = g[m], m = b + h[m] - 1; b <= m; b++) c[l[b]] -= p[b] * a;
  }
  function Ie(a, b, c) {
    var d = a.i,
      e = a.Fc,
      f = a.Ec,
      g = a.Af,
      h = a.Dc,
      k = a.Cc,
      l = a.kb,
      p = a.me,
      m = a.wb,
      q = a.xb,
      r = a.ze,
      n,
      t,
      y;
    a.valid || w("luf_v_solve: LU-factorization is not valid");
    for (a = 1; a <= d; a++) (r[a] = c[a]), (c[a] = 0);
    if (b)
      for (a = 1; a <= d; a++) {
        if (((n = l[a]), (t = p[a]), (b = r[t]), 0 != b))
          for (c[n] = b /= g[n], y = e[n], n = y + f[n] - 1; y <= n; y++)
            r[m[y]] -= q[y] * b;
      }
    else
      for (a = d; 1 <= a; a--)
        if (((n = l[a]), (t = p[a]), (b = r[n]), 0 != b))
          for (c[t] = b /= g[n], y = h[t], n = y + k[t] - 1; y <= n; y++)
            r[m[y]] -= q[y] * b;
  }
  var $g = -1,
    ah = 101,
    bh = 102,
    ch = 103,
    dh = 104,
    eh = 105,
    fh = 106,
    gh = 107,
    hh = 108,
    ih = 109,
    jh = 110,
    kh = 111,
    lh = 112,
    mh = 113,
    nh = 114,
    oh = 115,
    ph = 116,
    qh = 117,
    N = 118,
    rh = 119,
    sh = 120,
    th = 121,
    uh = 122,
    vh = 123,
    T = 124,
    wh = 125,
    xh = 126,
    yh = 127,
    zh = 60,
    Ah = 201,
    Bh = 202,
    Ch = 203,
    Dh = 204,
    Eh = 205,
    Fh = 206,
    Gh = 207,
    Hh = 208,
    Ih = 209,
    Jh = 210,
    Kh = 211,
    Lh = 212,
    Mh = 213,
    Nh = 214,
    Oh = 215,
    Ph = 216,
    Qh = 217,
    Rh = 218,
    Sh = 219,
    Th = 220,
    Uh = 221,
    Vh = 222,
    Wh = 223,
    Xh = 224,
    Yh = 225,
    Zh = 226,
    $h = 227,
    ai = 228,
    bi = 229,
    ci = 230,
    di = 231,
    ei = 232,
    fi = 233,
    gi = 234,
    hi = 235,
    ii = 236,
    ji = 237,
    ki = 238,
    li = 239,
    mi = 240,
    ni = 241,
    oi = 242,
    pi = 243,
    qi = 244,
    ri = 245,
    si = 246,
    ti = 247,
    ui = 248,
    vi = 249,
    wi = 250,
    xi = 251,
    yi = 252,
    zi = 0,
    Ai = 1,
    Bi = 2,
    Ci = 3,
    Di = 4,
    Ei = 5,
    Fi = 301,
    Gi = 302,
    Hi = 303,
    Ii = 304,
    Ji = 305,
    Ki = 306,
    Li = 307,
    Mi = 308,
    Ni = 309,
    Oi = 310,
    Pi = 311,
    Qi = 312,
    Ri = 313,
    Si = 314,
    Ti = 315,
    Ui = 316,
    Vi = 317,
    Wi = 318,
    Xi = 319,
    Yi = 320,
    Zi = 321,
    $i = 322,
    aj = 323,
    bj = 324,
    cj = 325,
    dj = 326,
    ej = 327,
    fj = 328,
    gj = 329,
    hj = 330,
    ij = 331,
    jj = 332,
    kj = 333,
    lj = 334,
    mj = 335,
    nj = 336,
    oj = 337,
    pj = 338,
    qj = 339,
    rj = 340,
    sj = 341,
    tj = 342,
    uj = 343,
    vj = 344,
    wj = 345,
    xj = 346,
    yj = 347,
    zj = 348,
    Aj = 349,
    Bj = 350,
    Cj = 351,
    Dj = 352,
    Ej = 353,
    Fj = 354,
    Gj = 355,
    Hj = 356,
    Ij = 357,
    Jj = 358,
    Kj = 359,
    Lj = 360,
    Mj = 361,
    Nj = 362,
    Oj = 363,
    Pj = 364,
    Qj = 365,
    Rj = 366,
    Sj = 367,
    Tj = 368,
    Uj = 369,
    Vj = 370,
    Xj = 371,
    Yj = 372,
    Zj = 373,
    ak = 374,
    bk = 375,
    ck = 376,
    dk = 377,
    ek = 378,
    fk = 379,
    gk = 380,
    hk = 381,
    ik = 382,
    jk = 383,
    kk = 384,
    Wd = 401,
    Xd = 402,
    Yd = 403,
    Zd = 404,
    $d = 405,
    je = 412,
    ke = 413,
    ee = 422,
    fe = 423;
  function lk() {
    return { index: {}, S: {}, set: {}, t: {}, H: {}, a: {}, loop: {} };
  }
  function mk(a) {
    var b;
    b = a.b == Ah ? "_|_" : a.b == Eh ? "'...'" : a.h;
    a.Zb[a.lc++] = " ";
    a.lc == zh && (a.lc = 0);
    for (var c = 0; c < b.length; c++)
      (a.Zb[a.lc++] = b[c]), a.lc == zh && (a.lc = 0);
  }
  function nk(a) {
    var b;
    a.l != $g &&
      ("\n" == a.l && (a.bb++, (a.Vc = 0)),
      (b = a.cf()),
      0 > b && (b = $g),
      a.Vc++,
      b == $g
        ? "\n" == a.l
          ? a.bb--
          : ok(a, "final NL missing before end of file")
        : "\n" != b &&
          (0 <= " \t\n\v\f\r".indexOf(b)
            ? (b = " ")
            : ta(b) &&
              (mk(a), U(a, "control character " + b + " not allowed"))),
      (a.l = b));
  }
  function pk(a) {
    a.h += a.l;
    a.Bb++;
    nk(a);
  }
  function Y(a) {
    function b() {
      mk(a);
      U(a, "keyword s.t. incomplete");
    }
    function c() {
      mk(a);
      U(
        a,
        "cannot convert numeric literal " + a.h + " to floating-point number"
      );
    }
    function d() {
      if ("e" == a.l || "E" == a.l)
        for (
          pk(a),
            ("+" != a.l && "-" != a.l) || pk(a),
            wa(a.l) || (mk(a), U(a, "numeric literal " + a.h + " incomplete"));
          wa(a.l);

        )
          pk(a);
      if (ua(a.l) || "_" == a.l)
        mk(a), U(a, "symbol " + a.h + a.l + "... should be enclosed in quotes");
    }
    a.Hf = a.b;
    a.Gf = a.Bb;
    a.Ff = a.h;
    a.If = a.value;
    if (a.Ue)
      (a.Ue = 0), (a.b = a.Mf), (a.Bb = a.Lf), (a.h = a.Kf), (a.value = a.Nf);
    else {
      for (;;) {
        a.b = 0;
        a.Bb = 0;
        a.h = "";
        for (a.value = 0; " " == a.l || "\n" == a.l; ) nk(a);
        if (a.l == $g) a.b = Ah;
        else if ("#" == a.l) {
          for (; "\n" != a.l && a.l != $g; ) nk(a);
          continue;
        } else if (a.nc || (!ua(a.l) && "_" != a.l))
          if (!a.nc && wa(a.l)) {
            for (a.b = Dh; wa(a.l); ) pk(a);
            var e = !1;
            if ("." == a.l)
              if ((pk(a), "." == a.l))
                a.Bb--,
                  (a.h = a.h.substr(0, a.h.length - 1)),
                  (a.Te = 1),
                  (e = !0);
              else for (; wa(a.l); ) pk(a);
            e || d();
            vg(a.h, function (b) {
              a.value = b;
            }) && c();
          } else if ("'" == a.l || '"' == a.l) {
            var f = a.l,
              g = !1;
            a.b = Eh;
            nk(a);
            e = function () {
              for (;;) {
                if (("\n" == a.l && !g) || a.l == $g)
                  mk(a),
                    U(a, "unexpected end of line; string literal incomplete");
                if (a.l == f)
                  if ((nk(a), a.l == f)) {
                    if (g)
                      if ((nk(a), a.l == f)) {
                        nk(a);
                        break;
                      } else (a.h += '""'), (a.Bb += 2);
                  } else if (g) (a.h += '"'), a.Bb++;
                  else break;
                pk(a);
              }
            };
            a.l == f ? (nk(a), a.l == f && ((g = !0), nk(a), e())) : e();
          } else if (a.nc || "+" != a.l)
            if (a.nc || "-" != a.l)
              if ("*" == a.l)
                (a.b = $h), pk(a), "*" == a.l && ((a.b = bi), pk(a));
              else if ("/" == a.l) {
                if (((a.b = ai), pk(a), "*" == a.l)) {
                  for (nk(a); ; )
                    if (a.l == $g)
                      U(
                        a,
                        "unexpected end of file; comment sequence incomplete"
                      );
                    else if ("*" == a.l) {
                      if ((nk(a), "/" == a.l)) break;
                    } else nk(a);
                  nk(a);
                  continue;
                }
              } else if ("^" == a.l) (a.b = bi), pk(a);
              else if ("<" == a.l)
                (a.b = ci),
                  pk(a),
                  "=" == a.l
                    ? ((a.b = di), pk(a))
                    : ">" == a.l
                    ? ((a.b = hi), pk(a))
                    : "-" == a.l && ((a.b = yi), pk(a));
              else if ("=" == a.l) (a.b = ei), pk(a), "=" == a.l && pk(a);
              else if (">" == a.l)
                (a.b = gi),
                  pk(a),
                  "=" == a.l
                    ? ((a.b = fi), pk(a))
                    : ">" == a.l && ((a.b = wi), pk(a));
              else if ("!" == a.l)
                (a.b = Rh), pk(a), "=" == a.l && ((a.b = hi), pk(a));
              else if ("&" == a.l)
                (a.b = ii), pk(a), "&" == a.l && ((a.b = Fh), pk(a));
              else if ("|" == a.l)
                (a.b = ji), pk(a), "|" == a.l && ((a.b = Sh), pk(a));
              else if (a.nc || "." != a.l)
                if ("," == a.l) (a.b = li), pk(a);
                else if (":" == a.l)
                  (a.b = mi), pk(a), "=" == a.l && ((a.b = oi), pk(a));
                else if (";" == a.l) (a.b = ni), pk(a);
                else if ("(" == a.l) (a.b = qi), pk(a);
                else if (")" == a.l) (a.b = ri), pk(a);
                else if ("[" == a.l) (a.b = si), pk(a);
                else if ("]" == a.l) (a.b = ti), pk(a);
                else if ("{" == a.l) (a.b = ui), pk(a);
                else if ("}" == a.l) (a.b = vi), pk(a);
                else if ("~" == a.l) (a.b = xi), pk(a);
                else if (va(a.l) || 0 <= "+-._".indexOf(a.l)) {
                  for (a.b = Ch; va(a.l) || 0 <= "+-._".indexOf(a.l); ) pk(a);
                  switch (
                    vg(a.h, function (b) {
                      a.value = b;
                    })
                  ) {
                    case 0:
                      a.b = Dh;
                      break;
                    case 1:
                      c();
                  }
                } else mk(a), U(a, "character " + a.l + " not allowed");
              else if (((a.b = ki), pk(a), a.Te))
                (a.b = pi), (a.Bb = 2), (a.h = ".."), (a.Te = 0);
              else if ("." == a.l) (a.b = pi), pk(a);
              else {
                if (wa(a.l)) {
                  a.b = Dh;
                  for (pk(a); wa(a.l); ) pk(a);
                  d();
                  vg(a.h, function (b) {
                    a.value = b;
                  }) && c();
                }
              }
            else (a.b = Zh), pk(a);
          else (a.b = Yh), pk(a);
        else {
          for (a.b = Bh; va(a.l) || "_" == a.l; ) pk(a);
          "and" == a.h
            ? (a.b = Fh)
            : "by" == a.h
            ? (a.b = Gh)
            : "cross" == a.h
            ? (a.b = Hh)
            : "diff" == a.h
            ? (a.b = Ih)
            : "div" == a.h
            ? (a.b = Jh)
            : "else" == a.h
            ? (a.b = Kh)
            : "if" == a.h
            ? (a.b = Lh)
            : "in" == a.h
            ? (a.b = Mh)
            : "Infinity" == a.h
            ? (a.b = Nh)
            : "inter" == a.h
            ? (a.b = Oh)
            : "less" == a.h
            ? (a.b = Ph)
            : "mod" == a.h
            ? (a.b = Qh)
            : "not" == a.h
            ? (a.b = Rh)
            : "or" == a.h
            ? (a.b = Sh)
            : "s" == a.h && "." == a.l
            ? ((a.b = Th),
              pk(a),
              "t" != a.l && b(),
              pk(a),
              "." != a.l && b(),
              pk(a))
            : "symdiff" == a.h
            ? (a.b = Uh)
            : "then" == a.h
            ? (a.b = Vh)
            : "union" == a.h
            ? (a.b = Wh)
            : "within" == a.h && (a.b = Xh);
        }
        break;
      }
      mk(a);
      a.Pf = 0;
    }
  }
  function qk(a) {
    a.Ue = 1;
    a.Mf = a.b;
    a.Lf = a.Bb;
    a.Kf = a.h;
    a.Nf = a.value;
    a.b = a.Hf;
    a.Bb = a.Gf;
    a.h = a.Ff;
    a.value = a.If;
  }
  function rk(a, b) {
    return a.b == Bh && a.h == b;
  }
  function sk(a) {
    return (
      (a.b == Fh && "a" == a.h[0]) ||
      a.b == Gh ||
      a.b == Hh ||
      a.b == Ih ||
      a.b == Jh ||
      a.b == Kh ||
      a.b == Lh ||
      a.b == Mh ||
      a.b == Oh ||
      a.b == Ph ||
      a.b == Qh ||
      (a.b == Rh && "n" == a.h[0]) ||
      (a.b == Sh && "o" == a.h[0]) ||
      a.b == Uh ||
      a.b == Vh ||
      a.b == Wh ||
      a.b == Xh
    );
  }
  function tk(a, b, c, d) {
    var e = {};
    e.Ta = a;
    e.T = 0;
    e.a = lk();
    e.value = {};
    switch (a) {
      case Fi:
        e.a.Q = b.Q;
        break;
      case Gi:
        e.a.M = b.M;
        break;
      case Hi:
        e.a.index.ya = b.index.ya;
        e.a.index.e = b.index.e;
        break;
      case Ii:
      case Ji:
        for (a = b.S.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.S.S = b.S.S;
        e.a.S.list = b.S.list;
        break;
      case Ki:
        for (a = b.set.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.set.set = b.set.set;
        e.a.set.list = b.set.list;
        break;
      case Li:
        for (a = b.t.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.t.t = b.t.t;
        e.a.t.list = b.t.list;
        e.a.t.Ac = b.t.Ac;
        break;
      case Mi:
        for (a = b.H.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.H.H = b.H.H;
        e.a.H.list = b.H.list;
        e.a.H.Ac = b.H.Ac;
        break;
      case Ni:
      case Oi:
        for (a = b.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.list = b.list;
        break;
      case Pi:
        e.a.slice = b.slice;
        break;
      case Qi:
      case Ri:
      case Si:
      case Ti:
        e.T = 1;
        break;
      case Ui:
      case Vi:
      case Wi:
      case Xi:
      case Yi:
      case Zi:
      case $i:
      case aj:
      case bj:
      case cj:
      case dj:
      case ej:
      case fj:
      case gj:
      case hj:
      case ij:
      case jj:
      case kj:
      case lj:
      case mj:
      case nj:
      case oj:
        b.a.x.R = e;
        e.T |= b.a.x.T;
        e.a.a.x = b.a.x;
        break;
      case pj:
      case qj:
      case rj:
      case sj:
      case tj:
      case uj:
      case vj:
      case wj:
      case xj:
      case yj:
      case zj:
      case Aj:
        a == Aj && (e.T = 1);
      case Bj:
        a == Bj && (e.T = 1);
      case Cj:
      case Dj:
      case Ej:
      case Fj:
      case Gj:
      case Hj:
      case Ij:
      case Jj:
      case Kj:
      case Lj:
      case Mj:
      case Nj:
      case Oj:
      case Pj:
      case Qj:
      case Rj:
      case Sj:
      case Tj:
      case Uj:
      case Vj:
      case Xj:
        b.a.x.R = e;
        e.T |= b.a.x.T;
        b.a.y.R = e;
        e.T |= b.a.y.T;
        e.a.a.x = b.a.x;
        e.a.a.y = b.a.y;
        break;
      case Yj:
      case Zj:
      case ak:
        b.a.x.R = e;
        e.T |= b.a.x.T;
        b.a.y.R = e;
        e.T |= b.a.y.T;
        null != b.a.z && ((b.a.z.R = e), (e.T |= b.a.z.T));
        e.a.a.x = b.a.x;
        e.a.a.y = b.a.y;
        e.a.a.z = b.a.z;
        break;
      case bk:
      case ck:
        for (a = b.list; null != a; a = a.e) (a.x.R = e), (e.T |= a.x.T);
        e.a.list = b.list;
        break;
      case dk:
      case ek:
      case fk:
      case gk:
      case hk:
      case ik:
      case jk:
      case kk:
        a = b.loop.domain;
        null != a.code && ((a.code.R = e), (e.T |= a.code.T));
        for (a = a.list; null != a; a = a.e) (a.code.R = e), (e.T |= a.code.T);
        null != b.loop.x && ((b.loop.x.R = e), (e.T |= b.loop.x.T));
        e.a.loop.domain = b.loop.domain;
        e.a.loop.x = b.loop.x;
    }
    e.type = c;
    e.q = d;
    e.R = null;
    e.valid = 0;
    e.value = {};
    return e;
  }
  function Z(a, b, c, d) {
    var e = lk();
    e.a.x = b;
    return tk(a, e, c, d);
  }
  function uk(a, b, c, d, e) {
    var f = lk();
    f.a.x = b;
    f.a.y = c;
    return tk(a, f, d, e);
  }
  function vk(a, b, c, d, e, f) {
    var g = lk();
    g.a.x = b;
    g.a.y = c;
    g.a.z = d;
    return tk(a, g, e, f);
  }
  function wk(a, b) {
    var c = {},
      d;
    c.x = b;
    c.e = null;
    if (null == a) a = c;
    else {
      for (d = a; null != d.e; d = d.e);
      d.e = c;
    }
    return a;
  }
  function xk(a) {
    var b;
    for (b = 0; null != a; a = a.e) b++;
    return b;
  }
  function yk(a) {
    var b,
      c,
      d,
      e,
      f,
      g,
      h,
      k,
      l,
      p = lk(),
      m = a.V[a.h];
    null == m && U(a, a.h + " not defined");
    switch (m.type) {
      case kh:
        b = m.link;
        k = b.name;
        l = 0;
        break;
      case uh:
        c = m.link;
        k = c.name;
        l = c.q;
        0 == c.X && (c.X = 1);
        break;
      case sh:
        d = m.link;
        k = d.name;
        l = d.q;
        break;
      case yh:
        e = m.link;
        k = e.name;
        l = e.q;
        break;
      case ch:
        (f = m.link), (k = f.name), (l = f.q);
    }
    Y(a);
    if (a.b == si) {
      0 == l && U(a, k + " cannot be subscripted");
      Y(a);
      for (var q = null; ; )
        if (
          ((g = zk(a)),
          g.type == N && (g = Z(Vi, g, T, 0)),
          g.type != T && U(a, "subscript expression has invalid type"),
          (q = wk(q, g)),
          a.b == li)
        )
          Y(a);
        else if (a.b == ti) break;
        else U(a, "syntax error in subscript list");
      g = q;
      l != xk(g) &&
        U(
          a,
          k +
            " must have " +
            l +
            " subscript" +
            (1 == l ? "" : "s") +
            " rather than " +
            xk(g)
        );
      Y(a);
    } else 0 != l && U(a, k + " must be subscripted"), (g = null);
    l = a.Ob || m.type != yh ? Di : zi;
    a.b == ki &&
      (Y(a),
      a.b != Bh && U(a, "invalid use of period"),
      m.type != yh && m.type != ch && U(a, k + " cannot have a suffix"),
      "lb" == a.h
        ? (l = Ai)
        : "ub" == a.h
        ? (l = Bi)
        : "status" == a.h
        ? (l = Ci)
        : "val" == a.h
        ? (l = Di)
        : "dual" == a.h
        ? (l = Ei)
        : U(a, "suffix ." + a.h + " invalid"),
      Y(a));
    switch (m.type) {
      case kh:
        p.index.ya = b;
        p.index.e = b.list;
        h = tk(Hi, p, T, 0);
        b.list = h;
        break;
      case uh:
        p.set.set = c;
        p.set.list = g;
        h = tk(Ki, p, fh, c.X);
        break;
      case sh:
        p.S.S = d;
        p.S.list = g;
        h = d.type == T ? tk(Ji, p, T, 0) : tk(Ii, p, N, 0);
        break;
      case yh:
        a.Ob ||
          (l != Ci && l != Di && l != Ei) ||
          U(
            a,
            "invalid reference to status, primal value, or dual value of variable " +
              e.name +
              " above solve statement"
          );
        p.t.t = e;
        p.t.list = g;
        p.t.Ac = l;
        h = tk(Li, p, l == zi ? jh : N, 0);
        break;
      case ch:
        a.Ob ||
          (l != Ci && l != Di && l != Ei) ||
          U(
            a,
            "invalid reference to status, primal value, or dual value of " +
              (f.type == ch ? "constraint" : "objective") +
              " " +
              f.name +
              " above solve statement"
          ),
          (p.H.H = f),
          (p.H.list = g),
          (p.H.Ac = l),
          (h = tk(Mi, p, N, 0));
    }
    return h;
  }
  function Ak(a, b) {
    var c = zk(a);
    c.type == T && (c = Z(Ui, c, N, 0));
    c.type != N && U(a, "argument for " + b + " has invalid type");
    return c;
  }
  function Bk(a, b) {
    var c = zk(a);
    c.type == N && (c = Z(Vi, c, T, 0));
    c.type != T && U(a, "argument for " + b + " has invalid type");
    return c;
  }
  function Ck(a, b, c) {
    var d = {};
    d.name = b;
    d.code = c;
    d.value = null;
    d.list = null;
    d.e = null;
    if (null == a.list) a.list = d;
    else {
      for (a = a.list; null != a.e; a = a.e);
      a.e = d;
    }
  }
  function Dk(a) {
    var b,
      c = lk(),
      d = Array(21);
    ia(d, 0, 21);
    var e,
      f,
      g,
      h = 0;
    e = a.Pf;
    Y(a);
    for (g = 1; ; g++) {
      20 < g && U(a, "too many components within parentheses");
      var k = function () {
        b = Ek(a);
        if (a.b == li || 1 < g)
          b.type == N && (b = Z(Vi, b, T, 0)),
            b.type != T && U(a, "component expression has invalid type");
        d[g].name = null;
        d[g].code = b;
      };
      if (a.b == Bh)
        if (
          (Y(a),
          (f = a.b),
          qk(a),
          !e || (f != li && f != ri) || null != a.V[a.h])
        )
          k();
        else {
          for (f = 1; f < g; f++)
            null != d[f].name &&
              d[f].name == a.h &&
              U(a, "duplicate dummy index " + a.h + " not allowed");
          d[g].name = a.h;
          d[g].code = null;
          Y(a);
          h = 1;
          1 == g && a.b == ri && U(a, d[g].name + " not defined");
        }
      else k();
      if (a.b == li) Y(a);
      else if (a.b == ri) break;
      else U(a, "right parenthesis missing where expected");
    }
    if (1 != g || h)
      if (h) {
        c.slice = {};
        for (f = 1; f <= g; f++) Ck(c.slice, d[f].name, d[f].code);
        b = tk(Pi, c, xh, g);
      } else {
        c.list = null;
        for (f = 1; f <= g; f++) c.list = wk(c.list, d[f].code);
        b = tk(Ni, c, xh, g);
      }
    else b = d[1].code;
    Y(a);
    h && a.b != Mh && U(a, "keyword in missing where expected");
    e &&
      a.b == Mh &&
      !h &&
      (1 == g
        ? U(a, "syntax error in indexing expression")
        : U(a, "0-ary slice not allowed"));
    return b;
  }
  function Fk(a) {
    var b, c, d, e;
    Y(a);
    a.b == vi && U(a, "empty indexing expression not allowed");
    for (b = {}; ; ) {
      e = c = null;
      a.b == Bh
        ? (Y(a),
          (d = a.b),
          qk(a),
          d == Mh &&
            null == a.V[a.h] &&
            ((c = {}), (d = a.h), Ck(c, d, null), Y(a), Y(a)))
        : a.b == qi &&
          ((a.Pf = 1),
          (e = Gk(a)),
          e.Ta == Pi && ((c = e.a.slice), (e = null), Y(a)));
      null == e && (e = Gk(a));
      if (e.type != fh) {
        null != c && U(a, "domain expression has invalid type");
        d = a;
        var f = lk(),
          g = void 0;
        f.list = null;
        for (g = 1; ; g++) {
          e.type == N && (e = Z(Vi, e, T, 0));
          e.type == T && (e = Z(Xi, e, xh, 1));
          e.type != xh && U(d, "member expression has invalid type");
          null != f.list &&
            f.list.x.q != e.q &&
            U(
              d,
              "member " +
                (g - 1) +
                " has " +
                f.list.x.q +
                " component" +
                (1 == f.list.x.q ? "" : "s") +
                " while member " +
                g +
                " has " +
                e.q +
                " component" +
                (1 == e.q ? "" : "s")
            );
          f.list = wk(f.list, e);
          if (d.b == li) Y(d);
          else if (d.b == vi) break;
          else U(d, "syntax error in literal set");
          e = zk(d);
        }
        e = tk(Oi, f, fh, f.list.x.q);
      }
      if (null == c) for (c = {}, d = 1; d <= e.q; d++) Ck(c, null, null);
      f = 0;
      for (d = c.list; null != d; d = d.e) f++;
      f != e.q &&
        U(
          a,
          f +
            " " +
            (1 == f ? "index" : "indices") +
            " specified for set of dimension " +
            e.q
        );
      c.code = e;
      e = b;
      d = c;
      f = void 0;
      if (null == e.list) e.list = d;
      else {
        for (f = e.list; null != f.e; f = f.e);
        f.e = d;
      }
      for (d = c.list; null != d; d = d.e)
        null != d.name && (a.V[d.name] = { type: kh, link: d });
      if (a.b == li) Y(a);
      else if (a.b == mi || a.b == vi) break;
      else U(a, "syntax error in indexing expression");
    }
    a.b == mi &&
      (Y(a),
      (e = Ek(a)),
      e.type == T && (e = Z(Ui, e, N, 0)),
      e.type == N && (e = Z(Wi, e, nh, 0)),
      e.type != nh && U(a, "expression following colon has invalid type"),
      (b.code = e),
      a.b != vi && U(a, "syntax error in indexing expression"));
    Y(a);
    return b;
  }
  function Hk(a, b) {
    var c, d;
    for (c = b.list; null != c; c = c.e)
      for (d = c.list; null != d; d = d.e) null != d.name && delete a.V[d.name];
  }
  function Ik(a) {
    var b, c;
    for (b = a.a.loop.domain.list; null != b; b = b.e)
      for (c = b.list; null != c; c = c.e) null != c.code && (c.code.R = a);
  }
  function Jk(a) {
    function b() {
      U(a, "integrand following " + f + "{...} has invalid type");
    }
    var c,
      d = lk(),
      e,
      f;
    "sum" == a.h
      ? (e = dk)
      : "prod" == a.h
      ? (e = ek)
      : "min" == a.h
      ? (e = fk)
      : "max" == a.h
      ? (e = gk)
      : "forall" == a.h
      ? (e = hk)
      : "exists" == a.h
      ? (e = ik)
      : "setof" == a.h
      ? (e = jk)
      : U(a, "operator " + a.h + " unknown");
    f = a.h;
    Y(a);
    d.loop.domain = Fk(a);
    switch (e) {
      case dk:
      case ek:
      case fk:
      case gk:
        d.loop.x = Kk(a);
        d.loop.x.type == T && (d.loop.x = Z(Ui, d.loop.x, N, 0));
        d.loop.x.type == N || (e == dk && d.loop.x.type == jh) || b();
        c = tk(e, d, d.loop.x.type, 0);
        break;
      case hk:
      case ik:
        d.loop.x = Lk(a);
        d.loop.x.type == T && (d.loop.x = Z(Ui, d.loop.x, N, 0));
        d.loop.x.type == N && (d.loop.x = Z(Wi, d.loop.x, nh, 0));
        d.loop.x.type != nh && b();
        c = tk(e, d, nh, 0);
        break;
      case jk:
        (d.loop.x = zk(a)),
          d.loop.x.type == N && (d.loop.x = Z(Vi, d.loop.x, T, 0)),
          d.loop.x.type == T && (d.loop.x = Z(Xi, d.loop.x, xh, 1)),
          d.loop.x.type != xh && b(),
          (c = tk(e, d, fh, d.loop.x.q));
    }
    Hk(a, d.loop.domain);
    Ik(c);
    return c;
  }
  function Mk(a) {
    var b = 0;
    for (a = a.list; null != a; a = a.e)
      for (var c = a.list; null != c; c = c.e) null == c.code && b++;
    return b;
  }
  function Nk(a, b) {
    U(a, "operand preceding " + b + " has invalid type");
  }
  function Ok(a, b) {
    U(a, "operand following " + b + " has invalid type");
  }
  function Pk(a, b, c, d) {
    U(
      a,
      "operands preceding and following " +
        b +
        " have different dimensions " +
        c +
        " and " +
        d +
        ", respectively"
    );
  }
  function Qk(a) {
    var b, c;
    if (a.b == Dh)
      (b = lk()), (b.Q = a.value), (b = tk(Fi, b, N, 0)), Y(a), (c = b);
    else if (a.b == Nh) (b = lk()), (b.Q = s), (c = tk(Fi, b, N, 0)), Y(a);
    else if (a.b == Eh)
      (b = lk()), (b.M = a.h), (b = tk(Gi, b, T, 0)), Y(a), (c = b);
    else if (a.b == Bh)
      switch ((Y(a), (c = a.b), qk(a), c)) {
        case si:
          c = yk(a);
          break;
        case qi:
          c = lk();
          var d;
          "abs" == a.h
            ? (b = bj)
            : "ceil" == a.h
            ? (b = cj)
            : "floor" == a.h
            ? (b = dj)
            : "exp" == a.h
            ? (b = ej)
            : "log" == a.h
            ? (b = fj)
            : "log10" == a.h
            ? (b = gj)
            : "sqrt" == a.h
            ? (b = hj)
            : "sin" == a.h
            ? (b = ij)
            : "cos" == a.h
            ? (b = jj)
            : "atan" == a.h
            ? (b = kj)
            : "min" == a.h
            ? (b = bk)
            : "max" == a.h
            ? (b = ck)
            : "round" == a.h
            ? (b = lj)
            : "trunc" == a.h
            ? (b = mj)
            : "Irand224" == a.h
            ? (b = Qi)
            : "Uniform01" == a.h
            ? (b = Ri)
            : "Uniform" == a.h
            ? (b = Aj)
            : "Normal01" == a.h
            ? (b = Si)
            : "Normal" == a.h
            ? (b = Bj)
            : "card" == a.h
            ? (b = nj)
            : "length" == a.h
            ? (b = oj)
            : "substr" == a.h
            ? (b = Uj)
            : "str2time" == a.h
            ? (b = Vj)
            : "time2str" == a.h
            ? (b = Xj)
            : "gmtime" == a.h
            ? (b = Ti)
            : U(a, "function " + a.h + " unknown");
          d = a.h;
          Y(a);
          Y(a);
          if (b == bk || b == ck)
            for (c.list = null; ; )
              if (((c.list = wk(c.list, Ak(a, d))), a.b == li)) Y(a);
              else if (a.b == ri) break;
              else U(a, "syntax error in argument list for " + d);
          else if (b == Qi || b == Ri || b == Si || b == Ti)
            a.b != ri && U(a, d + " needs no arguments");
          else if (b == Aj || b == Bj)
            (c.a.x = Ak(a, d)),
              a.b != li &&
                (a.b == ri
                  ? U(a, d + " needs two arguments")
                  : U(a, "syntax error in argument for " + d)),
              Y(a),
              (c.a.y = Ak(a, d)),
              a.b == li
                ? U(a, d + " needs two argument")
                : a.b != ri && U(a, "syntax error in argument for " + d);
          else if (b == kj || b == lj || b == mj) {
            c.a.x = Ak(a, d);
            if (a.b == li) {
              switch (b) {
                case kj:
                  b = xj;
                  break;
                case lj:
                  b = yj;
                  break;
                case mj:
                  b = zj;
              }
              Y(a);
              c.a.y = Ak(a, d);
            }
            a.b == li
              ? U(a, d + " needs one or two arguments")
              : a.b != ri && U(a, "syntax error in argument for " + d);
          } else if (b == Uj)
            (c.a.x = Bk(a, d)),
              a.b != li &&
                (a.b == ri
                  ? U(a, d + " needs two or three arguments")
                  : U(a, "syntax error in argument for " + d)),
              Y(a),
              (c.a.y = Ak(a, d)),
              a.b == li && ((b = ak), Y(a), (c.a.z = Ak(a, d))),
              a.b == li
                ? U(a, d + " needs two or three arguments")
                : a.b != ri && U(a, "syntax error in argument for " + d);
          else if (b == Vj)
            (c.a.x = Bk(a, d)),
              a.b != li &&
                (a.b == ri
                  ? U(a, d + " needs two arguments")
                  : U(a, "syntax error in argument for " + d)),
              Y(a),
              (c.a.y = Bk(a, d)),
              a.b == li
                ? U(a, d + " needs two argument")
                : a.b != ri && U(a, "syntax error in argument for " + d);
          else if (b == Xj)
            (c.a.x = Ak(a, d)),
              a.b != li &&
                (a.b == ri
                  ? U(a, d + " needs two arguments")
                  : U(a, "syntax error in argument for " + d)),
              Y(a),
              (c.a.y = Bk(a, d)),
              a.b == li
                ? U(a, d + " needs two argument")
                : a.b != ri && U(a, "syntax error in argument for " + d);
          else {
            var e = c.a,
              f;
            b == nj
              ? ((f = Gk(a)),
                f.type != fh && U(a, "argument for " + d + " has invalid type"))
              : (f = b == oj ? Bk(a, d) : Ak(a, d));
            e.x = f;
            a.b == li
              ? U(a, d + " needs one argument")
              : a.b != ri && U(a, "syntax error in argument for " + d);
          }
          b = b == Uj || b == ak || b == Xj ? tk(b, c, T, 0) : tk(b, c, N, 0);
          Y(a);
          c = b;
          break;
        case ui:
          c = Jk(a);
          break;
        default:
          c = yk(a);
      }
    else if (a.b == qi) c = Dk(a);
    else if (a.b == ui)
      (b = lk()),
        Y(a),
        a.b == vi
          ? ((b.list = null), (b = tk(Oi, b, fh, 1)), Y(a))
          : (qk(a),
            (b.loop.domain = Fk(a)),
            (b.loop.x = null),
            Hk(a, b.loop.domain),
            (b = tk(kk, b, fh, Mk(b.loop.domain))),
            Ik(b)),
        (c = b);
    else if (a.b == Lh) {
      Y(a);
      b = Ek(a);
      b.type == T && (b = Z(Ui, b, N, 0));
      b.type == N && (b = Z(Wi, b, nh, 0));
      b.type != nh && U(a, "expression following if has invalid type");
      a.b != Vh && U(a, "keyword then missing where expected");
      Y(a);
      c = Gk(a);
      c.type != N &&
        c.type != T &&
        c.type != fh &&
        c.type != jh &&
        U(a, "expression following then has invalid type");
      if (a.b != Kh)
        c.type == fh && U(a, "keyword else missing where expected"), (d = null);
      else {
        Y(a);
        d = Gk(a);
        d.type != N &&
          d.type != T &&
          d.type != fh &&
          d.type != jh &&
          U(a, "expression following else has invalid type");
        if (c.type == jh || d.type == jh)
          c.type == T && (c = Z(Ui, c, N, 0)),
            c.type == N && (c = Z(Yi, c, jh, 0)),
            d.type == T && (d = Z(Ui, d, N, 0)),
            d.type == N && (d = Z(Yi, d, jh, 0));
        if (c.type == T || d.type == T)
          c.type == N && (c = Z(Vi, c, T, 0)),
            d.type == N && (d = Z(Vi, d, T, 0));
        c.type != d.type &&
          U(a, "expressions following then and else have incompatible types");
        c.q != d.q &&
          U(
            a,
            "expressions following then and else have different dimensions " +
              c.q +
              " and " +
              d.q +
              ", respectively"
          );
      }
      c = vk(Zj, b, c, d, c.type, c.q);
    } else
      sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "syntax error in expression");
    a.b == bi &&
      ((d = a.h),
      c.type == T && (c = Z(Ui, c, N, 0)),
      c.type != N && Nk(a, d),
      Y(a),
      (b = a.b == Yh || a.b == Zh ? Rk(a) : Qk(a)),
      b.type == T && (b = Z(Ui, b, N, 0)),
      b.type != N && Ok(a, d),
      (c = uk(wj, c, b, N, 0)));
    return c;
  }
  function Rk(a) {
    var b;
    a.b == Yh
      ? (Y(a),
        (b = Qk(a)),
        b.type == T && (b = Z(Ui, b, N, 0)),
        b.type != N && b.type != jh && Ok(a, "+"),
        (b = Z(Zi, b, b.type, 0)))
      : a.b == Zh
      ? (Y(a),
        (b = Qk(a)),
        b.type == T && (b = Z(Ui, b, N, 0)),
        b.type != N && b.type != jh && Ok(a, "-"),
        (b = Z($i, b, b.type, 0)))
      : (b = Qk(a));
    return b;
  }
  function Kk(a) {
    for (var b, c = Rk(a); ; )
      if (a.b == $h)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && c.type != jh && Nk(a, "*"),
          Y(a),
          (b = Rk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && b.type != jh && Ok(a, "*"),
          c.type == jh &&
            b.type == jh &&
            U(a, "multiplication of linear forms not allowed"),
          (c =
            c.type == N && b.type == N
              ? uk(sj, c, b, N, 0)
              : uk(sj, c, b, jh, 0));
      else if (a.b == ai)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && c.type != jh && Nk(a, "/"),
          Y(a),
          (b = Rk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && Ok(a, "/"),
          (c = c.type == N ? uk(tj, c, b, N, 0) : uk(tj, c, b, jh, 0));
      else if (a.b == Jh)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && Nk(a, "div"),
          Y(a),
          (b = Rk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && Ok(a, "div"),
          (c = uk(uj, c, b, N, 0));
      else if (a.b == Qh)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && Nk(a, "mod"),
          Y(a),
          (b = Rk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && Ok(a, "mod"),
          (c = uk(vj, c, b, N, 0));
      else break;
    return c;
  }
  function Sk(a) {
    for (var b, c = Kk(a); ; )
      if (a.b == Yh)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && c.type != jh && Nk(a, "+"),
          Y(a),
          (b = Kk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && b.type != jh && Ok(a, "+"),
          c.type == N && b.type == jh && (c = Z(Yi, c, jh, 0)),
          c.type == jh && b.type == N && (b = Z(Yi, b, jh, 0)),
          (c = uk(pj, c, b, c.type, 0));
      else if (a.b == Zh)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && c.type != jh && Nk(a, "-"),
          Y(a),
          (b = Kk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && b.type != jh && Ok(a, "-"),
          c.type == N && b.type == jh && (c = Z(Yi, c, jh, 0)),
          c.type == jh && b.type == N && (b = Z(Yi, b, jh, 0)),
          (c = uk(qj, c, b, c.type, 0));
      else if (a.b == Ph)
        c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && Nk(a, "less"),
          Y(a),
          (b = Kk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type != N && Ok(a, "less"),
          (c = uk(rj, c, b, N, 0));
      else break;
    return c;
  }
  function zk(a) {
    for (var b, c = Sk(a); ; )
      if (a.b == ii)
        c.type == N && (c = Z(Vi, c, T, 0)),
          c.type != T && Nk(a, "&"),
          Y(a),
          (b = Sk(a)),
          b.type == N && (b = Z(Vi, b, T, 0)),
          b.type != T && Ok(a, "&"),
          (c = uk(Cj, c, b, T, 0));
      else break;
    return c;
  }
  function Tk(a) {
    var b,
      c,
      d = zk(a);
    a.b == pi &&
      (d.type == T && (d = Z(Ui, d, N, 0)),
      d.type != N && Nk(a, ".."),
      Y(a),
      (b = zk(a)),
      b.type == T && (b = Z(Ui, b, N, 0)),
      b.type != N && Ok(a, ".."),
      a.b == Gh
        ? (Y(a),
          (c = zk(a)),
          c.type == T && (c = Z(Ui, c, N, 0)),
          c.type != N && Ok(a, "by"))
        : (c = null),
      (d = vk(Yj, d, b, c, fh, 1)));
    return d;
  }
  function Uk(a) {
    for (var b, c = Tk(a); ; )
      if (a.b == Hh)
        c.type != fh && Nk(a, "cross"),
          Y(a),
          (b = Tk(a)),
          b.type != fh && Ok(a, "cross"),
          (c = uk(Pj, c, b, fh, c.q + b.q));
      else break;
    return c;
  }
  function Vk(a) {
    for (var b, c = Uk(a); ; )
      if (a.b == Oh)
        c.type != fh && Nk(a, "inter"),
          Y(a),
          (b = Uk(a)),
          b.type != fh && Ok(a, "inter"),
          c.q != b.q && Pk(a, "inter", c.q, b.q),
          (c = uk(Oj, c, b, fh, c.q));
      else break;
    return c;
  }
  function Gk(a) {
    for (var b, c = Vk(a); ; )
      if (a.b == Wh)
        c.type != fh && Nk(a, "union"),
          Y(a),
          (b = Vk(a)),
          b.type != fh && Ok(a, "union"),
          c.q != b.q && Pk(a, "union", c.q, b.q),
          (c = uk(Lj, c, b, fh, c.q));
      else if (a.b == Ih)
        c.type != fh && Nk(a, "diff"),
          Y(a),
          (b = Vk(a)),
          b.type != fh && Ok(a, "diff"),
          c.q != b.q && Pk(a, "diff", c.q, b.q),
          (c = uk(Mj, c, b, fh, c.q));
      else if (a.b == Uh)
        c.type != fh && Nk(a, "symdiff"),
          Y(a),
          (b = Vk(a)),
          b.type != fh && Ok(a, "symdiff"),
          c.q != b.q && Pk(a, "symdiff", c.q, b.q),
          (c = uk(Nj, c, b, fh, c.q));
      else break;
    return c;
  }
  function Wk(a) {
    var b,
      c = -1,
      d = "",
      e = Gk(a);
    switch (a.b) {
      case ci:
        c = Dj;
        break;
      case di:
        c = Ej;
        break;
      case ei:
        c = Fj;
        break;
      case fi:
        c = Gj;
        break;
      case gi:
        c = Hj;
        break;
      case hi:
        c = Ij;
        break;
      case Mh:
        c = Qj;
        break;
      case Xh:
        c = Sj;
        break;
      case Rh:
        d = a.h;
        Y(a);
        a.b == Mh
          ? (c = Rj)
          : a.b == Xh
          ? (c = Tj)
          : U(a, "invalid use of " + d);
        d += " ";
        break;
      default:
        return e;
    }
    d += a.h;
    switch (c) {
      case Fj:
      case Ij:
      case Dj:
      case Ej:
      case Hj:
      case Gj:
        e.type != N && e.type != T && Nk(a, d);
        Y(a);
        b = Gk(a);
        b.type != N && b.type != T && Ok(a, d);
        e.type == N && b.type == T && (e = Z(Vi, e, T, 0));
        e.type == T && b.type == N && (b = Z(Vi, b, T, 0));
        e = uk(c, e, b, nh, 0);
        break;
      case Qj:
      case Rj:
        e.type == N && (e = Z(Vi, e, T, 0));
        e.type == T && (e = Z(Xi, e, xh, 1));
        e.type != xh && Nk(a, d);
        Y(a);
        b = Gk(a);
        b.type != fh && Ok(a, d);
        e.q != b.q && Pk(a, d, e.q, b.q);
        e = uk(c, e, b, nh, 0);
        break;
      case Sj:
      case Tj:
        e.type != fh && Nk(a, d),
          Y(a),
          (b = Gk(a)),
          b.type != fh && Ok(a, d),
          e.q != b.q && Pk(a, d, e.q, b.q),
          (e = uk(c, e, b, nh, 0));
    }
    return e;
  }
  function Xk(a) {
    var b, c;
    a.b == Rh
      ? ((c = a.h),
        Y(a),
        (b = Wk(a)),
        b.type == T && (b = Z(Ui, b, N, 0)),
        b.type == N && (b = Z(Wi, b, nh, 0)),
        b.type != nh && Ok(a, c),
        (b = Z(aj, b, nh, 0)))
      : (b = Wk(a));
    return b;
  }
  function Lk(a) {
    for (var b, c = "", d = Xk(a); ; )
      if (a.b == Fh)
        (c = a.h),
          d.type == T && (d = Z(Ui, d, N, 0)),
          d.type == N && (d = Z(Wi, d, nh, 0)),
          d.type != nh && Nk(a, c),
          Y(a),
          (b = Xk(a)),
          b.type == T && (b = Z(Ui, b, N, 0)),
          b.type == N && (b = Z(Wi, b, nh, 0)),
          b.type != nh && Ok(a, c),
          (d = uk(Jj, d, b, nh, 0));
      else break;
    return d;
  }
  function Ek(a) {
    for (var b, c = Lk(a); ; )
      if (a.b == Sh) {
        var d = a.h;
        c.type == T && (c = Z(Ui, c, N, 0));
        c.type == N && (c = Z(Wi, c, nh, 0));
        c.type != nh && Nk(a, d);
        Y(a);
        b = Lk(a);
        b.type == T && (b = Z(Ui, b, N, 0));
        b.type == N && (b = Z(Wi, b, nh, 0));
        b.type != nh && Ok(a, d);
        c = uk(Kj, c, b, nh, 0);
      } else break;
    return c;
  }
  function Yk(a) {
    function b() {
      U(a, "at most one := or default/data allowed");
    }
    function c() {
      U(a, a.h + " not a plain set");
    }
    function d() {
      U(a, "dimension of " + a.h + " too small");
    }
    function e() {
      U(a, "component number must be integer between 1 and " + k.set.X);
    }
    var f,
      g,
      h = 0,
      k;
    Y(a);
    a.b != Bh &&
      (sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "symbolic name missing where expected"));
    null != a.V[a.h] && U(a, a.h + " multiply declared");
    f = {};
    f.name = a.h;
    f.Ib = null;
    f.q = 0;
    f.domain = null;
    f.X = 0;
    f.Bf = null;
    f.assign = null;
    f.xa = null;
    f.Yc = null;
    f.data = 0;
    f.O = null;
    Y(a);
    a.b == Eh && ((f.Ib = a.h), Y(a));
    a.b == ui && ((f.domain = Fk(a)), (f.q = Mk(f.domain)));
    g = a.V[f.name] = {};
    g.type = uh;
    for (g.link = f; ; ) {
      if (a.b == li) Y(a);
      else if (a.b == ni) break;
      if (rk(a, "dimen")) {
        var l;
        Y(a);
        (a.b == Dh &&
          1 <= a.value &&
          20 >= a.value &&
          Math.floor(a.value) == a.value) ||
          U(a, "dimension must be integer between 1 and 20");
        l = (a.value + 0.5) | 0;
        h && U(a, "at most one dimension attribute allowed");
        0 < f.X &&
          U(
            a,
            "dimension " +
              l +
              " conflicts with dimension " +
              f.X +
              " already determined"
          );
        f.X = l;
        h = 1;
        Y(a);
      } else if (a.b == Xh || a.b == Mh) {
        a.b != Mh ||
          a.rg ||
          (ok(a, "keyword in understood as within"), (a.rg = 1));
        Y(a);
        l = { code: null, e: null };
        if (null == f.Bf) f.Bf = l;
        else {
          for (g = f.Bf; null != g.e; g = g.e);
          g.e = l;
        }
        l.code = Gk(a);
        l.code.type != fh &&
          U(a, "expression following within has invalid type");
        0 == f.X && (f.X = l.code.q);
        f.X != l.code.q &&
          U(
            a,
            "set expression following within must have dimension " +
              f.X +
              " rather than " +
              l.code.q
          );
      } else if (a.b == oi)
        (null == f.assign && null == f.xa && null == f.Yc) || b(),
          Y(a),
          (f.assign = Gk(a)),
          f.assign.type != fh &&
            U(a, "expression following := has invalid type"),
          0 == f.X && (f.X = f.assign.q),
          f.X != f.assign.q &&
            U(
              a,
              "set expression following := must have dimension " +
                f.X +
                " rather than " +
                f.assign.q
            );
      else if (rk(a, "default"))
        (null == f.assign && null == f.xa) || b(),
          Y(a),
          (f.xa = Gk(a)),
          f.xa.type != fh &&
            U(a, "expression following default has invalid type"),
          0 == f.X && (f.X = f.xa.q),
          f.X != f.xa.q &&
            U(
              a,
              "set expression following default must have dimension " +
                f.X +
                " rather than " +
                f.xa.q
            );
      else if (rk(a, "data")) {
        var p = 0;
        l = Array(20);
        (null == f.assign && null == f.Yc) || b();
        Y(a);
        f.Yc = k = {};
        a.b != Bh &&
          (sk(a)
            ? U(a, "invalid use of reserved keyword " + a.h)
            : U(a, "set name missing where expected"));
        g = a.V[a.h];
        null == g && U(a, a.h + " not defined");
        g.type != uh && c();
        k.set = g.link;
        0 != k.set.q && c();
        k.set == f && U(a, "set cannot be initialized by itself");
        f.q >= k.set.X && d();
        0 == f.X && (f.X = k.set.X - f.q);
        f.q + f.X > k.set.X
          ? d()
          : f.q + f.X < k.set.X && U(a, "dimension of " + a.h + " too big");
        Y(a);
        a.b == qi ? Y(a) : U(a, "left parenthesis missing where expected");
        for (g = 0; g < k.set.X; g++) l[g] = 0;
        for (g = 0; ; )
          if (
            (a.b != Dh && U(a, "component number missing where expected"),
            0 !=
              wg(a.h, function (a) {
                p = a;
              }) && e(),
            (1 <= p && p <= k.set.X) || e(),
            0 != l[p - 1] && U(a, "component " + p + " multiply specified"),
            (k.Z[g++] = p),
            (l[p - 1] = 1),
            Y(a),
            a.b == li)
          )
            Y(a);
          else if (a.b == ri) break;
          else U(a, "syntax error in data attribute");
        g < k.set.X &&
          U(a, "there are must be " + k.set.X + " components rather than " + g);
        Y(a);
      } else U(a, "syntax error in set statement");
    }
    null != f.domain && Hk(a, f.domain);
    0 == f.X && (f.X = 1);
    Y(a);
    return f;
  }
  function Zk(a) {
    function b() {
      g && U(a, "at most one binary allowed");
      d.type == T && U(a, "symbolic parameter cannot be binary");
      d.type = ah;
      g = 1;
      Y(a);
    }
    function c() {
      U(a, "at most one := or default allowed");
    }
    var d,
      e,
      f = 0,
      g = 0,
      h = 0;
    Y(a);
    a.b != Bh &&
      (sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "symbolic name missing where expected"));
    null != a.V[a.h] && U(a, a.h + " multiply declared");
    d = {};
    d.name = a.h;
    d.Ib = null;
    d.q = 0;
    d.domain = null;
    d.type = N;
    d.wd = null;
    d.qa = null;
    d.assign = null;
    d.xa = null;
    d.data = 0;
    d.Wc = null;
    d.O = null;
    Y(a);
    a.b == Eh && ((d.Ib = a.h), Y(a));
    a.b == ui && ((d.domain = Fk(a)), (d.q = Mk(d.domain)));
    e = a.V[d.name] = {};
    e.type = sh;
    for (e.link = d; ; ) {
      if (a.b == li) Y(a);
      else if (a.b == ni) break;
      if (rk(a, "integer"))
        f && U(a, "at most one integer allowed"),
          d.type == T && U(a, "symbolic parameter cannot be integer"),
          d.type != ah && (d.type = mh),
          (f = 1),
          Y(a);
      else if (rk(a, "binary")) b();
      else if (rk(a, "logical"))
        a.Ke || (ok(a, "keyword logical understood as binary"), (a.Ke = 1)),
          b();
      else if (rk(a, "symbolic"))
        h && U(a, "at most one symbolic allowed"),
          d.type != N && U(a, "integer or binary parameter cannot be symbolic"),
          (null == d.wd && null == d.qa && null == d.assign && null == d.xa) ||
            U(
              a,
              "keyword symbolic must precede any other parameter attributes"
            ),
          (d.type = T),
          (h = 1),
          Y(a);
      else if (
        a.b == ci ||
        a.b == di ||
        a.b == ei ||
        a.b == fi ||
        a.b == gi ||
        a.b == hi
      ) {
        var k,
          l = {};
        switch (a.b) {
          case ci:
            l.jd = Dj;
            k = a.h;
            break;
          case di:
            l.jd = Ej;
            k = a.h;
            break;
          case ei:
            l.jd = Fj;
            k = a.h;
            break;
          case fi:
            l.jd = Gj;
            k = a.h;
            break;
          case gi:
            l.jd = Hj;
            k = a.h;
            break;
          case hi:
            (l.jd = Ij), (k = a.h);
        }
        l.code = null;
        l.e = null;
        if (null == d.wd) d.wd = l;
        else {
          for (e = d.wd; null != e.e; e = e.e);
          e.e = l;
        }
        Y(a);
        l.code = zk(a);
        l.code.type != N &&
          l.code.type != T &&
          U(a, "expression following " + k + " has invalid type");
        d.type != T && l.code.type == T && (l.code = Z(Ui, l.code, N, 0));
        d.type == T && l.code.type != T && (l.code = Z(Vi, l.code, T, 0));
      } else if (a.b == Mh || a.b == Xh) {
        a.b != Xh ||
          a.qg ||
          (ok(a, "keyword within understood as in"), (a.qg = 1));
        Y(a);
        l = { code: null, e: null };
        if (null == d.qa) d.qa = l;
        else {
          for (e = d.qa; null != e.e; e = e.e);
          e.e = l;
        }
        l.code = Gk(a);
        l.code.type != fh && U(a, "expression following in has invalid type");
        1 != l.code.q &&
          U(
            a,
            "set expression following in must have dimension 1 rather than " +
              l.code.q
          );
      } else
        a.b == oi
          ? ((null == d.assign && null == d.xa) || c(),
            Y(a),
            (d.assign = zk(a)),
            d.assign.type != N &&
              d.assign.type != T &&
              U(a, "expression following := has invalid type"),
            d.type != T &&
              d.assign.type == T &&
              (d.assign = Z(Ui, d.assign, N, 0)),
            d.type == T &&
              d.assign.type != T &&
              (d.assign = Z(Vi, d.assign, T, 0)))
          : rk(a, "default")
          ? ((null == d.assign && null == d.xa) || c(),
            Y(a),
            (d.xa = zk(a)),
            d.xa.type != N &&
              d.xa.type != T &&
              U(a, "expression following default has invalid type"),
            d.type != T && d.xa.type == T && (d.xa = Z(Ui, d.xa, N, 0)),
            d.type == T && d.xa.type != T && (d.xa = Z(Vi, d.xa, T, 0)))
          : U(a, "syntax error in parameter statement");
    }
    null != d.domain && Hk(a, d.domain);
    Y(a);
    return d;
  }
  function $k(a) {
    function b() {
      d && U(a, "at most one binary allowed");
      e.type = ah;
      d = 1;
      Y(a);
    }
    var c = 0,
      d = 0;
    a.Ob && U(a, "variable statement must precede solve statement");
    Y(a);
    a.b != Bh &&
      (sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "symbolic name missing where expected"));
    null != a.V[a.h] && U(a, a.h + " multiply declared");
    var e = {};
    e.name = a.h;
    e.Ib = null;
    e.q = 0;
    e.domain = null;
    e.type = N;
    e.P = null;
    e.W = null;
    e.O = null;
    Y(a);
    a.b == Eh && ((e.Ib = a.h), Y(a));
    a.b == ui && ((e.domain = Fk(a)), (e.q = Mk(e.domain)));
    var f = (a.V[e.name] = {});
    f.type = yh;
    for (f.link = e; ; ) {
      if (a.b == li) Y(a);
      else if (a.b == ni) break;
      if (rk(a, "integer"))
        c && U(a, "at most one integer allowed"),
          e.type != ah && (e.type = mh),
          (c = 1),
          Y(a);
      else if (rk(a, "binary")) b();
      else if (rk(a, "logical"))
        a.Ke || (ok(a, "keyword logical understood as binary"), (a.Ke = 1)),
          b();
      else if (rk(a, "symbolic")) U(a, "variable cannot be symbolic");
      else if (a.b == fi)
        null != e.P &&
          (e.P == e.W
            ? U(a, "both fixed value and lower bound not allowed")
            : U(a, "at most one lower bound allowed")),
          Y(a),
          (e.P = zk(a)),
          e.P.type == T && (e.P = Z(Ui, e.P, N, 0)),
          e.P.type != N && U(a, "expression following >= has invalid type");
      else if (a.b == di)
        null != e.W &&
          (e.W == e.P
            ? U(a, "both fixed value and upper bound not allowed")
            : U(a, "at most one upper bound allowed")),
          Y(a),
          (e.W = zk(a)),
          e.W.type == T && (e.W = Z(Ui, e.W, N, 0)),
          e.W.type != N && U(a, "expression following <= has invalid type");
      else if (a.b == ei) {
        if (null != e.P || null != e.W)
          e.P == e.W
            ? U(a, "at most one fixed value allowed")
            : null != e.P
            ? U(a, "both lower bound and fixed value not allowed")
            : U(a, "both upper bound and fixed value not allowed");
        f = a.h;
        Y(a);
        e.P = zk(a);
        e.P.type == T && (e.P = Z(Ui, e.P, N, 0));
        e.P.type != N &&
          U(a, "expression following " + f + " has invalid type");
        e.W = e.P;
      } else
        a.b == ci || a.b == gi || a.b == hi
          ? U(a, "strict bound not allowed")
          : U(a, "syntax error in variable statement");
    }
    null != e.domain && Hk(a, e.domain);
    Y(a);
    return e;
  }
  function al(a) {
    function b() {
      U(a, "syntax error in constraint statement");
    }
    var c, d, e, f;
    a.Ob && U(a, "constraint statement must precede solve statement");
    rk(a, "subject")
      ? (Y(a), rk(a, "to") || U(a, "keyword subject to incomplete"), Y(a))
      : rk(a, "subj")
      ? (Y(a), rk(a, "to") || U(a, "keyword subj to incomplete"), Y(a))
      : a.b == Th && Y(a);
    a.b != Bh &&
      (sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "symbolic name missing where expected"));
    null != a.V[a.h] && U(a, a.h + " multiply declared");
    var g = {};
    g.name = a.h;
    g.Ib = null;
    g.q = 0;
    g.domain = null;
    g.type = ch;
    g.code = null;
    g.P = null;
    g.W = null;
    g.O = null;
    Y(a);
    a.b == Eh && ((g.Ib = a.h), Y(a));
    a.b == ui && ((g.domain = Fk(a)), (g.q = Mk(g.domain)));
    c = a.V[g.name] = {};
    c.type = ch;
    c.link = g;
    a.b != mi && U(a, "colon missing where expected");
    Y(a);
    c = zk(a);
    c.type == T && (c = Z(Ui, c, N, 0));
    c.type != N &&
      c.type != jh &&
      U(a, "expression following colon has invalid type");
    a.b == li && Y(a);
    switch (a.b) {
      case di:
      case fi:
      case ei:
        break;
      case ci:
      case gi:
      case hi:
        U(a, "strict inequality not allowed");
        break;
      case ni:
        U(a, "constraint must be equality or inequality");
        break;
      default:
        b();
    }
    f = a.b;
    e = a.h;
    Y(a);
    d = zk(a);
    d.type == T && (d = Z(Ui, d, N, 0));
    d.type != N &&
      d.type != jh &&
      U(a, "expression following " + e + " has invalid type");
    a.b == li && (Y(a), a.b == ni && b());
    a.b == ci || a.b == di || a.b == ei || a.b == fi || a.b == gi || a.b == hi
      ? ((f != ei && a.b == f) ||
          U(
            a,
            "double inequality must be ... <= ... <= ... or ... >= ... >= ..."
          ),
        c.type == jh &&
          U(
            a,
            "leftmost expression in double inequality cannot be linear form"
          ),
        Y(a),
        (e = zk(a)),
        e.type == T && (e = Z(Ui, d, N, 0)),
        e.type != N &&
          e.type != jh &&
          U(
            a,
            "rightmost expression in double inequality constraint has invalid type"
          ),
        e.type == jh &&
          U(
            a,
            "rightmost expression in double inequality cannot be linear form"
          ))
      : (e = null);
    null != g.domain && Hk(a, g.domain);
    c.type != jh && (c = Z(Yi, c, jh, 0));
    d.type != jh && (d = Z(Yi, d, jh, 0));
    null != e && (e = Z(Yi, e, jh, 0));
    if (null == e)
      switch (f) {
        case di:
          g.code = c;
          g.P = null;
          g.W = d;
          break;
        case fi:
          g.code = c;
          g.P = d;
          g.W = null;
          break;
        case ei:
          (g.code = c), (g.P = d), (g.W = d);
      }
    else
      switch (f) {
        case di:
          g.code = d;
          g.P = c;
          g.W = e;
          break;
        case fi:
          (g.code = d), (g.P = e), (g.W = c);
      }
    a.b != ni && b();
    Y(a);
    return g;
  }
  function bl(a) {
    var b,
      c = { domain: null };
    c.list = b = null;
    Y(a);
    a.b == ui && (c.domain = Fk(a));
    for (a.b == mi && Y(a); ; ) {
      var d = { v: {} },
        e = function () {
          d.type = hh;
          d.v.code = Ek(a);
        };
      d.type = 0;
      d.e = null;
      null == c.list ? (c.list = d) : (b.e = d);
      b = d;
      if (a.b == Bh) {
        var f;
        Y(a);
        f = a.b;
        qk(a);
        if (f != li && f != ni) e();
        else {
          e = a.V[a.h];
          null == e && U(a, a.h + " not defined");
          d.type = e.type;
          switch (e.type) {
            case kh:
              d.v.ya = e.link;
              break;
            case uh:
              d.v.set = e.link;
              break;
            case sh:
              d.v.S = e.link;
              break;
            case yh:
              d.v.t = e.link;
              a.Ob ||
                U(
                  a,
                  "invalid reference to variable " +
                    d.v.t.name +
                    " above solve statement"
                );
              break;
            case ch:
              (d.v.H = e.link),
                a.Ob ||
                  U(
                    a,
                    "invalid reference to " +
                      (d.v.H.type == ch ? "constraint" : "objective") +
                      " " +
                      d.v.H.name +
                      " above solve statement"
                  );
          }
          Y(a);
        }
      } else e();
      if (a.b == li) Y(a);
      else break;
    }
    null != c.domain && Hk(a, c.domain);
    a.b != ni && U(a, "syntax error in display statement");
    Y(a);
    return c;
  }
  function cl(a) {
    (!a.nc && rk(a, "end")) || (a.nc && dl(a, "end"))
      ? (Y(a),
        a.b == ni
          ? Y(a)
          : ok(
              a,
              "no semicolon following end statement; missing semicolon inserted"
            ))
      : ok(a, "unexpected end of file; missing end statement inserted");
    a.b != Ah && ok(a, "some text detected beyond end statement; text ignored");
  }
  function el(a, b) {
    var c = { v: {} };
    c.bb = a.bb;
    c.Vc = a.Vc;
    c.e = null;
    if (rk(a, "set"))
      b && U(a, "set statement not allowed here"),
        (c.type = uh),
        (c.v.set = Yk(a));
    else if (rk(a, "param"))
      b && U(a, "parameter statement not allowed here"),
        (c.type = sh),
        (c.v.S = Zk(a));
    else if (rk(a, "var"))
      b && U(a, "variable statement not allowed here"),
        (c.type = yh),
        (c.v.t = $k(a));
    else if (rk(a, "subject") || rk(a, "subj") || a.b == Th)
      b && U(a, "constraint statement not allowed here"),
        (c.type = ch),
        (c.v.H = al(a));
    else if (rk(a, "minimize") || rk(a, "maximize")) {
      b && U(a, "objective statement not allowed here");
      c.type = ch;
      var d = c.v,
        e,
        f;
      rk(a, "minimize") ? (f = ph) : rk(a, "maximize") && (f = oh);
      a.Ob && U(a, "objective statement must precede solve statement");
      Y(a);
      a.b != Bh &&
        (sk(a)
          ? U(a, "invalid use of reserved keyword " + a.h)
          : U(a, "symbolic name missing where expected"));
      null != a.V[a.h] && U(a, a.h + " multiply declared");
      e = {};
      e.name = a.h;
      e.Ib = null;
      e.q = 0;
      e.domain = null;
      e.type = f;
      e.code = null;
      e.P = null;
      e.W = null;
      e.O = null;
      Y(a);
      a.b == Eh && ((e.Ib = a.h), Y(a));
      a.b == ui && ((e.domain = Fk(a)), (e.q = Mk(e.domain)));
      f = a.V[e.name] = {};
      f.type = ch;
      f.link = e;
      a.b != mi && U(a, "colon missing where expected");
      Y(a);
      e.code = zk(a);
      e.code.type == T && (e.code = Z(Ui, e.code, N, 0));
      e.code.type == N && (e.code = Z(Yi, e.code, jh, 0));
      e.code.type != jh && U(a, "expression following colon has invalid type");
      null != e.domain && Hk(a, e.domain);
      a.b != ni && U(a, "syntax error in objective statement");
      Y(a);
      d.H = e;
    } else if (rk(a, "table")) {
      b && U(a, "table statement not allowed here");
      c.type = wh;
      var d = c.v,
        g,
        h,
        k;
      Y(a);
      a.b != Bh &&
        (sk(a)
          ? U(a, "invalid use of reserved keyword " + a.h)
          : U(a, "symbolic name missing where expected"));
      null != a.V[a.h] && U(a, a.h + " multiply declared");
      e = { v: { qa: {}, Oc: {} } };
      e.name = a.h;
      Y(a);
      a.b == Eh ? ((e.Ib = a.h), Y(a)) : (e.Ib = null);
      a.b == ui
        ? ((e.type = rh),
          (e.v.Oc.domain = Fk(a)),
          rk(a, "OUT") || U(a, "keyword OUT missing where expected"))
        : ((e.type = lh),
          rk(a, "IN") || U(a, "keyword IN missing where expected"));
      Y(a);
      for (e.a = f = null; ; )
        if (
          ((g = {}),
          (a.b != li && a.b != mi && a.b != ni) ||
            U(a, "argument expression missing where expected"),
          (g.code = zk(a)),
          g.code.type == N && (g.code = Z(Vi, g.code, T, 0)),
          g.code.type != T && U(a, "argument expression has invalid type"),
          (g.e = null),
          null == f ? (e.a = g) : (f.e = g),
          (f = g),
          a.b == li)
        )
          Y(a);
        else if (a.b == mi || a.b == ni) break;
      a.b == mi ? Y(a) : U(a, "colon missing where expected");
      switch (e.type) {
        case lh:
          a.b == Bh
            ? ((g = a.V[a.h]),
              null == g && U(a, a.h + " not defined"),
              g.type != uh && U(a, a.h + " not a set"),
              (e.v.qa.set = g.link),
              null != e.v.qa.set.assign && U(a, a.h + " needs no data"),
              0 != e.v.qa.set.q && U(a, a.h + " must be a simple set"),
              Y(a),
              a.b == yi ? Y(a) : U(a, "delimiter <- missing where expected"))
            : sk(a)
            ? U(a, "invalid use of reserved keyword " + a.h)
            : (e.v.qa.set = null);
          e.v.qa.We = g = null;
          f = 0;
          for (a.b == si ? Y(a) : U(a, "field list missing where expected"); ; )
            if (
              ((h = {}),
              a.b != Bh &&
                (sk(a)
                  ? U(a, "invalid use of reserved keyword " + a.h)
                  : U(a, "field name missing where expected")),
              (h.name = a.h),
              Y(a),
              (h.e = null),
              null == g ? (e.v.qa.We = h) : (g.e = h),
              (g = h),
              f++,
              a.b == li)
            )
              Y(a);
            else if (a.b == ti) break;
            else U(a, "syntax error in field list");
          null != e.v.qa.set &&
            e.v.qa.set.X != f &&
            U(
              a,
              "there must be " +
                e.v.qa.set.X +
                " field" +
                (1 == e.v.qa.set.X ? "" : "s") +
                " rather than " +
                f
            );
          Y(a);
          for (e.v.qa.list = h = null; a.b == li; )
            Y(a),
              (k = {}),
              a.b != Bh &&
                (sk(a)
                  ? U(a, "invalid use of reserved keyword " + a.h)
                  : U(a, "parameter name missing where expected")),
              (g = a.V[a.h]),
              null == g && U(a, a.h + " not defined"),
              g.type != sh && U(a, a.h + " not a parameter"),
              (k.S = g.link),
              k.S.q != f &&
                U(
                  a,
                  a.h +
                    " must have " +
                    f +
                    " subscript" +
                    (1 == f ? "" : "s") +
                    " rather than " +
                    k.S.q
                ),
              null != k.S.assign && U(a, a.h + " needs no data"),
              Y(a),
              a.b == xi
                ? (Y(a),
                  a.b != Bh &&
                    (sk(a)
                      ? U(a, "invalid use of reserved keyword " + a.h)
                      : U(a, "field name missing where expected")),
                  (g = a.h),
                  Y(a))
                : (g = k.S.name),
              (k.name = g),
              (k.e = null),
              null == h ? (e.v.qa.list = k) : (h.e = k),
              (h = k);
          break;
        case rh:
          for (e.v.Oc.list = f = null; ; )
            if (
              ((h = {}),
              (a.b != li && a.b != ni) ||
                U(a, "expression missing where expected"),
              (g = a.b == Bh ? a.h : ""),
              (h.code = zk(a)),
              a.b == xi &&
                (Y(a),
                a.b != Bh &&
                  (sk(a)
                    ? U(a, "invalid use of reserved keyword " + a.h)
                    : U(a, "field name missing where expected")),
                (g = a.h),
                Y(a)),
              "" == g && U(a, "field name required"),
              (h.name = g),
              (h.e = null),
              null == f ? (e.v.Oc.list = h) : (f.e = h),
              (f = h),
              a.b == li)
            )
              Y(a);
            else if (a.b == ni) break;
            else U(a, "syntax error in output list");
          Hk(a, e.v.Oc.domain);
      }
      a.b != ni && U(a, "syntax error in table statement");
      Y(a);
      d.nd = e;
    } else if (rk(a, "solve"))
      b && U(a, "solve statement not allowed here"),
        (c.type = vh),
        (d = c.v),
        a.Ob && U(a, "at most one solve statement allowed"),
        (a.Ob = 1),
        Y(a),
        a.b != ni && U(a, "syntax error in solve statement"),
        Y(a),
        (d.zh = null);
    else if (rk(a, "check"))
      (c.type = bh),
        (d = c.v),
        (e = { domain: null, code: null }),
        Y(a),
        a.b == ui && (e.domain = Fk(a)),
        a.b == mi && Y(a),
        (e.code = Ek(a)),
        e.code.type != nh && U(a, "expression has invalid type"),
        null != e.domain && Hk(a, e.domain),
        a.b != ni && U(a, "syntax error in check statement"),
        Y(a),
        (d.Rg = e);
    else if (rk(a, "display")) (c.type = dh), (c.v.Sg = bl(a));
    else if (rk(a, "printf")) {
      c.type = th;
      d = c.v;
      g = { domain: null, zd: null };
      g.list = f = null;
      Y(a);
      a.b == ui && (g.domain = Fk(a));
      a.b == mi && Y(a);
      g.zd = zk(a);
      g.zd.type == N && (g.zd = Z(Vi, g.zd, T, 0));
      for (
        g.zd.type != T && U(a, "format expression has invalid type");
        a.b == li;

      )
        Y(a),
          (e = { code: null, e: null }),
          null == g.list ? (g.list = e) : (f.e = e),
          (f = e),
          (e.code = Gk(a)),
          e.code.type != N &&
            e.code.type != T &&
            e.code.type != nh &&
            U(a, "only numeric, symbolic, or logical expression allowed");
      null != g.domain && Hk(a, g.domain);
      g.Ea = null;
      g.Mg = 0;
      if (a.b == gi || a.b == wi)
        (g.Mg = a.b == wi),
          Y(a),
          (g.Ea = zk(a)),
          g.Ea.type == N && (g.Ea = Z(Vi, g.Ea, T, 0)),
          g.Ea.type != T && U(a, "file name expression has invalid type");
      a.b != ni && U(a, "syntax error in printf statement");
      Y(a);
      d.nh = g;
    } else if (rk(a, "for")) {
      c.type = ih;
      d = c.v;
      g = { domain: null };
      g.list = f = null;
      Y(a);
      a.b != ui && U(a, "indexing expression missing where expected");
      g.domain = Fk(a);
      a.b == mi && Y(a);
      if (a.b != ui) g.list = el(a, 1);
      else {
        for (Y(a); a.b != vi; )
          (e = el(a, 1)), null == f ? (g.list = e) : (f.e = e), (f = e);
        Y(a);
      }
      Hk(a, g.domain);
      d.Ug = g;
    } else
      a.b == Bh
        ? (b && U(a, "constraint statement not allowed here"),
          (c.type = ch),
          (c.v.H = al(a)))
        : sk(a)
        ? U(a, "invalid use of reserved keyword " + a.h)
        : U(a, "syntax error in model section");
    return c;
  }
  function fl(a) {
    var b, c;
    for (c = null; a.b != Ah && !rk(a, "data") && !rk(a, "end"); )
      (b = el(a, 0)), null == c ? (a.uc = b) : (c.e = b), (c = b);
  }
  function gl(a, b) {
    var c,
      d = {};
    d.Y = b;
    d.e = null;
    if (null == a) a = d;
    else {
      for (c = a; null != c.e; c = c.e);
      c.e = d;
    }
    return a;
  }
  function hl(a) {
    for (var b = 0; null != a; a = a.e) b++;
    return b;
  }
  function il(a) {
    for (var b = 0; null != a; a = a.e) null == a.Y && b++;
    return b;
  }
  function jl(a) {
    for (var b = null; 0 < a--; ) b = gl(b, null);
    return b;
  }
  function kl(a) {
    return a.b == Dh || a.b == Ch || a.b == Eh;
  }
  function dl(a, b) {
    return kl(a) && a.h == b;
  }
  function ll(a) {
    var b;
    b = a.b == Dh ? ml(a.value) : nl(a.h);
    Y(a);
    return b;
  }
  function ol(a, b, c) {
    var d, e;
    switch (a.b) {
      case si:
        e = ti;
        break;
      case qi:
        e = ri;
    }
    0 == c && U(a, b + " cannot be subscripted");
    Y(a);
    for (d = null; ; )
      if (
        (kl(a)
          ? (d = gl(d, ll(a)))
          : a.b == $h
          ? ((d = gl(d, null)), Y(a))
          : U(a, "number, symbol, or asterisk missing where expected"),
        a.b == li)
      )
        Y(a);
      else if (a.b == e) break;
      else U(a, "syntax error in slice");
    if (hl(d) != c)
      switch (e) {
        case ti:
          U(
            a,
            b +
              " must have " +
              c +
              " subscript" +
              (1 == c ? "" : "s") +
              ", not " +
              hl(d)
          );
          break;
        case ri:
          U(a, b + " has dimension " + c + ", not " + hl(d));
      }
    Y(a);
    return d;
  }
  function pl(a, b) {
    var c;
    c = a.V[b];
    (null != c && c.type == uh) || U(a, b + " not a set");
    c = c.link;
    (null == c.assign && null == c.Yc) || U(a, b + " needs no data");
    c.data = 1;
    return c;
  }
  function ql(a, b, c) {
    var d,
      e,
      f = null;
    for (d = null; null != c; c = c.e)
      null == c.Y
        ? (kl(a) ||
            ((e = il(c)),
            1 == e
              ? U(a, "one item missing in data group beginning with " + rl(f))
              : U(
                  a,
                  e + " items missing in data group beginning with " + rl(f)
                )),
          (e = ll(a)),
          null == f && (f = e))
        : (e = sl(c.Y)),
        (d = tl(d, e)),
        null != c.e && a.b == li && Y(a);
    ul(a, b.value.set, d);
  }
  function vl(a, b, c, d) {
    var e, f, g, h, k;
    for (e = null; a.b != oi; )
      kl(a) || U(a, "number, symbol, or := missing where expected"),
        (e = gl(e, ll(a)));
    for (Y(a); kl(a); )
      for (k = ll(a), f = e; null != f; f = f.e) {
        var l = 0;
        if (!dl(a, "+"))
          if (dl(a, "-")) {
            Y(a);
            continue;
          } else
            (g = hl(f)),
              1 == g
                ? U(a, "one item missing in data group beginning with " + rl(k))
                : U(
                    a,
                    g + " items missing in data group beginning with " + rl(k)
                  );
        h = null;
        for (g = c; null != g; g = g.e)
          if (null == g.Y)
            switch (++l) {
              case 1:
                h = tl(h, sl(d ? f.Y : k));
                break;
              case 2:
                h = tl(h, sl(d ? k : f.Y));
            }
          else h = tl(h, sl(g.Y));
        ul(a, b.value.set, h);
        Y(a);
      }
  }
  function wl(a) {
    function b() {
      U(a, "slice currently used must specify 2 asterisks, not " + il(h));
    }
    function c() {
      U(a, "transpose indicator (tr) incomplete");
    }
    function d() {
      Y(a);
      dl(a, "tr") || c();
      2 != il(h) && b();
      Y(a);
      a.b != ri && c();
      Y(a);
      a.b == mi && Y(a);
      k = 1;
      vl(a, g, h, k);
    }
    var e,
      f,
      g,
      h,
      k = 0;
    Y(a);
    kl(a) || U(a, "set name missing where expected");
    e = pl(a, a.h);
    Y(a);
    f = null;
    if (a.b == si) {
      0 == e.q && U(a, e.name + " cannot be subscripted");
      for (Y(a); ; )
        if (
          (kl(a) || U(a, "number or symbol missing where expected"),
          (f = tl(f, ll(a))),
          a.b == li)
        )
          Y(a);
        else if (a.b == ti) break;
        else U(a, "syntax error in subscript list");
      e.q != xl(f) &&
        U(
          a,
          e.name +
            " must have " +
            e.q +
            " subscript" +
            (1 == e.q ? "" : "s") +
            " rather than " +
            xl(f)
        );
      Y(a);
    } else 0 != e.q && U(a, e.name + " must be subscripted");
    null != yl(a, e.O, f) && U(a, e.name + zl("[", f) + " already defined");
    g = Al(e.O, f);
    g.value.set = Bl(a, qh, e.X);
    for (h = jl(e.X); ; )
      if ((a.b == li && Y(a), a.b == oi)) Y(a);
      else if (a.b == qi)
        Y(a),
          (f = dl(a, "tr")),
          qk(a),
          f
            ? d()
            : ((h = ol(a, e.name, e.X)), (k = 0), 0 == il(h) && ql(a, g, h));
      else if (kl(a)) ql(a, g, h);
      else if (a.b == mi) 2 != il(h) && b(), Y(a), vl(a, g, h, k);
      else if (a.b == qi) d();
      else if (a.b == ni) {
        Y(a);
        break;
      } else U(a, "syntax error in set data block");
  }
  function Cl(a, b) {
    var c;
    c = a.V[b];
    (null != c && c.type == sh) || U(a, b + " not a parameter");
    c = c.link;
    null != c.assign && U(a, b + " needs no data");
    c.data && U(a, b + " already provided with data");
    c.data = 1;
    return c;
  }
  function Dl(a, b, c) {
    null != b.xa &&
      U(
        a,
        "default value for " + b.name + " already specified in model section"
      );
    b.Wc = c;
  }
  function El(a, b, c) {
    null != yl(a, b.O, c) && U(a, b.name + zl("[", c) + " already defined");
    c = Al(b.O, c);
    switch (b.type) {
      case N:
      case mh:
      case ah:
        a.b == Dh || U(a, b.name + " requires numeric data");
        b = c.value;
        c = a.value;
        Y(a);
        b.Q = c;
        break;
      case T:
        c.value.Y = ll(a);
    }
  }
  function Fl(a, b, c) {
    var d,
      e,
      f = null;
    for (d = null; null != c; c = c.e)
      null == c.Y
        ? (kl(a) ||
            U(
              a,
              il(c) + 1 + " items missing in data group beginning with " + rl(f)
            ),
          (e = ll(a)),
          null == f && (f = e))
        : (e = sl(c.Y)),
        (d = tl(d, e)),
        a.b == li && Y(a);
    kl(a) || U(a, "one item missing in data group beginning with " + rl(f));
    El(a, b, d);
  }
  function Gl(a, b, c, d) {
    var e, f, g, h, k;
    for (e = null; a.b != oi; )
      kl(a) || U(a, "number, symbol, or := missing where expected"),
        (e = gl(e, ll(a)));
    for (Y(a); kl(a); )
      for (k = ll(a), f = e; null != f; f = f.e) {
        var l = 0;
        if (dl(a, ".")) Y(a);
        else {
          h = null;
          for (g = c; null != g; g = g.e)
            if (null == g.Y)
              switch (++l) {
                case 1:
                  h = tl(h, sl(d ? f.Y : k));
                  break;
                case 2:
                  h = tl(h, sl(d ? k : f.Y));
              }
            else h = tl(h, sl(g.Y));
          kl(a) ||
            ((g = hl(f)),
            1 == g
              ? U(a, "one item missing in data group beginning with " + rl(k))
              : U(
                  a,
                  g + " items missing in data group beginning with " + rl(k)
                ));
          El(a, b, h);
        }
      }
  }
  function Hl(a, b) {
    var c = null,
      d,
      e,
      f,
      g,
      h = 0;
    g = null;
    kl(a) &&
      (Y(a),
      (e = a.b),
      qk(a),
      e == mi &&
        ((c = pl(a, a.h)),
        0 != c.q && U(a, c.name + " must be a simple set"),
        null != c.O.head && U(a, c.name + " already defined"),
        (Al(c.O, null).value.set = Bl(a, qh, c.X)),
        (g = c.name),
        (h = c.X),
        Y(a),
        Y(a)));
    for (e = null; a.b != oi; )
      kl(a) || U(a, "parameter name or := missing where expected"),
        (d = Cl(a, a.h)),
        0 == d.q && U(a, a.h + " not a subscripted parameter"),
        0 != h &&
          d.q != h &&
          U(
            a,
            g +
              " has dimension " +
              h +
              " while " +
              d.name +
              " has dimension " +
              d.q
          ),
        null != b && Dl(a, d, sl(b)),
        (e = gl(e, d)),
        (g = d.name),
        (h = d.q),
        Y(a),
        a.b == li && Y(a);
    0 == hl(e) && U(a, "at least one parameter name required");
    Y(a);
    for (a.b == li && Y(a); kl(a); ) {
      g = null;
      for (f = 1; f <= h; f++)
        kl(a) ||
          ((d = hl(e) + h - f + 1),
          U(a, d + " items missing in data group beginning with " + rl(g.Y))),
          (g = tl(g, ll(a))),
          f < h && a.b == li && Y(a);
      null != c && ul(a, c.O.head.value.set, Il(g));
      a.b == li && Y(a);
      for (f = e; null != f; f = f.e)
        dl(a, ".")
          ? Y(a)
          : (kl(a) ||
              ((d = hl(f)),
              1 == d
                ? U(
                    a,
                    "one item missing in data group beginning with " + rl(g.Y)
                  )
                : U(
                    a,
                    d + " items missing in data group beginning with " + rl(g.Y)
                  )),
            El(a, f.Y, Il(g)),
            null != f.e && a.b == li && Y(a));
      a.b == li && (Y(a), kl(a) || qk(a));
    }
    for (f = e; null != f; f = f.e) f.Y = null;
  }
  function Jl(a) {
    function b() {
      U(a, e.name + " not a subscripted parameter");
    }
    function c() {
      U(a, "slice currently used must specify 2 asterisks, not " + il(g));
    }
    function d() {
      U(a, "transpose indicator (tr) incomplete");
    }
    var e,
      f = null,
      g,
      h = 0;
    Y(a);
    dl(a, "default") &&
      (Y(a),
      kl(a) || U(a, "default value missing where expected"),
      (f = ll(a)),
      a.b != mi && U(a, "colon missing where expected"));
    if (a.b == mi)
      Y(a),
        a.b == li && Y(a),
        Hl(a, f),
        a.b != ni &&
          U(a, "symbol, number, or semicolon missing where expected"),
        Y(a);
    else
      for (
        kl(a) || U(a, "parameter name missing where expected"),
          e = Cl(a, a.h),
          Y(a),
          dl(a, "default") &&
            (Y(a),
            kl(a) || U(a, "default value missing where expected"),
            (f = ll(a)),
            Dl(a, e, f)),
          g = jl(e.q);
        ;

      )
        if ((a.b == li && Y(a), a.b == oi)) Y(a);
        else if (a.b == si) (g = ol(a, e.name, e.q)), (h = 0);
        else if (kl(a)) Fl(a, e, g);
        else if (a.b == mi)
          0 == e.q && b(), 2 != il(g) && c(), Y(a), Gl(a, e, g, h);
        else if (a.b == qi)
          Y(a),
            dl(a, "tr") || d(),
            0 == e.q && b(),
            2 != il(g) && c(),
            Y(a),
            a.b != ri && d(),
            Y(a),
            a.b == mi && Y(a),
            (h = 1),
            Gl(a, e, g, h);
        else if (a.b == ni) {
          Y(a);
          break;
        } else U(a, "syntax error in parameter data block");
  }
  function Kl(a) {
    for (; a.b != Ah && !dl(a, "end"); )
      dl(a, "set")
        ? wl(a)
        : dl(a, "param")
        ? Jl(a)
        : U(a, "syntax error in data section");
  }
  function Ll(a, b, c) {
    ((0 < b && 0 < c && b > 0.999 * s - c) ||
      (0 > b && 0 > c && b < -0.999 * s - c)) &&
      U(a, b + " + " + c + "; floating-point overflow");
    return b + c;
  }
  function Ml(a, b, c) {
    ((0 < b && 0 > c && b > 0.999 * s + c) ||
      (0 > b && 0 < c && b < -0.999 * s + c)) &&
      U(a, b + " - " + c + "; floating-point overflow");
    return b - c;
  }
  function Nl(a, b, c) {
    if (b < c) return 0;
    0 < b &&
      0 > c &&
      b > 0.999 * s + c &&
      U(a, b + " less " + c + "; floating-point overflow");
    return b - c;
  }
  function Ol(a, b, c) {
    1 < Math.abs(c) &&
      Math.abs(b) > (0.999 * s) / Math.abs(c) &&
      U(a, b + " * " + c + "; floating-point overflow");
    return b * c;
  }
  function Pl(a, b, c) {
    Math.abs(c) < aa && U(a, b + " / " + c + "; floating-point zero divide");
    1 > Math.abs(c) &&
      Math.abs(b) > 0.999 * s * Math.abs(c) &&
      U(a, b + " / " + c + "; floating-point overflow");
    return b / c;
  }
  function Ql(a, b, c) {
    Math.abs(c) < aa && U(a, b + " div " + c + "; floating-point zero divide");
    1 > Math.abs(c) &&
      Math.abs(b) > 0.999 * s * Math.abs(c) &&
      U(a, b + " div " + c + "; floating-point overflow");
    b /= c;
    return 0 < b ? Math.floor(b) : 0 > b ? Math.ceil(b) : 0;
  }
  function Rl(a, b) {
    var c;
    if (0 == a) c = 0;
    else if (0 == b) c = a;
    else if (
      ((c = Math.abs(a) % Math.abs(b)),
      0 != c && (0 > a && (c = -c), (0 < a && 0 > b) || (0 > a && 0 < b)))
    )
      c += b;
    return c;
  }
  function Sl(a, b, c) {
    ((0 == b && 0 >= c) || (0 > b && c != Math.floor(c))) &&
      U(a, b + " ** " + c + "; result undefined");
    0 == b
      ? (a = Math.pow(b, c))
      : (((1 < Math.abs(b) &&
          1 < c &&
          +Math.log(Math.abs(b)) > (0.999 * Math.log(s)) / c) ||
          (1 > Math.abs(b) &&
            -1 > c &&
            +Math.log(Math.abs(b)) < (0.999 * Math.log(s)) / c)) &&
          U(a, b + " ** " + c + "; floating-point overflow"),
        (a =
          (1 < Math.abs(b) &&
            -1 > c &&
            -Math.log(Math.abs(b)) < (0.999 * Math.log(s)) / c) ||
          (1 > Math.abs(b) &&
            1 < c &&
            -Math.log(Math.abs(b)) > (0.999 * Math.log(s)) / c)
            ? 0
            : Math.pow(b, c)));
    return a;
  }
  function Tl(a, b) {
    b > 0.999 * Math.log(s) && U(a, "exp(" + b + "); floating-point overflow");
    return Math.exp(b);
  }
  function Ul(a, b) {
    0 >= b && U(a, "log(" + b + "); non-positive argument");
    return Math.log(b);
  }
  function Vl(a, b) {
    0 >= b && U(a, "log10(" + b + "); non-positive argument");
    return Math.log(b) / Math.LN10;
  }
  function Wl(a, b) {
    0 > b && U(a, "sqrt(" + b + "); negative argument");
    return Math.sqrt(b);
  }
  function Xl(a, b) {
    (-1e6 <= b && 1e6 >= b) || U(a, "sin(" + b + "); argument too large");
    return Math.sin(b);
  }
  function Yl(a, b) {
    (-1e6 <= b && 1e6 >= b) || U(a, "cos(" + b + "); argument too large");
    return Math.cos(b);
  }
  function Zl(a) {
    return Math.atan(a);
  }
  function $l(a, b) {
    return Math.atan2(a, b);
  }
  function am(a, b, c) {
    c != Math.floor(c) &&
      U(a, "round(" + b + ", " + c + "); non-integer second argument");
    18 >= c &&
      ((a = Math.pow(10, c)),
      Math.abs(b) < (0.999 * s) / a &&
        ((b = Math.floor(b * a + 0.5)), 0 != b && (b /= a)));
    return b;
  }
  function bm(a, b, c) {
    c != Math.floor(c) &&
      U(a, "trunc(" + b + ", " + c + "); non-integer second argument");
    18 >= c &&
      ((a = Math.pow(10, c)),
      Math.abs(b) < (0.999 * s) / a &&
        ((b = 0 <= b ? Math.floor(b * a) : Math.ceil(b * a)),
        0 != b && (b /= a)));
    return b;
  }
  function cm(a, b, c) {
    var d;
    b >= c && U(a, "Uniform(" + b + ", " + c + "); invalid range");
    d = dm(a.Gd) / 2147483648;
    return (d = Ll(a, b * (1 - d), c * d));
  }
  function em(a) {
    var b, c;
    do
      (b = -1 + 2 * (dm(a.Gd) / 2147483648)),
        (c = -1 + 2 * (dm(a.Gd) / 2147483648)),
        (b = b * b + c * c);
    while (1 < b || 0 == b);
    return c * Math.sqrt((-2 * Math.log(b)) / b);
  }
  function fm(a, b, c) {
    return Ll(a, b, Ol(a, c, em(a)));
  }
  function ml(a) {
    var b = {};
    b.Q = a;
    b.M = null;
    return b;
  }
  function nl(a) {
    var b = { Q: 0 };
    b.M = a;
    return b;
  }
  function sl(a) {
    var b = {};
    null == a.M ? ((b.Q = a.Q), (b.M = null)) : ((b.Q = 0), (b.M = a.M));
    return b;
  }
  function gm(a, b) {
    var c;
    if (null == a.M && null == b.M) c = a.Q < b.Q ? -1 : a.Q > b.Q ? 1 : 0;
    else if (null == a.M) c = -1;
    else if (null == b.M) c = 1;
    else {
      c = a.M;
      var d = b.M;
      c = c == d ? 0 : c > d ? 1 : -1;
    }
    return c;
  }
  function rl(a) {
    var b;
    if (null == a.M) b = String(a.Q);
    else {
      var c,
        d,
        e = a.M;
      if (ua(e[0]) || "_" == e[0])
        for (a = !1, c = 1; c < e.length; c++) {
          if (!(va(e[c]) || 0 <= "+-._".indexOf(e[c]))) {
            a = !0;
            break;
          }
        }
      else a = !0;
      b = "";
      d = 0;
      var f = function (a) {
        255 > d && ((b += a), d++);
      };
      a && f("'");
      for (c = 0; c < e.length; c++) a && "'" == e[c] && f("'"), f(e[c]);
      a && f("'");
      255 == d && (b = b.slice(0, 252) + "...");
    }
    return b;
  }
  function tl(a, b) {
    var c,
      d = {};
    d.Y = b;
    d.e = null;
    if (null == a) a = d;
    else {
      for (c = a; null != c.e; c = c.e);
      c.e = d;
    }
    return a;
  }
  function xl(a) {
    for (var b = 0; null != a; a = a.e) b++;
    return b;
  }
  function Il(a) {
    var b, c;
    if (null == a) b = null;
    else {
      for (b = c = {}; null != a; a = a.e)
        (c.Y = sl(a.Y)), null != a.e && (c = c.e = {});
      c.e = null;
    }
    return b;
  }
  function hm(a, b) {
    var c, d, e;
    c = a;
    for (d = b; null != c; c = c.e, d = d.e)
      if (((e = gm(c.Y, d.Y)), 0 != e)) return e;
    return 0;
  }
  function im(a, b) {
    for (var c = null, d = 1, e = a; d <= b; d++, e = e.e) c = tl(c, sl(e.Y));
    return c;
  }
  function zl(a, b) {
    function c(a) {
      255 > f && (g += a);
      f++;
    }
    var d,
      e,
      f = 0,
      g = "",
      h = "",
      k = xl(b);
    "[" == a && 0 < k && c("[");
    "(" == a && 1 < k && c("(");
    for (d = b; null != d; d = d.e)
      for (d != b && c(","), h = rl(d.Y), e = 0; e < h.length; e++) c(h[e]);
    "[" == a && 0 < k && c("]");
    "(" == a && 1 < k && c(")");
    255 == f && (g = g.slice(0, 252) + "...");
    return g;
  }
  function jm(a, b) {
    Al(a, b).value.eh = null;
  }
  function ul(a, b, c) {
    null != yl(a, b, c) && U(a, "duplicate tuple " + zl("(", c) + " detected");
    jm(b, c);
  }
  function km(a, b) {
    var c, d;
    c = Bl(a, qh, b.q);
    for (d = b.head; null != d; d = d.e) jm(c, Il(d.w));
    return c;
  }
  function lm(a, b, c, d) {
    var e;
    0 == d && U(a, b + " .. " + c + " by " + d + "; zero stride not allowed");
    e =
      0 < c && 0 > b && c > 0.999 * s + b
        ? +s
        : 0 > c && 0 < b && c < -0.999 * s + b
        ? -s
        : c - b;
    1 > Math.abs(d) && Math.abs(e) > 0.999 * s * Math.abs(d)
      ? (e = (0 < e && 0 < d) || (0 > e && 0 > d) ? +s : 0)
      : ((e = Math.floor(e / d) + 1), 0 > e && (e = 0));
    2147483646 < e && U(a, b + " .. " + c + " by " + d + "; set too large");
    return (e + 0.5) | 0;
  }
  function mm(a, b, c, d, e) {
    1 <= e && lm(a, b, c, d);
    return b + (e - 1) * d;
  }
  function nm(a) {
    var b;
    0 == a ? (b = null) : ((b = {}), (b.u = a), (b.t = null), (b.e = null));
    return b;
  }
  function om(a) {
    var b, c;
    if (null == a) b = null;
    else {
      for (b = c = {}; null != a; a = a.e)
        (c.u = a.u), (c.t = a.t), null != a.e && (c = c.e = {});
      c.e = null;
    }
    return b;
  }
  function pm(a, b, c, d, e) {
    var f = null,
      g,
      h = 0;
    for (g = c; null != g; g = g.e)
      null == g.t
        ? (h = Ll(a, h, Ol(a, b, g.u)))
        : (g.t.ja = Ll(a, g.t.ja, Ol(a, b, g.u)));
    for (g = e; null != g; g = g.e)
      null == g.t
        ? (h = Ll(a, h, Ol(a, d, g.u)))
        : (g.t.ja = Ll(a, g.t.ja, Ol(a, d, g.u)));
    for (g = c; null != g; g = g.e)
      null != g.t &&
        0 != g.t.ja &&
        ((a = {}),
        (a.u = g.t.ja),
        (a.t = g.t),
        (a.e = f),
        (f = a),
        (g.t.ja = 0));
    for (g = e; null != g; g = g.e)
      null != g.t &&
        0 != g.t.ja &&
        ((a = {}),
        (a.u = g.t.ja),
        (a.t = g.t),
        (a.e = f),
        (f = a),
        (g.t.ja = 0));
    0 != h && ((a = {}), (a.u = h), (a.t = null), (a.e = f), (f = a));
    return f;
  }
  function qm(a, b, c) {
    for (var d = null, e, f = 0; null != b; )
      (e = b),
        (b = b.e),
        null == e.t ? (f = Ll(a, f, e.u)) : ((e.e = d), (d = e));
    c(f);
    return d;
  }
  function rm(a, b) {
    switch (a) {
      case qh:
        b.eh = null;
        break;
      case N:
        b.Q = 0;
        break;
      case T:
        b.Y = null;
        break;
      case nh:
        b.sg = 0;
        break;
      case xh:
        b.w = null;
        break;
      case fh:
        b.set = null;
        break;
      case gh:
        b.t = null;
        break;
      case jh:
        b.form = null;
        break;
      case eh:
        b.H = null;
    }
  }
  function Bl(a, b, c) {
    var d = {};
    d.type = b;
    d.q = c;
    d.size = 0;
    d.head = null;
    d.Xa = null;
    d.V = !1;
    d.ca = null;
    d.e = a.pg;
    null != d.e && (d.e.ca = d);
    return (a.pg = d);
  }
  function sm(a, b, c) {
    return hm(b, c);
  }
  function yl(a, b, c) {
    if (30 < b.size && !b.V) {
      var d = { root: null };
      d.yg = sm;
      d.info = a;
      d.size = 0;
      d.height = 0;
      b.V = d;
      for (a = b.head; null != a; a = a.e) qe(b.V, a.w).link = a;
    }
    if (b.V) {
      b = b.V;
      for (a = b.root; null != a; ) {
        d = b.yg(b.info, c, a.key);
        if (0 == d) break;
        a = 0 > d ? a.left : a.right;
      }
      c = a;
      a = null == c ? null : c.link;
    } else for (a = b.head; null != a && 0 != hm(a.w, c); a = a.e);
    return a;
  }
  function Al(a, b) {
    var c = {};
    c.w = b;
    c.e = null;
    c.value = {};
    a.size++;
    null == a.head ? (a.head = c) : (a.Xa.e = c);
    a.Xa = c;
    null != a.V && (qe(a.V, c.w).link = c);
    return c;
  }
  function tm(a) {
    var b;
    if (null != a.Le)
      for (b = a.list, a = a.Le; null != b; b = b.e, a = a.e)
        a: {
          var c = b,
            d = a.Y,
            e = void 0,
            f = void 0;
          if (null != c.value) {
            if (0 == gm(c.value, d)) break a;
            c.value = null;
          }
          for (e = c.list; null != e; e = e.a.index.e)
            for (f = e; null != f; f = f.R)
              f.valid && ((f.valid = 0), rm(f.type, f.value));
          c.value = sl(d);
        }
  }
  function um(a, b, c, d, e) {
    var f,
      g = 0;
    if (!vm(a, b.code, c)) return 1;
    f = b.Le;
    b.Le = c;
    tm(b);
    e(a, d);
    b.Le = f;
    tm(b);
    return g;
  }
  function wm(a, b) {
    if (null != b.Kc) {
      var c,
        d,
        e = null,
        f = null;
      c = b.Kc;
      b.Kc = c.e;
      for (d = c.list; null != d; d = d.e)
        null == e ? (e = f = {}) : (f = f.e = {}),
          null == d.code
            ? ((f.Y = b.w.Y), (b.w = b.w.e))
            : (f.Y = xm(a, d.code));
      f.e = null;
      um(a, c, e, b, wm) && (b.Ve = 1);
      for (d = c.list; null != d; d = d.e) e = e.e;
    } else
      null == b.domain.code || ym(a, b.domain.code)
        ? b.Zd(a, b.info)
        : (b.Ve = 2);
  }
  function zm(a, b, c, d, e) {
    var f = {};
    null == b
      ? (e(a, d), (f.Ve = 0))
      : ((f.domain = b),
        (f.Kc = b.list),
        (f.w = c),
        (f.info = d),
        (f.Zd = e),
        (f.Ve = 0),
        wm(a, f));
    return f.Ve;
  }
  function Am(a, b) {
    if (null != b.Kc) {
      var c, d, e;
      c = b.Kc;
      b.Kc = c.e;
      e = null;
      for (d = c.list; null != d; d = d.e)
        null != d.code && (e = tl(e, xm(a, d.code)));
      if (c.code.Ta == Yj) {
        var f, g, h, k;
        g = $(a, c.code.a.a.x);
        h = $(a, c.code.a.a.y);
        k = null == c.code.a.a.z ? 1 : $(a, c.code.a.a.z);
        e = lm(a, g, h, k);
        d = tl(null, ml(0));
        for (f = 1; f <= e && b.Wf; f++)
          (d.Y.Q = mm(a, g, h, k, f)), um(a, c, d, b, Am);
      } else
        for (f = Bm(a, c.code).head; null != f && b.Wf; f = f.e) {
          g = f.w;
          h = e;
          k = !1;
          for (d = c.list; null != d; d = d.e) {
            if (null != d.code) {
              if (0 != gm(g.Y, h.Y)) {
                k = !0;
                break;
              }
              h = h.e;
            }
            g = g.e;
          }
          k || um(a, c, f.w, b, Am);
        }
      b.Kc = c;
    } else if (null == b.domain.code || ym(a, b.domain.code))
      b.Wf = !b.Zd(a, b.info);
  }
  function Cm(a, b, c, d) {
    var e = {};
    null == b
      ? d(a, c)
      : ((e.domain = b),
        (e.Kc = b.list),
        (e.Wf = 1),
        (e.info = c),
        (e.Zd = d),
        Am(a, e));
  }
  function Dm(a, b, c) {
    U(a, b + zl("[", c) + " out of domain");
  }
  function Em(a) {
    var b = null;
    if (null != a)
      for (a = a.list; null != a; a = a.e)
        for (var c = a.list; null != c; c = c.e)
          null == c.code && (b = tl(b, sl(c.value)));
    return b;
  }
  function Fm(a, b, c, d) {
    for (var e = b.Bf, f = 1; null != e; e = e.e, f++)
      for (var g = d.head; null != g; g = g.e)
        if (!vm(a, e.code, g.w)) {
          var h = zl("(", g.w);
          U(
            a,
            b.name +
              zl("[", c) +
              " contains " +
              h +
              " which not within specified set; see (" +
              f +
              ")"
          );
        }
  }
  function Gm(a, b, c) {
    function d() {
      Fm(a, b, c, e);
      f = Al(b.O, Il(c));
      f.value.set = e;
    }
    var e,
      f = yl(a, b.O, c);
    null != f
      ? (e = f.value.set)
      : null != b.assign
      ? ((e = Bm(a, b.assign)), d())
      : null != b.xa
      ? ((e = Bm(a, b.xa)), d())
      : U(a, "no value for " + b.name + zl("[", c));
    return e;
  }
  function Hm(a, b) {
    null != b.ga
      ? Fm(a, b.set, b.ga.w, b.ga.value.set)
      : (b.oe = Gm(a, b.set, b.w));
  }
  function Im(a, b) {
    var c = b.Yc,
      d,
      e,
      f,
      g = Array(20);
    x("Generating " + b.name + "...");
    d = c.set;
    Cm(a, d.domain, d, Jm);
    for (d = c.set.O.head.value.set.head; null != d; d = d.e) {
      f = Il(d.w);
      for (e = 0; e < c.set.X; e++) g[e] = null;
      for (e = 0; null != f; f = f.e) g[c.Z[e++] - 1] = f;
      for (e = 0; e < c.set.X; e++) g[e].e = g[e + 1];
      0 == b.q ? (f = null) : ((f = g[0]), (g[b.q - 1].e = null));
      e = yl(a, b.O, f);
      null == e && ((e = Al(b.O, f)), (e.value.set = Bl(a, qh, b.X)));
      f = g[b.q];
      g[c.set.X - 1].e = null;
      jm(e.value.set, f);
    }
    b.data = 1;
  }
  function Km(a, b, c) {
    var d = {};
    d.set = b;
    d.w = c;
    null != b.Yc && 0 == b.data && Im(a, b);
    if (1 == b.data)
      for (
        c = b.O.Xa, b.data = 2, d.ga = b.O.head;
        null != d.ga &&
        (zm(a, b.domain, d.ga.w, d, Hm) && Dm(a, b.name, d.ga.w), d.ga != c);
        d.ga = d.ga.e
      );
    d.ga = null;
    zm(a, d.set.domain, d.w, d, Hm) && Dm(a, b.name, d.w);
    return d.oe;
  }
  function Jm(a, b) {
    var c = Em(b.domain);
    Km(a, b, c);
    return 0;
  }
  function Lm(a, b, c, d) {
    var e, f;
    switch (b.type) {
      case mh:
        d != Math.floor(d) &&
          U(a, b.name + zl("[", c) + " = " + d + " not integer");
        break;
      case ah:
        0 != d &&
          1 != d &&
          U(a, b.name + zl("[", c) + " = " + d + " not binary");
    }
    e = b.wd;
    for (f = 1; null != e; e = e.e, f++) {
      var g;
      g = $(a, e.code);
      var h = function (e) {
        U(
          a,
          b.name +
            zl("[", c) +
            " = " +
            d +
            " not " +
            e +
            " " +
            g +
            "; see (" +
            f +
            ")"
        );
      };
      switch (e.jd) {
        case Dj:
          d < g || h("<");
          break;
        case Ej:
          d <= g || h("<=");
          break;
        case Fj:
          d != g && h("=");
          break;
        case Gj:
          d >= g || h(">=");
          break;
        case Hj:
          d > g || h(">");
          break;
        case Ij:
          d == g && h("<>");
      }
    }
    f = 1;
    for (e = b.qa; null != e; e = e.e, f++)
      (h = tl(null, ml(d))),
        vm(a, e.code, h) ||
          U(
            a,
            b.name +
              zl("[", c) +
              " = " +
              d +
              " not in specified set; see (" +
              f +
              ")"
          );
  }
  function Mm(a, b, c) {
    function d(d) {
      Lm(a, b, c, d);
      e = Al(b.O, Il(c));
      return (e.value.Q = d);
    }
    var e = yl(a, b.O, c);
    return null != e
      ? e.value.Q
      : null != b.assign
      ? d($(a, b.assign))
      : null != b.xa
      ? d($(a, b.xa))
      : null != b.Wc
      ? (null != b.Wc.M &&
          U(a, "cannot convert " + rl(b.Wc) + " to floating-point number"),
        d(b.Wc.Q))
      : U(a, "no value for " + b.name + zl("[", c));
  }
  function Nm(a, b) {
    null != b.ga
      ? Lm(a, b.S, b.ga.w, b.ga.value.Q)
      : (b.value = Mm(a, b.S, b.w));
  }
  function Om(a, b, c) {
    var d = {};
    d.S = b;
    d.w = c;
    if (1 == b.data)
      for (
        c = b.O.Xa, b.data = 2, d.ga = b.O.head;
        null != d.ga &&
        (zm(a, b.domain, d.ga.w, d, Nm) && Dm(a, b.name, d.ga.w), d.ga != c);
        d.ga = d.ga.e
      );
    d.ga = null;
    zm(a, d.S.domain, d.w, d, Nm) && Dm(a, b.name, d.w);
    return d.value;
  }
  function Pm(a, b, c, d) {
    var e,
      f = 1;
    for (e = b.wd; null != e; e = e.e, f++) {
      var g;
      g = xm(a, e.code);
      switch (e.jd) {
        case Dj:
          0 > gm(d, g) ||
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " = " + rl(d) + " not < " + g));
          break;
        case Ej:
          0 >= gm(d, g) ||
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " = " + rl(d) + " not <= " + g));
          break;
        case Fj:
          0 != gm(d, g) &&
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " = " + rl(d) + " not = " + g));
          break;
        case Gj:
          0 <= gm(d, g) ||
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " = " + rl(d) + " not >= " + g));
          break;
        case Hj:
          0 < gm(d, g) ||
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " = " + rl(d) + " not > " + g));
          break;
        case Ij:
          0 == gm(d, g) &&
            ((g = rl(g)),
            U(a, b.name + zl("[", c) + " <> " + rl(d) + " not > " + g));
      }
    }
    f = 1;
    for (e = b.qa; null != e; e = e.e, f++)
      (g = tl(null, sl(d))),
        vm(a, e.code, g) ||
          U(
            a,
            b.name,
            zl("[", c) +
              " = " +
              rl(d) +
              " not in specified set; see (" +
              f +
              ")"
          );
  }
  function Qm(a, b, c) {
    function d(d) {
      Pm(a, b, c, d);
      e = Al(b.O, Il(c));
      e.value.Y = sl(d);
      return d;
    }
    var e = yl(a, b.O, c);
    return null != e
      ? sl(e.value.Y)
      : null != b.assign
      ? d(xm(a, b.assign))
      : null != b.xa
      ? d(xm(a, b.xa))
      : null != b.Wc
      ? sl(b.Wc)
      : U(a, "no value for " + b.name + zl("[", c));
  }
  function Rm(a, b) {
    null != b.ga
      ? Pm(a, b.S, b.ga.w, b.ga.value.Y)
      : (b.value = Qm(a, b.S, b.w));
  }
  function Sm(a, b, c) {
    var d = {};
    d.S = b;
    d.w = c;
    if (1 == b.data)
      for (
        c = b.O.Xa, b.data = 2, d.ga = b.O.head;
        null != d.ga &&
        (zm(a, b.domain, d.ga.w, d, Rm) && Dm(a, b.name, d.ga.w), d.ga != c);
        d.ga = d.ga.e
      );
    d.ga = null;
    zm(a, d.S.domain, d.w, d, Rm) && Dm(a, b.name, d.w);
    return d.value;
  }
  function Tm(a, b) {
    var c = Em(b.domain);
    switch (b.type) {
      case N:
      case mh:
      case ah:
        Om(a, b, c);
        break;
      case T:
        Sm(a, b, c);
    }
    return 0;
  }
  function Um(a, b) {
    var c = b.t,
      d = b.w,
      e = yl(a, c.O, d);
    null != e
      ? (d = e.value.t)
      : ((e = Al(c.O, Il(d))),
        (d = e.value.t = {}),
        (d.C = 0),
        (d.t = c),
        (d.ga = e),
        (d.P = null == c.P ? 0 : $(a, c.P)),
        (d.W = null == c.W ? 0 : c.W == c.P ? d.P : $(a, c.W)),
        (d.ja = 0),
        (d.m = 0),
        (d.r = d.J = 0));
    b.oe = d;
  }
  function Vm(a, b, c) {
    var d = {};
    d.t = b;
    d.w = c;
    zm(a, d.t.domain, d.w, d, Um) && Dm(a, b.name, d.w);
    return d.oe;
  }
  function Wm(a, b, c) {
    var d = null,
      e = yl(a, b.O, c);
    if (null != e) c = e.value.H;
    else {
      e = Al(b.O, Il(c));
      c = e.value.H = {};
      c.ea = 0;
      c.H = b;
      c.ga = e;
      c.form = Xm(a, b.code);
      if (null == b.P && null == b.W)
        (c.form = qm(a, c.form, function (a) {
          d = a;
        })),
          (c.P = c.W = -d);
      else if (null != b.P && null == b.W)
        (c.form = pm(a, 1, c.form, -1, Xm(a, b.P))),
          (c.form = qm(a, c.form, function (a) {
            d = a;
          })),
          (c.P = -d),
          (c.W = 0);
      else if (null == b.P && null != b.W)
        (c.form = pm(a, 1, c.form, -1, Xm(a, b.W))),
          (c.form = qm(a, c.form, function (a) {
            d = a;
          })),
          (c.P = 0),
          (c.W = -d);
      else if (b.P == b.W)
        (c.form = pm(a, 1, c.form, -1, Xm(a, b.P))),
          (c.form = qm(a, c.form, function (a) {
            d = a;
          })),
          (c.P = c.W = -d);
      else {
        var f = null,
          g = null;
        c.form = qm(a, c.form, function (a) {
          d = a;
        });
        qm(a, Xm(a, b.P), function (a) {
          f = a;
        });
        qm(a, Xm(a, b.W), function (a) {
          g = a;
        });
        c.P = Ml(a, f, d);
        c.W = Ml(a, g, d);
      }
      c.m = 0;
      c.r = c.J = 0;
    }
    return c;
  }
  function Ym(a, b) {
    b.oe = Wm(a, b.H, b.w);
  }
  function Zm(a, b, c) {
    var d = {};
    d.H = b;
    d.w = c;
    zm(a, d.H.domain, d.w, d, Ym) && Dm(a, b.name, d.w);
    return d.oe;
  }
  function $m(a, b) {
    var c = Em(b.domain);
    Zm(a, b, c);
    return 0;
  }
  function an(a, b) {
    var c = $(a, b.code.a.loop.x);
    switch (b.code.Ta) {
      case dk:
        b.value = Ll(a, b.value, c);
        break;
      case ek:
        b.value = Ol(a, b.value, c);
        break;
      case fk:
        b.value > c && (b.value = c);
        break;
      case gk:
        b.value < c && (b.value = c);
    }
    return 0;
  }
  function $(a, b) {
    var c, d, e;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return b.value.Q;
    switch (b.Ta) {
      case Fi:
        c = b.a.Q;
        break;
      case Ii:
        d = null;
        for (e = b.a.S.list; null != e; e = e.e) d = tl(d, xm(a, e.x));
        c = Om(a, b.a.S.S, d);
        break;
      case Li:
        d = null;
        for (e = b.a.t.list; null != e; e = e.e) d = tl(d, xm(a, e.x));
        e = Vm(a, b.a.t.t, d);
        switch (b.a.t.Ac) {
          case Ai:
            c = null == e.t.P ? -s : e.P;
            break;
          case Bi:
            c = null == e.t.W ? +s : e.W;
            break;
          case Ci:
            c = e.m;
            break;
          case Di:
            c = e.r;
            break;
          case Ei:
            c = e.J;
        }
        break;
      case Mi:
        d = null;
        for (e = b.a.H.list; null != e; e = e.e) d = tl(d, xm(a, e.x));
        e = Zm(a, b.a.H.H, d);
        switch (b.a.H.Ac) {
          case Ai:
            c = null == e.H.P ? -s : e.P;
            break;
          case Bi:
            c = null == e.H.W ? +s : e.W;
            break;
          case Ci:
            c = e.m;
            break;
          case Di:
            c = e.r;
            break;
          case Ei:
            c = e.J;
        }
        break;
      case Qi:
        c = bn(a.Gd);
        break;
      case Ri:
        c = dm(a.Gd) / 2147483648;
        break;
      case Si:
        c = em(a);
        break;
      case Ti:
        c = Math.round(Date.now() / 1e3);
        break;
      case Ui:
        e = xm(a, b.a.a.x);
        null == e.M
          ? (c = e.Q)
          : vg(e.M, function (a) {
              c = a;
            }) && U(a, "cannot convert " + rl(e) + " to floating-point number");
        break;
      case Zi:
        c = +$(a, b.a.a.x);
        break;
      case $i:
        c = -$(a, b.a.a.x);
        break;
      case bj:
        c = Math.abs($(a, b.a.a.x));
        break;
      case cj:
        c = Math.ceil($(a, b.a.a.x));
        break;
      case dj:
        c = Math.floor($(a, b.a.a.x));
        break;
      case ej:
        c = Tl(a, $(a, b.a.a.x));
        break;
      case fj:
        c = Ul(a, $(a, b.a.a.x));
        break;
      case gj:
        c = Vl(a, $(a, b.a.a.x));
        break;
      case hj:
        c = Wl(a, $(a, b.a.a.x));
        break;
      case ij:
        c = Xl(a, $(a, b.a.a.x));
        break;
      case jj:
        c = Yl(a, $(a, b.a.a.x));
        break;
      case kj:
        c = Zl($(a, b.a.a.x));
        break;
      case xj:
        c = $l($(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case lj:
        c = am(a, $(a, b.a.a.x), 0);
        break;
      case yj:
        c = am(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case mj:
        c = bm(a, $(a, b.a.a.x), 0);
        break;
      case zj:
        c = bm(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case pj:
        c = Ll(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case qj:
        c = Ml(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case rj:
        c = Nl(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case sj:
        c = Ol(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case tj:
        c = Pl(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case uj:
        c = Ql(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case vj:
        c = Rl($(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case wj:
        c = Sl(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case Aj:
        c = cm(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case Bj:
        c = fm(a, $(a, b.a.a.x), $(a, b.a.a.y));
        break;
      case nj:
        c = Bm(a, b.a.a.x).size;
        break;
      case oj:
        e = xm(a, b.a.a.x);
        d = null == e.M ? String(e.Q) : e.M;
        c = d.length;
        break;
      case Vj:
        e = xm(a, b.a.a.x);
        d = null == e.M ? String(e.Q) : e.M;
        e = xm(a, b.a.a.y);
        c = cn(a, d, null == e.M ? String(e.Q) : e.M);
        break;
      case Zj:
        c = ym(a, b.a.a.x)
          ? $(a, b.a.a.y)
          : null == b.a.a.z
          ? 0
          : $(a, b.a.a.z);
        break;
      case bk:
        c = +s;
        for (e = b.a.list; null != e; e = e.e)
          (d = $(a, e.x)), c > d && (c = d);
        break;
      case ck:
        c = -s;
        for (e = b.a.list; null != e; e = e.e)
          (d = $(a, e.x)), c < d && (c = d);
        break;
      case dk:
        e = {};
        e.code = b;
        e.value = 0;
        Cm(a, b.a.loop.domain, e, an);
        c = e.value;
        break;
      case ek:
        e = {};
        e.code = b;
        e.value = 1;
        Cm(a, b.a.loop.domain, e, an);
        c = e.value;
        break;
      case fk:
        e = {};
        e.code = b;
        e.value = +s;
        Cm(a, b.a.loop.domain, e, an);
        e.value == +s && U(a, "min{} over empty set; result undefined");
        c = e.value;
        break;
      case gk:
        (e = {}),
          (e.code = b),
          (e.value = -s),
          Cm(a, b.a.loop.domain, e, an),
          e.value == -s && U(a, "max{} over empty set; result undefined"),
          (c = e.value);
    }
    b.valid = 1;
    return (b.value.Q = c);
  }
  function xm(a, b) {
    var c;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return sl(b.value.Y);
    switch (b.Ta) {
      case Gi:
        c = nl(b.a.M);
        break;
      case Hi:
        c = sl(b.a.index.ya.value);
        break;
      case Ji:
        var d;
        d = null;
        for (c = b.a.S.list; null != c; c = c.e) d = tl(d, xm(a, c.x));
        c = Sm(a, b.a.S.S, d);
        break;
      case Vi:
        c = ml($(a, b.a.a.x));
        break;
      case Cj:
        d = xm(a, b.a.a.x);
        c = xm(a, b.a.a.y);
        c = nl(
          (null == d.M ? String(d.Q) : d.M) + (null == c.M ? String(c.Q) : c.M)
        );
        break;
      case Zj:
        c = ym(a, b.a.a.x)
          ? xm(a, b.a.a.y)
          : null == b.a.a.z
          ? ml(0)
          : xm(a, b.a.a.z);
        break;
      case Uj:
      case ak:
        var e;
        c = xm(a, b.a.a.x);
        c = null == c.M ? String(c.Q) : c.M;
        b.Ta == Uj
          ? ((e = $(a, b.a.a.y)),
            e != Math.floor(e) &&
              U(a, "substr('...', " + e + "); non-integer second argument"),
            (1 > e || e > c.length + 1) &&
              U(a, "substr('...', " + e + "); substring out of range"))
          : ((e = $(a, b.a.a.y)),
            (d = $(a, b.a.a.z)),
            (e == Math.floor(e) && d == Math.floor(d)) ||
              U(
                a,
                "substr('...', " +
                  e +
                  ", " +
                  d +
                  "); non-integer second and/or third argument"
              ),
            (1 > e || 0 > d || e + d > c.length + 1) &&
              U(
                a,
                "substr('...', " + e + ", " + d + "); substring out of range"
              ));
        c = nl(c.slice(e - 1, e + d - 1));
        break;
      case Xj:
        (d = $(a, b.a.a.x)),
          (c = xm(a, b.a.a.y)),
          (c = dn(a, d, null == c.M ? String(c.Q) : c.M)),
          (c = nl(c));
    }
    b.valid = 1;
    b.value.Y = sl(c);
    return c;
  }
  function en(a, b) {
    var c = 0;
    switch (b.code.Ta) {
      case hk:
        b.value &= ym(a, b.code.a.loop.x);
        b.value || (c = 1);
        break;
      case ik:
        (b.value |= ym(a, b.code.a.loop.x)), b.value && (c = 1);
    }
    return c;
  }
  function ym(a, b) {
    var c, d;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return b.value.sg;
    switch (b.Ta) {
      case Wi:
        c = 0 != $(a, b.a.a.x);
        break;
      case aj:
        c = !ym(a, b.a.a.x);
        break;
      case Dj:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) < $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 > gm(c, d)));
        break;
      case Ej:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) <= $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 >= gm(c, d)));
        break;
      case Fj:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) == $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 == gm(c, d)));
        break;
      case Gj:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) >= $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 <= gm(c, d)));
        break;
      case Hj:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) > $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 < gm(c, d)));
        break;
      case Ij:
        b.a.a.x.type == N
          ? (c = $(a, b.a.a.x) != $(a, b.a.a.y))
          : ((c = xm(a, b.a.a.x)), (d = xm(a, b.a.a.y)), (c = 0 != gm(c, d)));
        break;
      case Jj:
        c = ym(a, b.a.a.x) && ym(a, b.a.a.y);
        break;
      case Kj:
        c = ym(a, b.a.a.x) || ym(a, b.a.a.y);
        break;
      case Qj:
        c = fn(a, b.a.a.x);
        c = vm(a, b.a.a.y, c);
        break;
      case Rj:
        c = fn(a, b.a.a.x);
        c = !vm(a, b.a.a.y, c);
        break;
      case Sj:
        d = Bm(a, b.a.a.x);
        c = 1;
        for (d = d.head; null != d; d = d.e)
          if (!vm(a, b.a.a.y, d.w)) {
            c = 0;
            break;
          }
        break;
      case Tj:
        d = Bm(a, b.a.a.x);
        c = 1;
        for (d = d.head; null != d; d = d.e)
          if (vm(a, b.a.a.y, d.w)) {
            c = 0;
            break;
          }
        break;
      case hk:
        c = {};
        c.code = b;
        c.value = 1;
        Cm(a, b.a.loop.domain, c, en);
        c = c.value;
        break;
      case ik:
        (c = {}),
          (c.code = b),
          (c.value = 0),
          Cm(a, b.a.loop.domain, c, en),
          (c = c.value);
    }
    b.valid = 1;
    return (b.value.sg = c);
  }
  function fn(a, b) {
    var c;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return Il(b.value.w);
    switch (b.Ta) {
      case Ni:
        c = null;
        for (var d = b.a.list; null != d; d = d.e) c = tl(c, xm(a, d.x));
        break;
      case Xi:
        c = tl(null, xm(a, b.a.a.x));
    }
    b.valid = 1;
    b.value.w = Il(c);
    return c;
  }
  function gn(a, b) {
    var c;
    switch (b.code.Ta) {
      case jk:
        c = fn(a, b.code.a.loop.x);
        null == yl(a, b.value, c) && jm(b.value, c);
        break;
      case kk:
        jm(b.value, Em(b.code.a.loop.domain));
    }
    return 0;
  }
  function Bm(a, b) {
    var c, d;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return km(a, b.value.set);
    switch (b.Ta) {
      case Ki:
        c = null;
        for (d = b.a.set.list; null != d; d = d.e) c = tl(c, xm(a, d.x));
        c = km(a, Km(a, b.a.set.set, c));
        break;
      case Oi:
        c = Bl(a, qh, b.q);
        for (d = b.a.list; null != d; d = d.e) ul(a, c, fn(a, d.x));
        break;
      case Lj:
        d = Bm(a, b.a.a.x);
        for (c = Bm(a, b.a.a.y).head; null != c; c = c.e)
          null == yl(a, d, c.w) && jm(d, Il(c.w));
        c = d;
        break;
      case Mj:
        var e = Bm(a, b.a.a.x);
        d = Bm(a, b.a.a.y);
        c = Bl(a, qh, e.q);
        for (e = e.head; null != e; e = e.e)
          null == yl(a, d, e.w) && jm(c, Il(e.w));
        break;
      case Nj:
        d = Bm(a, b.a.a.x);
        c = Bm(a, b.a.a.y);
        for (var f = Bl(a, qh, d.q), e = d.head; null != e; e = e.e)
          null == yl(a, c, e.w) && jm(f, Il(e.w));
        for (e = c.head; null != e; e = e.e)
          null == yl(a, d, e.w) && jm(f, Il(e.w));
        c = f;
        break;
      case Oj:
        e = Bm(a, b.a.a.x);
        d = Bm(a, b.a.a.y);
        c = Bl(a, qh, e.q);
        for (e = e.head; null != e; e = e.e)
          null != yl(a, d, e.w) && jm(c, Il(e.w));
        break;
      case Pj:
        e = Bm(a, b.a.a.x);
        d = Bm(a, b.a.a.y);
        var g, h;
        c = Bl(a, qh, e.q + d.q);
        for (e = e.head; null != e; e = e.e)
          for (f = d.head; null != f; f = f.e) {
            g = Il(e.w);
            for (h = f.w; null != h; h = h.e) g = tl(g, sl(h.Y));
            jm(c, g);
          }
        break;
      case Yj:
        d = $(a, b.a.a.x);
        c = $(a, b.a.a.y);
        e = null == b.a.a.z ? 1 : $(a, b.a.a.z);
        f = Bl(a, qh, 1);
        g = lm(a, d, c, e);
        for (h = 1; h <= g; h++) jm(f, tl(null, ml(mm(a, d, c, e, h))));
        c = f;
        break;
      case Zj:
        c = ym(a, b.a.a.x) ? Bm(a, b.a.a.y) : Bm(a, b.a.a.z);
        break;
      case jk:
        d = {};
        d.code = b;
        d.value = Bl(a, qh, b.q);
        Cm(a, b.a.loop.domain, d, gn);
        c = d.value;
        break;
      case kk:
        (d = {}),
          (d.code = b),
          (d.value = Bl(a, qh, b.q)),
          Cm(a, b.a.loop.domain, d, gn),
          (c = d.value);
    }
    b.valid = 1;
    b.value.set = km(a, c);
    return c;
  }
  function hn() {}
  function vm(a, b, c) {
    var d, e, f;
    switch (b.Ta) {
      case Ki:
        f = null;
        for (e = b.a.set.list; null != e; e = e.e) f = tl(f, xm(a, e.x));
        b = Km(a, b.a.set.set, f);
        f = im(c, b.q);
        d = null != yl(a, b, f);
        break;
      case Oi:
        d = 0;
        f = im(c, b.q);
        for (
          e = b.a.list;
          null != e && !((c = fn(a, e.x)), (d = 0 == hm(f, c)));
          e = e.e
        );
        break;
      case Lj:
        d = vm(a, b.a.a.x, c) || vm(a, b.a.a.y, c);
        break;
      case Mj:
        d = vm(a, b.a.a.x, c) && !vm(a, b.a.a.y, c);
        break;
      case Nj:
        f = vm(a, b.a.a.x, c);
        a = vm(a, b.a.a.y, c);
        d = (f && !a) || (!f && a);
        break;
      case Oj:
        d = vm(a, b.a.a.x, c) && vm(a, b.a.a.y, c);
        break;
      case Pj:
        if ((d = vm(a, b.a.a.x, c))) {
          for (f = 1; f <= b.a.a.x.q; f++) c = c.e;
          d = vm(a, b.a.a.y, c);
        }
        break;
      case Yj:
        f = $(a, b.a.a.x);
        e = $(a, b.a.a.y);
        b = null == b.a.a.z ? 1 : $(a, b.a.a.z);
        lm(a, f, e, b);
        if (null != c.Y.M) {
          d = 0;
          break;
        }
        c = c.Y.Q;
        if ((0 < b && !(f <= c && c <= e)) || (0 > b && !(e <= c && c <= f))) {
          d = 0;
          break;
        }
        d = mm(a, f, e, b, (((c - f) / b + 0.5) | 0) + 1) == c;
        break;
      case Zj:
        d = ym(a, b.a.a.x) ? vm(a, b.a.a.y, c) : vm(a, b.a.a.z, c);
        break;
      case jk:
        U(a, "implementation restriction; in/within setof{} not allowed");
        break;
      case kk:
        (f = im(c, b.q)), (d = 0 == zm(a, b.a.loop.domain, f, null, hn));
    }
    return d;
  }
  function jn(a, b) {
    switch (b.code.Ta) {
      case dk:
        var c;
        c = Xm(a, b.code.a.loop.x);
        for (null == b.value ? (b.value = c) : (b.Xa.e = c); null != c; c = c.e)
          b.Xa = c;
    }
    return 0;
  }
  function Xm(a, b) {
    var c;
    b.T && b.valid && ((b.valid = 0), rm(b.type, b.value));
    if (b.valid) return om(b.value.form);
    switch (b.Ta) {
      case Li:
        var d = null;
        for (c = b.a.t.list; null != c; c = c.e) d = tl(d, xm(a, c.x));
        c = Vm(a, b.a.t.t, d);
        d = { u: 1 };
        d.t = c;
        d.e = null;
        c = d;
        break;
      case Yi:
        c = nm($(a, b.a.a.x));
        break;
      case Zi:
        c = pm(a, 0, nm(0), 1, Xm(a, b.a.a.x));
        break;
      case $i:
        c = pm(a, 0, nm(0), -1, Xm(a, b.a.a.x));
        break;
      case pj:
        c = pm(a, 1, Xm(a, b.a.a.x), 1, Xm(a, b.a.a.y));
        break;
      case qj:
        c = pm(a, 1, Xm(a, b.a.a.x), -1, Xm(a, b.a.a.y));
        break;
      case sj:
        c =
          b.a.a.x.type == N
            ? pm(a, $(a, b.a.a.x), Xm(a, b.a.a.y), 0, nm(0))
            : pm(a, $(a, b.a.a.y), Xm(a, b.a.a.x), 0, nm(0));
        break;
      case tj:
        c = pm(a, Pl(a, 1, $(a, b.a.a.y)), Xm(a, b.a.a.x), 0, nm(0));
        break;
      case Zj:
        c = ym(a, b.a.a.x)
          ? Xm(a, b.a.a.y)
          : null == b.a.a.z
          ? nm(0)
          : Xm(a, b.a.a.z);
        break;
      case dk:
        c = {};
        c.code = b;
        c.value = nm(0);
        c.Xa = null;
        Cm(a, b.a.loop.domain, c, jn);
        c = c.value;
        for (var e, f = 0, d = c; null != d; d = d.e)
          null == d.t ? (f = Ll(a, f, d.u)) : (d.t.ja = Ll(a, d.t.ja, d.u));
        e = c;
        c = null;
        for (d = e; null != d; d = e)
          (e = d.e),
            null == d.t && 0 != f
              ? ((d.u = f), (f = 0), (d.e = c), (c = d))
              : null != d.t &&
                0 != d.t.ja &&
                ((d.u = d.t.ja), (d.t.ja = 0), (d.e = c), (c = d));
    }
    b.valid = 1;
    b.value.form = om(c);
    return c;
  }
  var kn = (exports.mpl_tab_num_args = function (a) {
      return a.ff;
    }),
    ln = (exports.mpl_tab_get_arg = function (a, b) {
      return a.a[b];
    });
  exports.mpl_tab_get_args = function (a) {
    return a.a;
  };
  var mn = (exports.mpl_tab_num_flds = function (a) {
      return a.Wa;
    }),
    nn = (exports.mpl_tab_get_name = function (a, b) {
      return a.name[b];
    }),
    on = (exports.mpl_tab_get_type = function (a, b) {
      return a.type[b];
    }),
    pn = (exports.mpl_tab_get_num = function (a, b) {
      return a.Q[b];
    }),
    qn = (exports.mpl_tab_get_str = function (a, b) {
      return a.M[b];
    }),
    rn = (exports.mpl_tab_set_num = function (a, b, c) {
      a.type[b] = "N";
      a.Q[b] = c;
    }),
    sn = (exports.mpl_tab_set_str = function (a, b, c) {
      a.type[b] = "S";
      a.M[b] = c;
    });
  function tn(a, b) {
    var c = a.Mc,
      d,
      e,
      f;
    f = 0;
    for (d = b.v.Oc.list; null != d; d = d.e)
      switch ((f++, d.code.type)) {
        case N:
          c.type[f] = "N";
          c.Q[f] = $(a, d.code);
          c.M[f][0] = "\x00";
          break;
        case T:
          (e = xm(a, d.code)),
            null == e.M
              ? ((c.type[f] = "N"), (c.Q[f] = e.Q), (c.M[f][0] = "\x00"))
              : ((c.type[f] = "S"), (c.Q[f] = 0), (c.M[f] = e.M));
      }
    c = a.Mc;
    c.link.writeRecord(c) &&
      U(a, "error on writing data to table " + a.ib.v.nd.name);
    return 0;
  }
  function un(a, b) {
    ym(a, b.code) || U(a, "check" + zl("[", Em(b.domain)) + " failed");
    return 0;
  }
  function vn(a, b, c) {
    var d = c.value.set;
    wn(a, b.name + zl("[", c.w) + (null == d.head ? " is empty" : ":"));
    for (b = d.head; null != b; b = b.e) wn(a, "   " + zl("(", b.w));
  }
  function xn(a, b, c) {
    switch (b.type) {
      case N:
      case mh:
      case ah:
        wn(a, b.name + zl("[", c.w) + " = " + c.value.Q);
        break;
      case T:
        wn(a, b.name + zl("[", c.w) + " = " + rl(c.value.Y));
    }
  }
  function yn(a, b, c, d) {
    d == zi || d == Di
      ? wn(a, b.name + zl("[", c.w) + ".val = " + c.value.t.r)
      : d == Ai
      ? wn(
          a,
          b.name +
            zl("[", c.w) +
            ".lb = " +
            (null == c.value.t.t.P ? -s : c.value.t.P)
        )
      : d == Bi
      ? wn(
          a,
          b.name +
            zl("[", c.w) +
            ".ub = " +
            (null == c.value.t.t.W ? +s : c.value.t.W)
        )
      : d == Ci
      ? wn(a, b.name + zl("[", c.w) + ".status = " + c.value.t.m)
      : d == Ei && wn(a, b.name + zl("[", c.w) + ".dual = " + c.value.t.J);
  }
  function zn(a, b, c, d) {
    d == zi || d == Di
      ? wn(a, b.name + zl("[", c.w) + ".val = " + c.value.H.r)
      : d == Ai
      ? wn(
          a,
          b.name +
            zl("[", c.w) +
            ".lb = " +
            (null == c.value.H.H.P ? -s : c.value.H.P)
        )
      : d == Bi
      ? wn(
          a,
          b.name +
            zl("[", c.w) +
            ".ub = " +
            (null == c.value.H.H.W ? +s : c.value.H.W)
        )
      : d == Ci
      ? wn(a, b.name + zl("[", c.w) + ".status = " + c.value.H.m)
      : d == Ei && wn(a, b.name + zl("[", c.w) + ".dual = " + c.value.H.J);
  }
  function An(a, b) {
    for (var c, d = b.list; null != d; d = d.e)
      if (d.type == kh) (c = d.v.ya), wn(a, c.name + " = " + rl(c.value));
      else if (d.type == uh) {
        var e = d.v.set;
        null != e.assign
          ? Cm(a, e.domain, e, Jm)
          : (null != e.Yc && 0 == e.data && Im(a, e),
            null != e.O.head && Km(a, e, e.O.head.w));
        null == e.O.head && wn(a, e.name + " has empty content");
        for (c = e.O.head; null != c; c = c.e) vn(a, e, c);
      } else if (d.type == sh)
        for (
          e = d.v.S,
            null != e.assign
              ? Cm(a, e.domain, e, Tm)
              : null != e.O.head &&
                (e.type != T ? Om(a, e, e.O.head.w) : Sm(a, e, e.O.head.w)),
            null == e.O.head && wn(a, e.name + " has empty content"),
            c = e.O.head;
          null != c;
          c = c.e
        )
          xn(a, e, c);
      else if (d.type == yh)
        for (
          e = d.v.t,
            null == e.O.head && wn(a, e.name + " has empty content"),
            c = e.O.head;
          null != c;
          c = c.e
        )
          yn(a, e, c, zi);
      else if (d.type == ch)
        for (
          e = d.v.H,
            null == e.O.head && wn(a, e.name + " has empty content"),
            c = e.O.head;
          null != c;
          c = c.e
        )
          zn(a, e, c, zi);
      else if (d.type == hh)
        if (
          ((e = d.v.code),
          e.Ta == Ii || e.Ta == Ji || e.Ta == Ki || e.Ta == Li || e.Ta == Mi)
        ) {
          c = a;
          var f = { value: {} },
            g = void 0;
          f.w = null;
          for (g = e.a.S.list || e.a.t.list; null != g; g = g.e)
            f.w = tl(f.w, xm(c, g.x));
          switch (e.Ta) {
            case Ii:
              f.value.Q = Om(c, e.a.S.S, f.w);
              xn(c, e.a.S.S, f);
              break;
            case Ji:
              f.value.Y = Sm(c, e.a.S.S, f.w);
              xn(c, e.a.S.S, f);
              break;
            case Ki:
              f.value.set = Km(c, e.a.set.set, f.w);
              vn(c, e.a.set.set, f);
              break;
            case Li:
              f.value.t = Vm(c, e.a.t.t, f.w);
              yn(c, e.a.t.t, f, e.a.t.Ac);
              break;
            case Mi:
              (f.value.H = Zm(c, e.a.H.H, f.w)), zn(c, e.a.H.H, f, e.a.H.Ac);
          }
        } else
          switch (((c = a), e.type)) {
            case N:
              e = $(c, e);
              wn(c, String(e));
              break;
            case T:
              e = xm(c, e);
              wn(c, rl(e));
              break;
            case nh:
              e = ym(c, e);
              wn(c, e ? "true" : "false");
              break;
            case xh:
              e = fn(c, e);
              wn(c, zl("(", e));
              break;
            case fh:
              e = Bm(c, e);
              0 == e.head && wn(c, "set is empty");
              for (e = e.head; null != e; e = e.e) wn(c, "   " + zl("(", e.w));
              break;
            case jh:
              for (
                f = void 0,
                  e = Xm(c, e),
                  null == e && wn(c, "linear form is empty"),
                  f = e;
                null != f;
                f = f.e
              )
                null == f.t
                  ? wn(c, "   " + f.u)
                  : wn(c, "   " + f.u + " " + f.t.t.name + zl("[", f.t.ga.w));
          }
    return 0;
  }
  function Bn(a, b) {
    null == a.Hg
      ? "\n" == b
        ? (a.he(a.dd, a.le), (a.dd = ""))
        : (a.dd += b)
      : a.Hg(b);
  }
  function Cn(a, b) {
    for (var c = 0; c < b.length; c++) Bn(a, b[c]);
  }
  function Dn(a, b) {
    var c,
      d,
      e,
      f,
      g,
      h = xm(a, b.zd);
    d = null == h.M ? String(h.Q) : h.M;
    c = b.list;
    for (f = 0; f < d.length; f++)
      if ("%" == d[f])
        if (((e = f++), "%" == d[f])) Bn(a, "%");
        else {
          if (null == c) break;
          for (
            ;
            "-" == d[f] ||
            "+" == d[f] ||
            " " == d[f] ||
            "#" == d[f] ||
            "0" == d[f];

          )
            f++;
          for (; wa(d[f]); ) f++;
          if ("." == d[f]) for (f++; wa(d[f]); ) f++;
          if (
            "d" == d[f] ||
            "i" == d[f] ||
            "e" == d[f] ||
            "E" == d[f] ||
            "f" == d[f] ||
            "F" == d[f] ||
            "g" == d[f] ||
            "G" == d[f]
          ) {
            switch (c.code.type) {
              case N:
                g = $(a, c.code);
                break;
              case T:
                h = xm(a, c.code);
                null != h.M &&
                  U(a, "cannot convert " + rl(h) + " to floating-point number");
                g = h.Q;
                break;
              case nh:
                g = ym(a, c.code) ? 1 : 0;
            }
            "d" == d[f] || "i" == d[f]
              ? ((-2147483647 <= g && 2147483647 >= g) ||
                  U(a, "cannot convert " + g + " to integer"),
                Cn(a, xa(d.slice(e, f + 1), Math.floor(g + 0.5) | 0)))
              : Cn(a, xa(d.slice(e, f + 1), g));
          } else if ("s" == d[f]) {
            switch (c.code.type) {
              case N:
                g = String($(a, c.code));
                break;
              case nh:
                g = ym(a, c.code) ? "T" : "F";
                break;
              case T:
                (h = xm(a, c.code)), (g = null == h.M ? String(h.Q) : h.M);
            }
            Cn(a, xa(d.slice(e, f + 1), g));
          } else U(a, "format specifier missing or invalid");
          c = c.e;
        }
      else
        "\\" == d[f]
          ? (f++,
            "t" == d[f]
              ? Bn(a, "\t")
              : "n" == d[f]
              ? Bn(a, "\n")
              : "\x00" == d[f]
              ? U(
                  a,
                  "invalid use of escape character \\ in format control string"
                )
              : Bn(a, d[f]))
          : Bn(a, d[f]);
    return 0;
  }
  function En(a, b) {
    for (var c = a.ib, d = b.list; null != d; d = d.e) Fn(a, d);
    a.ib = c;
    return 0;
  }
  function Fn(a, b) {
    a.ib = b;
    switch (b.type) {
      case ch:
        x("Generating " + b.v.H.name + "...");
        var c = b.v.H;
        Cm(a, c.domain, c, $m);
        break;
      case wh:
        switch (b.v.nd.type) {
          case lh:
            x("Reading " + b.v.nd.name + "...");
            break;
          case rh:
            x("Writing " + b.v.nd.name + "...");
        }
        var c = b.v.nd,
          d,
          e,
          f,
          g;
        a.Mc = f = {};
        f.id = 0;
        f.link = null;
        f.ff = 0;
        f.a = null;
        f.Wa = 0;
        f.name = null;
        f.type = null;
        f.Q = null;
        f.M = null;
        for (d = c.a; null != d; d = d.e) f.ff++;
        f.a = Array(1 + f.ff);
        for (g = 1; g <= f.ff; g++) f.a[g] = null;
        g = 0;
        for (d = c.a; null != d; d = d.e)
          g++,
            (e = xm(a, d.code)),
            (e = null == e.M ? String(e.Q) : e.M),
            (f.a[g] = e);
        switch (c.type) {
          case lh:
            g = c.v.qa.set;
            null != g &&
              (g.data && U(a, g.name + " already provided with data"),
              (Al(g.O, null).value.set = Bl(a, qh, g.X)),
              (g.data = 1));
            for (e = c.v.qa.list; null != e; e = e.e)
              e.S.data && U(a, e.S.name + " already provided with data"),
                (e.S.data = 1);
            for (e = c.v.qa.We; null != e; e = e.e) f.Wa++;
            for (e = c.v.qa.list; null != e; e = e.e) f.Wa++;
            f.name = Array(1 + f.Wa);
            f.type = Array(1 + f.Wa);
            f.Q = new Float64Array(1 + f.Wa);
            f.M = Array(1 + f.Wa);
            g = 0;
            for (e = c.v.qa.We; null != e; e = e.e)
              g++,
                (f.name[g] = e.name),
                (f.type[g] = "?"),
                (f.Q[g] = 0),
                (f.M[g] = "");
            for (e = c.v.qa.list; null != e; e = e.e)
              g++,
                (f.name[g] = e.name),
                (f.type[g] = "?"),
                (f.Q[g] = 0),
                (f.M[g] = "");
            for (Gn(a, "R"); ; ) {
              for (g = 1; g <= f.Wa; g++) f.type[g] = "?";
              g = a;
              d = g.Mc;
              d = d.link.readRecord(d);
              0 < d &&
                U(g, "error on reading data from table " + g.ib.v.nd.name);
              if (d) break;
              for (g = 1; g <= f.Wa; g++)
                "?" == f.type[g] &&
                  U(a, "field " + f.name[g] + " missing in input table");
              d = null;
              g = 0;
              for (e = c.v.qa.We; null != e; e = e.e)
                switch ((g++, f.type[g])) {
                  case "N":
                    d = tl(d, ml(f.Q[g]));
                    break;
                  case "S":
                    d = tl(d, nl(f.M[g]));
                }
              null != c.v.qa.set && ul(a, c.v.qa.set.O.head.value.set, Il(d));
              for (e = c.v.qa.list; null != e; e = e.e) {
                var h;
                g++;
                null != yl(a, e.S.O, d) &&
                  U(a, e.S.name + zl("[", d) + " already defined");
                h = Al(e.S.O, Il(d));
                switch (e.S.type) {
                  case N:
                  case mh:
                  case ah:
                    "N" != f.type[g] &&
                      U(a, e.S.name + " requires numeric data");
                    h.value.Q = f.Q[g];
                    break;
                  case T:
                    switch (f.type[g]) {
                      case "N":
                        h.value.Y = ml(f.Q[g]);
                        break;
                      case "S":
                        h.value.Y = nl(f.M[g]);
                    }
                }
              }
            }
            a.Mc = null;
            break;
          case rh:
            for (d = c.v.Oc.list; null != d; d = d.e) f.Wa++;
            f.name = Array(1 + f.Wa);
            f.type = Array(1 + f.Wa);
            f.Q = new Float64Array(1 + f.Wa);
            f.M = Array(1 + f.Wa);
            g = 0;
            for (d = c.v.Oc.list; null != d; d = d.e)
              g++,
                (f.name[g] = d.name),
                (f.type[g] = "?"),
                (f.Q[g] = 0),
                (f.M[g] = "");
            Gn(a, "W");
            Cm(a, c.v.Oc.domain, c, tn);
            c = a.Mc;
            c.link.flush(c);
            a.Mc = null;
        }
        break;
      case bh:
        x("Checking (line " + b.bb + ")...");
        c = b.v.Rg;
        Cm(a, c.domain, c, un);
        break;
      case dh:
        wn(a, "Display statement at line " + b.bb);
        c = b.v.Sg;
        Cm(a, c.domain, c, An);
        break;
      case th:
        c = b.v.nh;
        null == c.Ea
          ? (a.le = null)
          : ((f = xm(a, c.Ea)), (a.le = null == f.M ? f.Q : f.M));
        Cm(a, c.domain, c, Dn);
        break;
      case ih:
        (c = b.v.Ug), Cm(a, c.domain, c, En);
    }
  }
  function Hn(a) {
    var b;
    for (b = a.uc; null != b; b = b.e)
      switch (b.type) {
        case uh:
          b.v.set.O = Bl(a, fh, b.v.set.q);
          break;
        case sh:
          switch (b.v.S.type) {
            case N:
            case mh:
            case ah:
              b.v.S.O = Bl(a, N, b.v.S.q);
              break;
            case T:
              b.v.S.O = Bl(a, T, b.v.S.q);
          }
          break;
        case yh:
          b.v.t.O = Bl(a, gh, b.v.t.q);
          break;
        case ch:
          b.v.H.O = Bl(a, eh, b.v.H.q);
      }
  }
  function In(a, b, c) {
    a.bb = 0;
    a.Vc = 0;
    a.l = "\n";
    a.b = 0;
    a.Bb = 0;
    a.h = "";
    a.value = 0;
    a.Hf = Ah;
    a.Gf = 0;
    a.Ff = "";
    a.If = 0;
    a.Te = 0;
    a.Ue = 0;
    a.Mf = 0;
    a.Lf = 0;
    a.Kf = "";
    a.Nf = 0;
    ha(a.Zb, 0, " ", zh);
    a.lc = 0;
    a.cf = c;
    a.bf = b || "input";
    nk(a);
    Y(a);
  }
  function Jn(a, b, c) {
    null == c
      ? (a.he = function (a) {
          x(a);
        })
      : ((a.he = c), (a.jh = b));
    a.dd = "";
  }
  function wn(a, b) {
    a.he(b, a.le);
  }
  function Kn(a) {
    0 < a.dd.length && (a.he(a.dd, a.le), (a.dd = ""));
  }
  function U(a, b) {
    var c;
    switch (a.D) {
      case 1:
      case 2:
        c = Error(a.bf + ":" + a.bb + ": " + b);
        c.line = a.bb;
        c.column = a.Vc;
        for (var d; 0 < a.lc; )
          a.lc--,
            (d = a.Zb[0]),
            ga(a.Zb, 0, a.Zb, 1, zh - 1),
            (a.Zb[zh - 1] = d);
        x(
          "Context: " +
            a.bb +
            " > " +
            (" " == a.Zb[0] ? "" : "...") +
            a.Zb.join("").trim()
        );
        break;
      case 3:
        d = null == a.ib ? 0 : a.ib.bb;
        var e = null == a.ib ? 0 : a.ib.Vc;
        c = Error(d + ": " + b);
        c.line = d;
        c.column = e;
    }
    a.D = 4;
    throw c;
  }
  function ok(a, b) {
    switch (a.D) {
      case 1:
      case 2:
        x(a.bf + ":" + a.bb + ": warning: " + b);
        break;
      case 3:
        x(a.Yf + ":" + (null == a.ib ? 0 : a.ib.bb) + ": warning: " + b);
    }
  }
  var Ld = (exports.mpl_initialize = function () {
      var a = {
        bb: 0,
        Vc: 0,
        l: 0,
        b: 0,
        Bb: 0,
        h: "",
        value: 0,
        Hf: 0,
        Gf: 0,
        Ff: "",
        If: 0,
        Te: 0,
        Ue: 0,
        Mf: 0,
        Lf: 0,
        Kf: "",
        Nf: 0,
      };
      a.Zb = Array(zh);
      ha(a.Zb, 0, " ", zh);
      a.lc = 0;
      a.nc = 0;
      a.V = {};
      a.uc = null;
      a.Pf = 0;
      a.rg = 0;
      a.qg = 0;
      a.Ke = 0;
      a.Ob = 0;
      a.pg = null;
      a.Ah = "";
      a.Bh = "";
      a.Gd = sg();
      a.Of = 0;
      a.ib = null;
      a.Mc = null;
      a.g = 0;
      a.i = 0;
      a.n = null;
      a.f = null;
      a.cf = null;
      a.bf = null;
      a.he = null;
      a.jh = null;
      a.Hg = null;
      a.le = null;
      a.D = 0;
      a.Yf = null;
      a.xh = "";
      return a;
    }),
    Nd = (exports.mpl_read_model = function (a, b, c, d) {
      function e() {
        x(a.bb + " line" + (1 == a.bb ? "" : "s") + " were read");
        a.cf = null;
        return a.D;
      }
      0 != a.D && w("mpl_read_model: invalid call sequence");
      null == c && w("mpl_read_model: no input specified");
      a.D = 1;
      x("Reading model section from " + b + " ...");
      In(a, b, c);
      fl(a);
      null == a.uc && U(a, "empty model section not allowed");
      a.Yf = a.bf;
      Hn(a);
      if (rk(a, "data")) {
        if (d) return ok(a, "data section ignored"), e();
        a.nc = 1;
        Y(a);
        a.b != ni && U(a, "semicolon missing where expected");
        Y(a);
        a.D = 2;
        x("Reading data section from " + b + " ...");
        Kl(a);
      }
      cl(a);
      return e();
    }),
    Pd = (exports.mpl_read_data = function (a, b, c) {
      1 != a.D && 2 != a.D && w("mpl_read_data: invalid call sequence");
      null == c && w("mpl_read_data: no input specified");
      a.D = 2;
      x("Reading data section from " + b + " ...");
      a.nc = 1;
      In(a, b, c);
      dl(a, "data") &&
        (Y(a), a.b != ni && U(a, "semicolon missing where expected"), Y(a));
      Kl(a);
      cl(a);
      x(a.bb + " line" + (1 == a.bb ? "" : "s") + " were read");
      a.cf = null;
      return a.D;
    }),
    Rd = (exports.mpl_generate = function (a, b, c, d) {
      1 != a.D && 2 != a.D && w("mpl_generate: invalid call sequence");
      a.D = 3;
      a.ve = d;
      Jn(a, b, c);
      for (b = a.uc; null != b && (Fn(a, b), a.ib.type != vh); b = b.e);
      a.ib = b;
      Kn(a);
      for (b = a.uc; null != b; b = b.e)
        if (b.type == yh) for (c = b.v.t, c = c.O.head; null != c; c = c.e);
      for (b = a.uc; null != b; b = b.e)
        if (b.type == ch)
          for (c = b.v.H, c = c.O.head; null != c; c = c.e)
            for (c.value.H.ea = ++a.g, d = c.value.H.form; null != d; d = d.e)
              d.t.ga.value.t.C = -1;
      for (b = a.uc; null != b; b = b.e)
        if (b.type == yh)
          for (c = b.v.t, c = c.O.head; null != c; c = c.e)
            0 != c.value.t.C && (c.value.t.C = ++a.i);
      a.n = Array(1 + a.g);
      for (d = 1; d <= a.g; d++) a.n[d] = null;
      for (b = a.uc; null != b; b = b.e)
        if (b.type == ch)
          for (c = b.v.H, c = c.O.head; null != c; c = c.e)
            (d = c.value.H.ea), (a.n[d] = c.value.H);
      for (d = 1; d <= a.g; d++);
      a.f = Array(1 + a.i);
      for (d = 1; d <= a.i; d++) a.f[d] = null;
      for (b = a.uc; null != b; b = b.e)
        if (b.type == yh)
          for (c = b.v.t, c = c.O.head; null != c; c = c.e)
            (d = c.value.t.C), 0 != d && (a.f[d] = c.value.t);
      for (d = 1; d <= a.i; d++);
      x("Model has been successfully generated");
      return a.D;
    }),
    Sd = (exports.mpl_get_prob_name = function (a) {
      return a.Yf;
    }),
    Td = (exports.mpl_get_num_rows = function (a) {
      3 != a.D && w("mpl_get_num_rows: invalid call sequence");
      return a.g;
    }),
    be = (exports.mpl_get_num_cols = function (a) {
      3 != a.D && w("mpl_get_num_cols: invalid call sequence");
      return a.i;
    }),
    Ud = (exports.mpl_get_row_name = function (a, b) {
      3 != a.D && w("mpl_get_row_name: invalid call sequence");
      (1 <= b && b <= a.g) ||
        w("mpl_get_row_name: i = " + b + "; row number out of range");
      var c = a.n[b].H.name,
        c = c + zl("[", a.n[b].ga.w).slice(0, 255);
      255 == c.length && (c = c.slice(0, 252) + "...");
      return c;
    }),
    ie = (exports.mpl_get_row_kind = function (a, b) {
      var c;
      3 != a.D && w("mpl_get_row_kind: invalid call sequence");
      (1 <= b && b <= a.g) ||
        w("mpl_get_row_kind: i = " + b + "; row number out of range");
      switch (a.n[b].H.type) {
        case ch:
          c = 411;
          break;
        case ph:
          c = je;
          break;
        case oh:
          c = ke;
      }
      return c;
    }),
    Vd = (exports.mpl_get_row_bnds = function (a, b, c) {
      var d;
      3 != a.D && w("mpl_get_row_bnds: invalid call sequence");
      (1 <= b && b <= a.g) ||
        w("mpl_get_row_bnds: i = " + b + "; row number out of range");
      d = a.n[b];
      a = null == d.H.P ? -s : d.P;
      b = null == d.H.W ? +s : d.W;
      a == -s && b == +s
        ? ((d = Wd), (a = b = 0))
        : b == +s
        ? ((d = Xd), (b = 0))
        : a == -s
        ? ((d = Yd), (a = 0))
        : (d = d.H.P != d.H.W ? Zd : $d);
      c(a, b);
      return d;
    }),
    he = (exports.mpl_get_mat_row = function (a, b, c, d) {
      var e = 0;
      3 != a.D && w("mpl_get_mat_row: invalid call sequence");
      (1 <= b && b <= a.g) ||
        w("mpl_get_mat_row: i = " + b + "; row number out of range");
      for (a = a.n[b].form; null != a; a = a.e)
        e++, null != c && (c[e] = a.t.C), null != d && (d[e] = a.u);
      return e;
    }),
    ae = (exports.mpl_get_row_c0 = function (a, b) {
      var c;
      3 != a.D && w("mpl_get_row_c0: invalid call sequence");
      (1 <= b && b <= a.g) ||
        w("mpl_get_row_c0: i = " + b + "; row number out of range");
      c = a.n[b];
      return null == c.H.P && null == c.H.W ? -c.P : 0;
    }),
    ce = (exports.mpl_get_col_name = function (a, b) {
      3 != a.D && w("mpl_get_col_name: invalid call sequence");
      (1 <= b && b <= a.i) ||
        w("mpl_get_col_name: j = " + b + "; column number out of range");
      var c = a.f[b].t.name,
        c = c + zl("[", a.f[b].ga.w);
      255 == c.length && (c = c.slice(0, 252) + "...");
      return c;
    }),
    de = (exports.mpl_get_col_kind = function (a, b) {
      var c;
      3 != a.D && w("mpl_get_col_kind: invalid call sequence");
      (1 <= b && b <= a.i) ||
        w("mpl_get_col_kind: j = " + b + "; column number out of range");
      switch (a.f[b].t.type) {
        case N:
          c = 421;
          break;
        case mh:
          c = ee;
          break;
        case ah:
          c = fe;
      }
      return c;
    }),
    ge = (exports.mpl_get_col_bnds = function (a, b, c) {
      var d;
      3 != a.D && w("mpl_get_col_bnds: invalid call sequence");
      (1 <= b && b <= a.i) ||
        w("mpl_get_col_bnds: j = " + b + "; column number out of range");
      d = a.f[b];
      a = null == d.t.P ? -s : d.P;
      b = null == d.t.W ? +s : d.W;
      a == -s && b == +s
        ? ((d = Wd), (a = b = 0))
        : b == +s
        ? ((d = Xd), (b = 0))
        : a == -s
        ? ((d = Yd), (a = 0))
        : (d = d.t.P != d.t.W ? Zd : $d);
      c(a, b);
      return d;
    }),
    me = (exports.mpl_has_solve_stmt = function (a) {
      3 != a.D && w("mpl_has_solve_stmt: invalid call sequence");
      return a.Ob;
    }),
    ne = (exports.mpl_put_row_soln = function (a, b, c, d, e) {
      a.n[b].m = c;
      a.n[b].r = d;
      a.n[b].J = e;
    }),
    oe = (exports.mpl_put_col_soln = function (a, b, c, d, e) {
      a.f[b].m = c;
      a.f[b].r = d;
      a.f[b].J = e;
    }),
    pe = (exports.mpl_postsolve = function (a) {
      (3 != a.D || a.Of) && w("mpl_postsolve: invalid call sequence");
      var b;
      a.Of = 1;
      for (b = a.ib; null != b; b = b.e) Fn(a, b);
      a.ib = null;
      Kn(a);
      x("Model has been successfully processed");
      return a.D;
    }),
    Ln = "Monday Tuesday Wednesday Thursday Friday Saturday Sunday".split(" "),
    Mn =
      "January February March April May June July August September October November December".split(
        " "
      );
  function Nn(a) {
    for (var b = ""; 0 < a; ) (b += "^"), a--;
    return b;
  }
  function On(a, b, c, d, e, f) {
    x("Input string passed to str2time:");
    x(b);
    x(Nn(c + 1));
    x("Format string passed to str2time:\n");
    x(d);
    x(Nn(e + 1));
    U(a, f);
  }
  function cn(a, b, c) {
    function d() {
      On(a, b, n, c, t, "time zone offset value incomplete or invalid");
    }
    function e() {
      On(a, b, n, c, t, "time zone offset value out of range");
    }
    function f() {
      b[n] != c[t] && On(a, b, n, c, t, "character mismatch");
      n++;
    }
    var g, h, k, l, p, m, q, r, n, t;
    h = k = l = p = m = q = -1;
    r = 2147483647;
    for (t = n = 0; t < c.length; t++)
      if ("%" == c[t])
        if ((t++, "b" == c[t] || "h" == c[t])) {
          var y;
          for (
            0 <= k && On(a, b, n, c, t, "month multiply specified");
            " " == b[n];

          )
            n++;
          for (k = 1; 12 >= k; k++) {
            y = Mn[k - 1];
            var E = !1;
            for (g = 0; 2 >= g; g++)
              if (n[g].toUpperCase() != y[g].toUpperCase()) {
                E = !0;
                break;
              }
            if (!E) {
              n += 3;
              for (
                g = 3;
                "\x00" != y[g] && b[n].toUpperCase() == y[g].toUpperCase();
                g++
              )
                n++;
              break;
            }
          }
          12 < k &&
            On(a, b, n, c, t, "abbreviated month name missing or invalid");
        } else if ("d" == c[t]) {
          for (
            0 <= l && On(a, b, n, c, t, "day multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "day missing or invalid");
          l = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (l = 10 * l + (b[n++] - 0));
          (1 <= l && 31 >= l) || On(a, b, n, c, t, "day out of range");
        } else if ("H" == c[t]) {
          for (
            0 <= p && On(a, b, n, c, t, "hour multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "hour missing or invalid");
          p = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (p = 10 * p + (b[n++] - 0));
          (0 <= p && 23 >= p) || On(a, b, n, c, t, "hour out of range");
        } else if ("m" == c[t]) {
          for (
            0 <= k && On(a, b, n, c, t, "month multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "month missing or invalid");
          k = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (k = 10 * k + (b[n++] - 0));
          (1 <= k && 12 >= k) || On(a, b, n, c, t, "month out of range");
        } else if ("M" == c[t]) {
          for (
            0 <= m && On(a, b, n, c, t, "minute multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "minute missing or invalid");
          m = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (m = 10 * m + (b[n++] - 0));
          (0 <= m && 59 >= m) || On(a, b, n, c, t, "minute out of range");
        } else if ("S" == c[t]) {
          for (
            0 <= q && On(a, b, n, c, t, "second multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "second missing or invalid");
          q = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (q = 10 * q + (b[n++] - 0));
          (0 <= q && 60 >= q) || On(a, b, n, c, t, "second out of range");
        } else if ("y" == c[t]) {
          for (
            0 <= h && On(a, b, n, c, t, "year multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "year missing or invalid");
          h = b[n++] - 0;
          "0" <= b[n] && "9" >= b[n] && (h = 10 * h + (b[n++] - 0));
          h += 69 <= h ? 1900 : 2e3;
        } else if ("Y" == c[t]) {
          for (
            0 <= h && On(a, b, n, c, t, "year multiply specified");
            " " == b[n];

          )
            n++;
          ("0" <= b[n] && "9" >= b[n]) ||
            On(a, b, n, c, t, "year missing or invalid");
          h = 0;
          for (g = 1; 4 >= g && "0" <= b[n] && "9" >= b[n]; g++)
            h = 10 * h + (b[n++] - 0);
          (1 <= h && 4e3 >= h) || On(a, b, n, c, t, "year out of range");
        } else if ("z" == c[t]) {
          var C;
          for (
            2147483647 != r &&
            On(a, b, n, c, t, "time zone offset multiply specified");
            " " == b[n];

          )
            n++;
          if ("Z" == b[n]) (C = p = m = 0), n++;
          else {
            "+" == b[n]
              ? ((C = 1), n++)
              : "-" == b[n]
              ? ((C = -1), n++)
              : On(a, b, n, c, t, "time zone offset sign missing");
            p = 0;
            for (g = 1; 2 >= g; g++)
              ("0" <= b[n] && "9" >= b[n]) || d(), (p = 10 * p + (b[n++] - 0));
            23 < p && e();
            ":" == b[n] && (n++, ("0" <= b[n] && "9" >= b[n]) || d());
            m = 0;
            if ("0" <= b[n] && "9" >= b[n]) {
              for (g = 1; 2 >= g; g++)
                ("0" <= b[n] && "9" >= b[n]) || d(),
                  (m = 10 * m + (b[n++] - 0));
              59 < m && e();
            }
          }
          r = C * (60 * p + m);
        } else
          "%" == c[t] ? f() : On(a, b, n, c, t, "invalid conversion specifier");
      else " " != c[t] && f();
    0 > h && (h = 1970);
    0 > k && (k = 1);
    0 > l && (l = 1);
    0 > p && (p = 0);
    0 > m && (m = 0);
    0 > q && (q = 0);
    2147483647 == r && (r = 0);
    g = xg(l, k, h);
    return 60 * (60 * (24 * (g - xg(1, 1, 1970)) + p) + m) + q - 60 * r;
  }
  function Pn(a, b, c) {
    x("Format string passed to time2str:");
    x(b);
    x(Nn(c));
    U(a, "invalid conversion specifier");
  }
  function Qn(a) {
    return ((a + xg(1, 1, 1970)) % 7) + 1;
  }
  function Rn(a) {
    a = xg(1, 1, a) - xg(1, 1, 1970);
    switch (Qn(a)) {
      case 1:
        a += 0;
        break;
      case 2:
        a -= 1;
        break;
      case 3:
        a -= 2;
        break;
      case 4:
        a -= 3;
        break;
      case 5:
        a += 3;
        break;
      case 6:
        a += 2;
        break;
      case 7:
        a += 1;
    }
    Qn(a);
    return a;
  }
  function dn(a, b, c) {
    var d,
      e = 0,
      f = 0,
      g = 0,
      h,
      k,
      l,
      p = "",
      m;
    (-62135596800 <= b && 64092211199 >= b) ||
      U(a, "time2str(" + b + ",...); argument out of range");
    b = Math.floor(b + 0.5);
    h = Math.abs(b) / 86400;
    d = Math.floor(h);
    0 > b && (d = h == Math.floor(h) ? -d : -(d + 1));
    yg(d + xg(1, 1, 1970), function (a, b, c) {
      g = a;
      f = b;
      e = c;
    });
    k = (b - 86400 * d) | 0;
    h = k / 60;
    k %= 60;
    b = h / 60;
    h %= 60;
    for (l = 0; l < c.length; l++)
      "%" == c[l]
        ? (l++,
          "a" == c[l]
            ? (m = Ln[Qn(d) - 1].slice(0, 3))
            : "A" == c[l]
            ? (m = Ln[Qn(d) - 1])
            : "b" == c[l] || "h" == c[l]
            ? (m = Mn[f - 1].slice(0, 3))
            : "B" == c[l]
            ? (m = Mn[f - 1])
            : "C" == c[l]
            ? (m = String(Math.floor(e / 100)))
            : "d" == c[l]
            ? (m = String(g))
            : "D" == c[l]
            ? (m = f + "/" + g + "/" + (e % 100))
            : "e" == c[l]
            ? (m = String(g))
            : "F" == c[l]
            ? xa(m, e + "-" + f + "-" + g)
            : "g" == c[l]
            ? ((m = d < Rn(e) ? e - 1 : d < Rn(e + 1) ? e : e + 1),
              (m = String(m % 100)))
            : "G" == c[l]
            ? ((m = d < Rn(e) ? e - 1 : d < Rn(e + 1) ? e : e + 1),
              (m = String(m)))
            : "H" == c[l]
            ? (m = String(b))
            : "I" == c[l]
            ? (m = String(0 == b ? 12 : 12 >= b ? b : b - 12))
            : "j" == c[l]
            ? (m = String(xg(g, f, e) - xg(1, 1, e) + 1))
            : "k" == c[l]
            ? (m = String(b))
            : "l" == c[l]
            ? (m = String(0 == b ? 12 : 12 >= b ? b : b - 12))
            : "m" == c[l]
            ? (m = String(f))
            : "M" == c[l]
            ? (m = String(h))
            : "p" == c[l]
            ? (m = 11 >= b ? "AM" : "PM")
            : "P" == c[l]
            ? (m = 11 >= b ? "am" : "pm")
            : "r" == c[l]
            ? (m =
                (0 == b ? 12 : 12 >= b ? b : b - 12) +
                ":" +
                h +
                ":" +
                k +
                " " +
                (11 >= b ? "AM" : "PM"))
            : "R" == c[l]
            ? (m = b + ":" + h)
            : "S" == c[l]
            ? (m = String(k))
            : "T" == c[l]
            ? (m = b + ":" + h + ":" + k)
            : "u" == c[l]
            ? (m = String(Qn(d)))
            : "U" == c[l]
            ? ((m = xg(1, 1, e) - xg(1, 1, 1970)),
              (m += 7 - Qn(m)),
              (m = String((d + 7 - m) / 7)))
            : "V" == c[l]
            ? ((m =
                d < Rn(e)
                  ? d - Rn(e - 1)
                  : d < Rn(e + 1)
                  ? d - Rn(e)
                  : d - Rn(e + 1)),
              (m = String(m / 7 + 1)))
            : "w" == c[l]
            ? (m = String(Qn(d) % 7))
            : "W" == c[l]
            ? ((m = xg(1, 1, e) - xg(1, 1, 1970)),
              (m += (8 - Qn(m)) % 7),
              (m = String((d + 7 - m) / 7)))
            : "y" == c[l]
            ? (m = String(e % 100))
            : "Y" == c[l]
            ? (m = String(e))
            : "%" == c[l]
            ? (m = "%")
            : Pn(a, c, l))
        : (m = c[l]),
        (p += m);
    return p;
  }
  var Sn = {};
  function Gn(a, b) {
    var c = a.Mc,
      d = Sn[c.a[1].toLowerCase()];
    d
      ? (c.link = new d(c, b, a.ve))
      : U(a, "Invalid table driver '" + c.a[1] + "'");
    null == c.link && U(a, "error on opening table " + a.ib.v.nd.name);
  }
  var Tn = (exports.mpl_tab_drv_register = function (a, b) {
    Sn[a.toLowerCase()] = b;
  });
  function Un(a, b, c) {
    this.mode = b;
    this.Ea = null;
    this.count = 0;
    this.l = "\n";
    this.Hb = 0;
    this.pb = "";
    this.Wa = 0;
    this.Oa = [];
    this.ve = c;
    this.ng = 0;
    this.Ae = 1;
    this.og = 2;
    this.Be = 3;
    2 > kn(a) && w("csv_driver: file name not specified\n");
    this.Ea = ln(a, 2);
    if ("R" == b) {
      c
        ? ((this.data = c(a.a, b)), (this.cursor = 0))
        : w("csv_driver: unable to open " + this.Ea);
      this.Eg = 0;
      for (Vn(this); ; ) {
        Vn(this);
        if (this.Hb == this.Ae) break;
        this.Hb != this.Be &&
          w(this.Ea + ":" + this.count + ": invalid field name\n");
        this.Wa++;
        for (b = mn(a); 1 <= b && nn(a, b) != this.pb; b--);
        this.Oa[this.Wa] = b;
      }
      for (b = mn(a); 1 <= b && "RECNO" != nn(a, b); b--);
      this.Oa[0] = b;
    } else if ("W" == b) {
      this.data = "";
      c = mn(a);
      for (b = 1; b <= c; b++) this.data += nn(a, b) + (b < c ? "," : "\n");
      this.count++;
    }
  }
  function Vn(a) {
    if (-1 == a.l) (a.Hb = a.ng), (a.pb = "EOF");
    else if ("\n" == a.l) {
      if (
        ((a.Hb = a.Ae),
        (a.pb = "EOR"),
        Wn(a),
        "," == a.l && w(a.Ea + ":" + a.count + ": empty field not allowed\n"),
        "\n" == a.l && w(a.Ea + ":" + a.count + ": empty record not allowed\n"),
        "#" == a.l && 1 == a.count)
      )
        for (; "#" == a.l; ) {
          for (; "\n" != a.l; ) Wn(a);
          Wn(a);
          a.Eg++;
        }
    } else if (("," == a.l && Wn(a), "'" == a.l || '"' == a.l)) {
      var b = a.l;
      a.pb = "";
      a.Hb = a.Be;
      for (Wn(a); ; ) {
        if (a.l == b && (Wn(a), a.l != b))
          if ("," == a.l || "\n" == a.l) break;
          else w(a.Ea + ":" + a.count + ": invalid field");
        a.pb += a.l;
        Wn(a);
      }
      0 == a.pb.length && w(a.Ea + ":" + a.count + ": empty field not allowed");
    } else {
      a.pb = "";
      for (a.Hb = a.og; "," != a.l && "\n" != a.l; )
        ("'" != a.l && '"' != a.l) ||
          w(
            a.Ea +
              ":" +
              a.count +
              ": invalid use of single or double quote within field"
          ),
          (a.pb += a.l),
          Wn(a);
      0 == a.pb.length && w(a.Ea + ":" + a.count + ": empty field not allowed");
      vg(a.pb, function () {}) && (a.Hb = a.Be);
    }
  }
  function Wn(a) {
    var b;
    for ("\n" == a.l && a.count++; ; )
      if (
        ((b = a.cursor < a.data.length ? a.data[a.cursor++] : -1), "\r" != b)
      ) {
        "\n" != b &&
          ta(b) &&
          w(a.Ea + ":" + a.count + ": invalid control character " + b);
        break;
      }
    a.l = b;
  }
  Un.prototype.readRecord = function (a) {
    var b;
    0 < this.Oa[0] && rn(a, this.Oa[0], this.count - this.Eg - 1);
    for (b = 1; b <= this.Wa; b++) {
      Vn(this);
      if (this.Hb == this.ng) return -1;
      if (this.Hb == this.Ae) {
        var c = this.Wa - b + 1;
        1 == c
          ? w(this.Ea + ":" + this.count + ": one field missing")
          : w(this.Ea + ":" + this.count + ": " + c + " fields missing");
      } else if (this.Hb == this.og) {
        if (0 < this.Oa[b]) {
          var d = 0;
          vg(this.pb, function (a) {
            d = a;
          });
          rn(a, this.Oa[b], d);
        }
      } else this.Hb == this.Be && 0 < this.Oa[b] && sn(a, this.Oa[b], this.pb);
    }
    Vn(this);
    this.Hb != this.Ae && w(this.Ea + ":" + this.count + ": too many fields");
    return 0;
  };
  Un.prototype.writeRecord = function (a) {
    var b, c, d, e;
    c = mn(a);
    for (b = 1; b <= c; b++) {
      switch (on(a, b)) {
        case "N":
          this.data += pn(a, b);
          break;
        case "S":
          this.data += '"';
          d = qn(a, b);
          for (e = 0; d.length > e; e++)
            this.data = '"' == d[e] ? this.data + '""' : this.data + d[e];
          this.data += '"';
      }
      this.data += b < c ? "," : "\n";
    }
    this.count++;
    return 0;
  };
  Un.prototype.flush = function (a) {
    this.ve(a.a, this.mode, this.data);
  };
  Tn("CSV", Un);
  function Xn(a, b, c) {
    this.mode = b;
    this.Ea = null;
    2 > kn(a) && w("json driver: file name not specified");
    this.Ea = ln(a, 2);
    if ("R" == b)
      for (
        this.Oa = {},
          c
            ? ((this.data = c(a.a, b)),
              "string" == typeof this.data &&
                (this.data = JSON.parse(this.data)),
              (this.cursor = 1))
            : w("json driver: unable to open " + this.Ea),
          a = 0,
          b = this.data[0];
        a < b.length;
        a++
      )
        this.Oa[b[a]] = a;
    else if ("W" == b) {
      this.ve = c;
      c = [];
      this.data = [c];
      var d = mn(a);
      for (b = 1; b <= d; b++) c.push(nn(a, b));
    }
  }
  Xn.prototype.writeRecord = function (a) {
    var b,
      c = mn(a),
      d = [];
    for (b = 1; b <= c; b++)
      switch (on(a, b)) {
        case "N":
          d.push(pn(a, b));
          break;
        case "S":
          d.push(qn(a, b));
      }
    this.data.push(d);
    return 0;
  };
  Xn.prototype.readRecord = function (a) {
    var b = this.data[this.cursor++];
    if (null == b) return -1;
    for (var c = 1; c <= mn(a); c++) {
      var d = this.Oa[nn(a, c)];
      if (null != d)
        switch (((d = b[d]), typeof d)) {
          case "number":
            rn(a, c, d);
            break;
          case "boolean":
            rn(a, c, Number(d));
            break;
          case "string":
            sn(a, c, d);
            break;
          default:
            w("Unexpected data type " + d + " in " + this.Ea);
        }
    }
    return 0;
  };
  Xn.prototype.flush = function (a) {
    this.ve(a.a, this.mode, this.data);
  };
  Tn("JSON", Xn);
  function Xb() {
    var a = { $b: 0 };
    a.wc = a.gh = a.hh = 0;
    a.name = a.eb = null;
    a.ha = 0;
    a.ee = a.de = 0;
    a.Db = a.Pc = null;
    a.Kb = a.td = null;
    a.top = null;
    a.g = a.i = a.L = 0;
    a.uf = a.Sd = null;
    a.da = a.ue = 0;
    a.nf = a.xg = a.Kg = a.zg = 0;
    a.la = null;
    a.Aa = null;
    a.ka = null;
    a.Pa = null;
    return a;
  }
  function Yn(a, b, c) {
    0 == c
      ? ((b.ca = null),
        (b.e = a.Db),
        null == b.e ? (a.Pc = b) : (b.e.ca = b),
        (a.Db = b))
      : ((b.ca = a.Pc),
        (b.e = null),
        null == b.ca ? (a.Db = b) : (b.ca.e = b),
        (a.Pc = b));
  }
  function Zn(a, b) {
    null == b.ca ? (a.Db = b.e) : (b.ca.e = b.e);
    null == b.e ? (a.Pc = b.ca) : (b.e.ca = b.ca);
  }
  function $n(a, b) {
    b.ja || ((b.ja = 1), Zn(a, b), Yn(a, b, 0));
  }
  function ao(a, b, c) {
    0 == c
      ? ((b.ca = null),
        (b.e = a.Kb),
        null == b.e ? (a.td = b) : (b.e.ca = b),
        (a.Kb = b))
      : ((b.ca = a.td),
        (b.e = null),
        null == b.ca ? (a.Kb = b) : (b.ca.e = b),
        (a.td = b));
  }
  function bo(a, b) {
    null == b.ca ? (a.Kb = b.e) : (b.ca.e = b.e);
    null == b.e ? (a.td = b.ca) : (b.e.ca = b.ca);
  }
  function co(a, b) {
    b.ja || ((b.ja = 1), bo(a, b), ao(a, b, 0));
  }
  function eo(a) {
    var b = {};
    b.ea = ++a.ee;
    b.name = null;
    b.c = -s;
    b.d = +s;
    b.k = null;
    b.ja = 0;
    Yn(a, b, 1);
    return b;
  }
  function fo(a) {
    var b = {};
    b.C = ++a.de;
    b.name = null;
    b.Ra = 0;
    b.c = b.d = b.u = 0;
    b.k = null;
    b.ja = 0;
    b.qb = {};
    b.ub = {};
    ao(a, b, 1);
    return b;
  }
  function go(a, b, c) {
    var d = {};
    d.n = a;
    d.f = b;
    d.j = c;
    d.ua = null;
    d.B = a.k;
    d.ra = null;
    d.I = b.k;
    null != d.B && (d.B.ua = d);
    null != d.I && (d.I.ra = d);
    a.k = b.k = d;
  }
  function ho(a, b) {
    var c;
    c = {};
    c.Zd = b;
    c.info = {};
    c.link = a.top;
    a.top = c;
    return c.info;
  }
  function io(a) {
    for (var b; null != a.k; )
      (b = a.k),
        (a.k = b.B),
        null == b.ra ? (b.f.k = b.I) : (b.ra.I = b.I),
        null != b.I && (b.I.ra = b.ra);
  }
  function jo(a, b) {
    io(b);
    Zn(a, b);
  }
  function ko(a, b) {
    for (var c; null != b.k; )
      (c = b.k),
        (b.k = c.I),
        null == c.ua ? (c.n.k = c.B) : (c.ua.B = c.B),
        null != c.B && (c.B.ua = c.ua);
    bo(a, b);
  }
  function Yb(a, b, c) {
    var d = cb,
      e = cb,
      f = b.g,
      g = b.i,
      h,
      k,
      l;
    a.$b = b.dir;
    a.$b == za ? (l = 1) : a.$b == Ea && (l = -1);
    a.wc = f;
    a.gh = g;
    a.hh = b.L;
    d && null != b.name && (a.name = b.name);
    d && null != b.eb && (a.eb = b.eb);
    a.ha = l * b.ha;
    h = Array(1 + f);
    for (k = 1; k <= f; k++) {
      var p = b.n[k],
        m;
      h[k] = m = eo(a);
      d && null != p.name && (m.name = p.name);
      if (e) {
        var q = p.ma;
        p.type == Ka
          ? ((m.c = -s), (m.d = +s))
          : p.type == Sa
          ? ((m.c = p.c * q), (m.d = +s))
          : p.type == Ta
          ? ((m.c = -s), (m.d = p.d * q))
          : p.type == I
          ? ((m.c = p.c * q), (m.d = p.d * q))
          : p.type == B && (m.c = m.d = p.c * q);
      } else
        p.type == Ka
          ? ((m.c = -s), (m.d = +s))
          : p.type == Sa
          ? ((m.c = p.c), (m.d = +s))
          : p.type == Ta
          ? ((m.c = -s), (m.d = p.d))
          : p.type == I
          ? ((m.c = p.c), (m.d = p.d))
          : p.type == B && (m.c = m.d = p.c);
    }
    for (f = 1; f <= g; f++)
      if (
        ((m = b.f[f]),
        (k = fo(a)),
        d && null != m.name && (k.name = m.name),
        c == Sc && (k.Ra = Number(m.kind == Fc)),
        e)
      )
        for (
          p = m.va,
            m.type == Ka
              ? ((k.c = -s), (k.d = +s))
              : m.type == Sa
              ? ((k.c = m.c / p), (k.d = +s))
              : m.type == Ta
              ? ((k.c = -s), (k.d = m.d / p))
              : m.type == I
              ? ((k.c = m.c / p), (k.d = m.d / p))
              : m.type == B && (k.c = k.d = m.c / p),
            k.u = l * m.u * p,
            m = m.k;
          null != m;
          m = m.I
        )
          go(h[m.n.ea], k, m.n.ma * m.j * p);
      else
        for (
          m.type == Ka
            ? ((k.c = -s), (k.d = +s))
            : m.type == Sa
            ? ((k.c = m.c), (k.d = +s))
            : m.type == Ta
            ? ((k.c = -s), (k.d = m.d))
            : m.type == I
            ? ((k.c = m.c), (k.d = m.d))
            : m.type == B && (k.c = k.d = m.c),
            k.u = l * m.u,
            m = m.k;
          null != m;
          m = m.I
        )
          go(h[m.n.ea], k, m.j);
    a.da = c;
    a.ue = e;
  }
  function cc(a, b) {
    var c, d, e, f, g, h, k;
    db(b);
    Ca(b, a.name);
    Da(b, a.eb);
    Fa(b, a.$b);
    a.$b == za ? (h = 1) : a.$b == Ea && (h = -1);
    Xa(b, 0, h * a.ha);
    for (c = a.Db; null != c; c = c.e)
      (c.ja = e = La(b, 1)),
        Pa(b, e, c.name),
        (d =
          c.c == -s && c.d == +s
            ? Ka
            : c.d == +s
            ? Sa
            : c.c == -s
            ? Ta
            : c.c != c.d
            ? I
            : B),
        Va(b, e, d, c.c, c.d);
    g = new Int32Array(1 + b.g);
    k = new Float64Array(1 + b.g);
    for (c = a.Kb; null != c; c = c.e) {
      e = Oa(b, 1);
      Qa(b, e, c.name);
      Hc(b, e, c.Ra ? Fc : Ma);
      d =
        c.c == -s && c.d == +s
          ? Ka
          : c.d == +s
          ? Sa
          : c.c == -s
          ? Ta
          : c.c != c.d
          ? I
          : B;
      Wa(b, e, d, c.c, c.d);
      Xa(b, e, h * c.u);
      f = 0;
      for (d = c.k; null != d; d = d.I) f++, (g[f] = d.n.ja), (k[f] = d.j);
      Za(b, e, f, g, k);
    }
    a.g = b.g;
    a.i = b.i;
    a.L = b.L;
    a.uf = new Int32Array(1 + a.g);
    a.Sd = new Int32Array(1 + a.i);
    c = a.Db;
    for (e = 0; null != c; c = c.e) a.uf[++e] = c.ea;
    c = a.Kb;
    for (e = 0; null != c; c = c.e) a.Sd[++e] = c.C;
    a.name = a.eb = null;
    a.ha = 0;
    a.Db = a.Pc = null;
    a.Kb = a.td = null;
  }
  function Ub(a, b) {
    var c, d, e, f;
    a.$b == za ? (d = 1) : a.$b == Ea && (d = -1);
    a.da == Zb
      ? ((a.nf = b.na), (a.xg = b.sa))
      : a.da == le
      ? (a.Kg = b.df)
      : a.da == Sc && (a.zg = b.za);
    if (a.da == Zb) {
      null == a.la && (a.la = new Int8Array(1 + a.ee));
      for (f = 1; f <= a.ee; f++) a.la[f] = 0;
      null == a.ka && (a.ka = new Int8Array(1 + a.de));
      for (c = 1; c <= a.de; c++) a.ka[c] = 0;
    }
    null == a.Pa && (a.Pa = new Float64Array(1 + a.de));
    for (c = 1; c <= a.de; c++) a.Pa[c] = s;
    if (a.da != Sc)
      for (
        null == a.Aa && (a.Aa = new Float64Array(1 + a.ee)), f = 1;
        f <= a.ee;
        f++
      )
        a.Aa[f] = s;
    if (a.da == Zb) {
      for (f = 1; f <= a.g; f++)
        (c = b.n[f]), (e = a.uf[f]), (a.la[e] = c.m), (a.Aa[e] = d * c.J);
      for (c = 1; c <= a.i; c++)
        (d = b.f[c]), (e = a.Sd[c]), (a.ka[e] = d.m), (a.Pa[e] = d.r);
    } else if (a.da == le) {
      for (f = 1; f <= a.g; f++)
        (c = b.n[f]), (e = a.uf[f]), (a.Aa[e] = d * c.mc);
      for (c = 1; c <= a.i; c++) (d = b.f[c]), (e = a.Sd[c]), (a.Pa[e] = d.dc);
    } else if (a.da == Sc)
      for (c = 1; c <= a.i; c++) (d = b.f[c]), (e = a.Sd[c]), (a.Pa[e] = d.Sa);
    for (e = a.top; null != e; e = e.link) e.Zd(a, e.info);
  }
  function Vb(a, b) {
    var c, d, e, f;
    a.$b == za ? (e = 1) : a.$b == Ea && (e = -1);
    if (a.da == Zb) {
      b.valid = 0;
      b.na = a.nf;
      b.sa = a.xg;
      b.aa = b.ha;
      b.some = 0;
      for (d = 1; d <= b.g; d++)
        (c = b.n[d]),
          (c.m = a.la[d]),
          (c.J = a.ue ? e * a.Aa[d] * c.ma : e * a.Aa[d]),
          c.m == A
            ? (c.J = 0)
            : c.m == G
            ? (c.r = c.c)
            : c.m == Ua
            ? (c.r = c.d)
            : c.m == Ra
            ? (c.r = 0)
            : c.m == Na && (c.r = c.c);
      for (d = 1; d <= b.i; d++)
        (c = b.f[d]),
          (c.m = a.ka[d]),
          (c.r = a.ue ? a.Pa[d] * c.va : a.Pa[d]),
          c.m == A
            ? (c.J = 0)
            : c.m == G
            ? (c.r = c.c)
            : c.m == Ua
            ? (c.r = c.d)
            : c.m == Ra
            ? (c.r = 0)
            : c.m == Na && (c.r = c.c),
          (b.aa += c.u * c.r);
      for (d = 1; d <= b.g; d++)
        if (((c = b.n[d]), c.m == A)) {
          f = 0;
          for (e = c.k; null != e; e = e.B) f += e.j * e.f.r;
          c.r = f;
        }
      for (d = 1; d <= b.i; d++)
        if (((c = b.f[d]), c.m != A)) {
          f = c.u;
          for (e = c.k; null != e; e = e.I) f -= e.j * e.n.J;
          c.J = f;
        }
    } else if (a.da == le) {
      b.df = a.Kg;
      b.ae = b.ha;
      for (d = 1; d <= b.g; d++)
        (c = b.n[d]), (c.mc = a.ue ? e * a.Aa[d] * c.ma : e * a.Aa[d]);
      for (d = 1; d <= b.i; d++)
        (c = b.f[d]),
          (c.dc = a.ue ? a.Pa[d] * c.va : a.Pa[d]),
          (b.ae += c.u * c.dc);
      for (d = 1; d <= b.g; d++) {
        c = b.n[d];
        f = 0;
        for (e = c.k; null != e; e = e.B) f += e.j * e.f.dc;
        c.dc = f;
      }
      for (d = 1; d <= b.i; d++) {
        c = b.f[d];
        f = c.u;
        for (e = c.k; null != e; e = e.I) f -= e.j * e.n.mc;
        c.mc = f;
      }
    } else if (a.da == Sc) {
      b.za = a.zg;
      b.ta = b.ha;
      for (d = 1; d <= b.i; d++)
        (c = b.f[d]), (c.Sa = a.Pa[d]), (b.ta += c.u * c.Sa);
      for (d = 1; d <= b.g; d++) {
        c = b.n[d];
        f = 0;
        for (e = c.k; null != e; e = e.B) f += e.j * e.f.Sa;
        c.Sa = f;
      }
    }
  }
  function lo(a, b) {
    ho(a, function (a, b) {
      a.da == Zb && (a.la[b.s] = A);
      a.da != Sc && (a.Aa[b.s] = 0);
      return 0;
    }).s = b.ea;
    jo(a, b);
  }
  function mo(a, b) {
    var c, d;
    c = ho(a, function (a, b) {
      if (a.da == Zb)
        if (a.ka[b.F] == A || a.ka[b.F] == G || a.ka[b.F] == Ua)
          a.ka[b.F] = a.ka[b.F];
        else return 1;
      a.Pa[b.F] = b.Ng + a.Pa[b.F];
      return 0;
    });
    c.F = b.C;
    c.Ng = b.c;
    a.ha += b.u * b.c;
    for (d = b.k; null != d; d = d.I)
      (c = d.n),
        c.c == c.d
          ? (c.d = c.c -= d.j * b.c)
          : (c.c != -s && (c.c -= d.j * b.c), c.d != +s && (c.d -= d.j * b.c));
    b.d != +s && (b.d -= b.c);
    b.c = 0;
  }
  function no(a, b) {
    var c, d;
    c = ho(a, function (a, b) {
      a.da == Zb && (a.ka[b.F] = Na);
      a.Pa[b.F] = b.ph;
      return 0;
    });
    c.F = b.C;
    c.ph = b.c;
    a.ha += b.u * b.c;
    for (d = b.k; null != d; d = d.I)
      (c = d.n),
        c.c == c.d
          ? (c.d = c.c -= d.j * b.c)
          : (c.c != -s && (c.c -= d.j * b.c), c.d != +s && (c.d -= d.j * b.c));
    ko(a, b);
  }
  function oo(a, b) {
    var c, d, e;
    d = 1e-9 + 1e-12 * Math.abs(b.c);
    b.d - b.c > d ||
      ((ho(a, function (a, b) {
        if (a.da == Zb)
          if (a.la[b.s] == A) a.la[b.s] = A;
          else if (a.la[b.s] == Na) a.la[b.s] = 0 <= a.Aa[b.s] ? G : Ua;
          else return 1;
        return 0;
      }).s = b.ea),
      (c = 0.5 * (b.d + b.c)),
      (e = Math.floor(c + 0.5)),
      Math.abs(c - e) <= d && (c = e),
      (b.c = b.d = c));
  }
  function po(a, b) {
    var c, d, e, f;
    f = 1e-9 + 1e-12 * Math.abs(b.c);
    if (b.d - b.c > f) return 0;
    c = ho(a, function (a, b) {
      var c, d;
      if (a.da == Zb)
        if (a.ka[b.F] == A) a.ka[b.F] = A;
        else if (a.ka[b.F] == Na) {
          d = b.l;
          for (c = b.k; null != c; c = c.e) d -= c.j * a.Aa[c.Oa];
          a.ka[b.F] = 0 <= d ? G : Ua;
        } else return 1;
      return 0;
    });
    c.F = b.C;
    c.l = b.u;
    c.k = null;
    if (a.da == Zb)
      for (d = b.k; null != d; d = d.I)
        (e = {}), (e.Oa = d.n.ea), (e.j = d.j), (e.e = c.k), (c.k = e);
    c = 0.5 * (b.d + b.c);
    d = Math.floor(c + 0.5);
    Math.abs(c - d) <= f && (c = d);
    b.c = b.d = c;
    return 1;
  }
  function qo(a, b) {
    if (0.001 < b.c || -0.001 > b.d) return 1;
    b.c = -s;
    b.d = +s;
    lo(a, b);
    return 0;
  }
  function ro(a, b) {
    function c() {
      e.m = G;
      b.d = b.c;
    }
    function d() {
      e.m = Ua;
      b.c = b.d;
    }
    var e;
    if ((0.001 < b.u && b.c == -s) || (-0.001 > b.u && b.d == +s)) return 1;
    e = ho(a, function (a, b) {
      a.da == Zb && (a.ka[b.F] = b.m);
      return 0;
    });
    e.F = b.C;
    b.c == -s && b.d == +s
      ? ((e.m = Ra), (b.c = b.d = 0))
      : b.d == +s
      ? c()
      : b.c == -s
      ? d()
      : b.c != b.d
      ? 2.220446049250313e-16 <= b.u
        ? c()
        : -2.220446049250313e-16 >= b.u
        ? d()
        : Math.abs(b.c) <= Math.abs(b.d)
        ? c()
        : d()
      : (e.m = Na);
    no(a, b);
    return 0;
  }
  function so(a, b) {
    var c;
    if (a.Ra)
      if (((c = Math.floor(b + 0.5)), 1e-5 >= Math.abs(b - c))) b = c;
      else return 2;
    if (a.c != -s) {
      c = a.Ra ? 1e-5 : 1e-5 + 1e-8 * Math.abs(a.c);
      if (b < a.c - c) return 1;
      if (b < a.c + 0.001 * c) return (a.d = a.c), 0;
    }
    if (a.d != +s) {
      c = a.Ra ? 1e-5 : 1e-5 + 1e-8 * Math.abs(a.d);
      if (b > a.d + c) return 1;
      if (b > a.d - 0.001 * c) return (a.c = a.d), 0;
    }
    a.c = a.d = b;
    return 0;
  }
  function to(a, b) {
    var c, d, e;
    e = b.k;
    d = e.f;
    c = so(d, b.c / e.j);
    if (0 != c) return c;
    c = ho(a, function (a, b) {
      var c, d;
      if (a.da == Zb) {
        if (a.ka[b.F] != Na) return 1;
        a.la[b.s] = Na;
        a.ka[b.F] = A;
      }
      if (a.da != Sc) {
        d = b.l;
        for (c = b.k; null != c; c = c.e) d -= c.j * a.Aa[c.Oa];
        a.Aa[b.s] = d / b.Da;
      }
      return 0;
    });
    c.s = b.ea;
    c.F = d.C;
    c.Da = e.j;
    c.l = d.u;
    c.k = null;
    if (a.da != Sc)
      for (e = d.k; null != e; e = e.I)
        e.n != b &&
          ((d = {}), (d.Oa = e.n.ea), (d.j = e.j), (d.e = c.k), (c.k = d));
    jo(a, b);
    return 0;
  }
  function uo(a, b) {
    var c;
    a.Ra &&
      ((c = Math.floor(b + 0.5)),
      (b = 1e-5 >= Math.abs(b - c) ? c : Math.ceil(b)));
    if (
      a.c != -s &&
      ((c = a.Ra ? 0.001 : 0.001 + 1e-6 * Math.abs(a.c)), b < a.c + c)
    )
      return 0;
    if (a.d != +s) {
      c = a.Ra ? 1e-5 : 1e-5 + 1e-8 * Math.abs(a.d);
      if (b > a.d + c) return 4;
      if (b > a.d - 0.001 * c) return (a.c = a.d), 3;
    }
    c =
      a.c == -s
        ? 2
        : a.Ra && b > a.c + 0.5
        ? 2
        : b > a.c + 0.3 * (1 + Math.abs(a.c))
        ? 2
        : 1;
    a.c = b;
    return c;
  }
  function vo(a, b) {
    var c;
    a.Ra &&
      ((c = Math.floor(b + 0.5)),
      (b = 1e-5 >= Math.abs(b - c) ? c : Math.floor(b)));
    if (
      a.d != +s &&
      ((c = a.Ra ? 0.001 : 0.001 + 1e-6 * Math.abs(a.d)), b > a.d - c)
    )
      return 0;
    if (a.c != -s) {
      c = a.Ra ? 1e-5 : 1e-5 + 1e-8 * Math.abs(a.c);
      if (b < a.c - c) return 4;
      if (b < a.c + 0.001 * c) return (a.d = a.c), 3;
    }
    c =
      a.d == +s
        ? 2
        : a.Ra && b < a.d - 0.5
        ? 2
        : b < a.d - 0.3 * (1 + Math.abs(a.d))
        ? 2
        : 1;
    a.d = b;
    return c;
  }
  function wo(a, b) {
    var c, d, e, f, g, h;
    e = b.k;
    d = e.f;
    0 < e.j
      ? ((g = b.c == -s ? -s : b.c / e.j), (c = b.d == +s ? +s : b.d / e.j))
      : ((g = b.d == +s ? -s : b.d / e.j), (c = b.c == -s ? +s : b.c / e.j));
    if (g == -s) g = 0;
    else if (((g = uo(d, g)), 4 == g)) return 4;
    if (c == +s) h = 0;
    else if (3 == g) h = 0;
    else if (((h = vo(d, c)), 4 == h)) return 4;
    if (!g && !h) return (b.c = -s), (b.d = +s), lo(a, b), 0;
    c = ho(a, function (a, b) {
      var c, d;
      if (a.da == Sc) return 0;
      d = b.l;
      for (c = b.k; null != c; c = c.e) d -= c.j * a.Aa[c.Oa];
      if (a.da == Zb) {
        c = function () {
          b.Vf
            ? ((a.la[b.s] = 0 < b.Da ? G : Ua),
              (a.ka[b.F] = A),
              (a.Aa[b.s] = d / b.Da))
            : ((a.la[b.s] = A), (a.Aa[b.s] = 0));
          return 0;
        };
        var e = function () {
          b.gg
            ? ((a.la[b.s] = 0 < b.Da ? Ua : G),
              (a.ka[b.F] = A),
              (a.Aa[b.s] = d / b.Da))
            : ((a.la[b.s] = A), (a.Aa[b.s] = 0));
          return 0;
        };
        if (a.ka[b.F] == A) (a.la[b.s] = A), (a.Aa[b.s] = 0);
        else if (a.ka[b.F] == G) c();
        else if (a.ka[b.F] == Ua) e();
        else if (a.ka[b.F] == Na) {
          if (
            1e-7 < d &&
            ((0 < b.Da && b.c != -s) || (0 > b.Da && b.d != +s) || !b.Vf)
          )
            return (a.ka[b.F] = G), c();
          if (
            -1e-7 > d &&
            ((0 < b.Da && b.d != +s) || (0 > b.Da && b.c != -s) || !b.gg)
          )
            return (a.ka[b.F] = Ua), e();
          if (b.c != -s && b.d == +s) a.la[b.s] = G;
          else if (b.c == -s && b.d != +s) a.la[b.s] = Ua;
          else if (b.c != -s && b.d != +s)
            a.la[b.s] = b.Da * a.Pa[b.F] <= 0.5 * (b.c + b.d) ? G : Ua;
          else return 1;
          a.ka[b.F] = A;
          a.Aa[b.s] = d / b.Da;
        } else return 1;
      }
      a.da == le &&
        (a.Aa[b.s] =
          (2.220446049250313e-16 < d && b.Vf) ||
          (-2.220446049250313e-16 > d && b.gg)
            ? d / b.Da
            : 0);
      return 0;
    });
    c.s = b.ea;
    c.F = d.C;
    c.Da = e.j;
    c.l = d.u;
    c.c = b.c;
    c.d = b.d;
    c.Vf = g;
    c.gg = h;
    c.k = null;
    if (a.da != Sc)
      for (d = d.k; null != d; d = d.I)
        d != e &&
          ((f = {}), (f.Oa = d.n.ea), (f.j = d.j), (f.e = c.k), (c.k = f));
    jo(a, b);
    return g >= h ? g : h;
  }
  function xo(a, b) {
    var c, d, e, f;
    e = b.k;
    d = e.n;
    c = ho(a, function (a, b) {
      var c, d;
      if (a.da == Zb) {
        if (a.la[b.s] == A || a.la[b.s] == Ra) a.ka[b.F] = a.la[b.s];
        else if (a.la[b.s] == G) a.ka[b.F] = 0 < b.Da ? Ua : G;
        else if (a.la[b.s] == Ua) a.ka[b.F] = 0 < b.Da ? G : Ua;
        else return 1;
        a.la[b.s] = Na;
      }
      a.da != Sc && (a.Aa[b.s] += b.l / b.Da);
      c = b.qd;
      for (d = b.k; null != d; d = d.e) c -= d.j * a.Pa[d.Oa];
      a.Pa[b.F] = c / b.Da;
      return 0;
    });
    c.s = d.ea;
    c.F = b.C;
    c.Da = e.j;
    c.qd = d.c;
    c.l = b.u;
    c.k = null;
    for (e = d.k; null != e; e = e.B)
      e.f != b &&
        ((f = {}),
        (f.Oa = e.f.C),
        (f.j = e.j),
        (f.e = c.k),
        (c.k = f),
        (e.f.u -= (e.j / c.Da) * c.l));
    a.ha += (c.qd / c.Da) * c.l;
    0 < c.Da
      ? ((d.c = b.d == +s ? -s : c.qd - c.Da * b.d),
        (d.d = b.c == -s ? +s : c.qd - c.Da * b.c))
      : ((d.c = b.c == -s ? -s : c.qd - c.Da * b.c),
        (d.d = b.d == +s ? +s : c.qd - c.Da * b.d));
    ko(a, b);
  }
  function yo(a, b) {
    function c() {
      e.m = G;
      f.d = f.c;
    }
    function d() {
      e.m = Ua;
      f.c = f.d;
    }
    var e, f, g, h, k, l;
    g = b.k;
    f = g.n;
    k = f.c;
    if (k != -s)
      for (h = f.k; null != h; h = h.B)
        if (h != g)
          if (0 < h.j) {
            if (h.f.d == +s) {
              k = -s;
              break;
            }
            k -= h.j * h.f.d;
          } else {
            if (h.f.c == -s) {
              k = -s;
              break;
            }
            k -= h.j * h.f.c;
          }
    l = f.d;
    if (l != +s)
      for (h = f.k; null != h; h = h.B)
        if (h != g)
          if (0 < h.j) {
            if (h.f.c == -s) {
              l = +s;
              break;
            }
            l -= h.j * h.f.c;
          } else {
            if (h.f.d == +s) {
              l = +s;
              break;
            }
            l -= h.j * h.f.d;
          }
    h = 0 < g.j ? (k == -s ? -s : k / g.j) : l == +s ? -s : l / g.j;
    k = 0 < g.j ? (l == +s ? +s : l / g.j) : k == -s ? +s : k / g.j;
    if (
      (b.c != -s && ((l = 1e-9 + 1e-12 * Math.abs(b.c)), h < b.c - l)) ||
      (b.d != +s && ((l = 1e-9 + 1e-12 * Math.abs(b.d)), k > b.d + l))
    )
      return 1;
    b.c = -s;
    b.d = +s;
    e = ho(a, function (a, b) {
      if (a.da == Zb)
        if (a.la[b.s] == A) a.la[b.s] = A;
        else if (a.la[b.s] == Na) a.la[b.s] = b.m;
        else return 1;
      return 0;
    });
    e.s = f.ea;
    e.m = -1;
    g = b.u / g.j;
    if (2.220446049250313e-16 < g)
      if (f.c != -s) c();
      else {
        if (1e-5 < g) return 2;
        d();
      }
    else if (-2.220446049250313e-16 > g)
      if (f.d != +s) d();
      else {
        if (-1e-5 > g) return 2;
        c();
      }
    else
      f.d == +s
        ? c()
        : f.c == -s
        ? d()
        : Math.abs(f.c) <= Math.abs(f.d)
        ? c()
        : d();
    return 0;
  }
  function zo(a, b, c) {
    var d,
      e = null,
      f,
      g,
      h;
    d = 1;
    for (g = b.k; null != g; g = g.B) d < Math.abs(g.j) && (d = Math.abs(g.j));
    for (g = b.k; null != g; g = g.B) if (Math.abs(g.j) < 1e-7 * d) return 1;
    d = ho(a, function (a, b) {
      var c, d, e, f, g;
      if (a.da == Sc) return 0;
      if (a.da == Zb) {
        if (a.la[b.s] != A) return 1;
        for (c = b.k; null != c; c = c.e) {
          if (a.ka[c.C] != Na) return 1;
          a.ka[c.C] = c.m;
        }
      }
      for (c = b.k; null != c; c = c.e) {
        e = c.l;
        for (d = c.k; null != d; d = d.e) e -= d.j * a.Aa[d.Oa];
        c.l = e;
      }
      d = null;
      f = 0;
      for (c = b.k; null != c; c = c.e)
        if (((e = c.l), (g = Math.abs(e / c.Jc)), c.m == G))
          0 > e && f < g && ((d = c), (f = g));
        else if (c.m == Ua) 0 < e && f < g && ((d = c), (f = g));
        else return 1;
      null != d &&
        (a.da == Zb && ((a.la[b.s] = b.m), (a.ka[d.C] = A)),
        (a.Aa[b.s] = d.l / d.Jc));
      return 0;
    });
    d.s = b.ea;
    d.m = b.c == b.d ? Na : 0 == c ? G : Ua;
    d.k = null;
    for (g = b.k; null != g; g = g.B)
      if (
        ((f = g.f),
        a.da != Sc &&
          ((e = {}),
          (e.C = f.C),
          (e.m = -1),
          (e.Jc = g.j),
          (e.l = f.u),
          (e.k = null),
          (e.e = d.k),
          (d.k = e)),
        (0 == c && 0 > g.j) || (0 != c && 0 < g.j)
          ? (a.da != Sc && (e.m = G), (f.d = f.c))
          : (a.da != Sc && (e.m = Ua), (f.c = f.d)),
        a.da != Sc)
      )
        for (f = f.k; null != f; f = f.I)
          f != g &&
            ((h = {}), (h.Oa = f.n.ea), (h.j = f.j), (h.e = e.k), (e.k = h));
    b.c = -s;
    b.d = +s;
    return 0;
  }
  function Ao(a) {
    var b,
      c = 0,
      d,
      e;
    d = 0;
    for (b = a.k; null != b; b = b.B)
      if (0 < b.j) {
        if (b.f.c == -s) {
          d = -s;
          break;
        }
        d += b.j * b.f.c;
      } else {
        if (b.f.d == +s) {
          d = -s;
          break;
        }
        d += b.j * b.f.d;
      }
    e = 0;
    for (b = a.k; null != b; b = b.B)
      if (0 < b.j) {
        if (b.f.d == +s) {
          e = +s;
          break;
        }
        e += b.j * b.f.d;
      } else {
        if (b.f.c == -s) {
          e = +s;
          break;
        }
        e += b.j * b.f.c;
      }
    if (
      (a.c != -s && ((b = 0.001 + 1e-6 * Math.abs(a.c)), a.c - b > e)) ||
      (a.d != +s && ((b = 0.001 + 1e-6 * Math.abs(a.d)), a.d + b < d))
    )
      return 51;
    a.c != -s &&
      ((b = 1e-9 + 1e-12 * Math.abs(a.c)),
      a.c - b > d && (c = a.c + b <= e ? c | 1 : c | 2));
    a.d != +s &&
      ((b = 1e-9 + 1e-12 * Math.abs(a.d)),
      a.d + b < e && (c = a.d - b >= d ? c | 16 : c | 32));
    return c;
  }
  function Bo(a, b, c) {
    a.da == Zb &&
      ((a = ho(a, function (a, b) {
        if (a.da != Zb) return 1;
        a.la[b.s] = a.la[b.s] == A ? A : b.m;
        return 0;
      })),
      (a.s = b.ea),
      (a.m =
        b.d == +s ? G : b.c == -s ? Ua : b.c != b.d ? (0 == c ? Ua : G) : Na));
    0 == c ? (b.c = -s) : 1 == c && (b.d = +s);
  }
  function Co(a) {
    var b, c, d, e, f, g, h, k, l, p, m;
    k = l = p = m = 0;
    for (d = a.td; null != d; d = d.ca)
      if (d.Ra && d.c != d.d && (0 != d.c || 1 != d.d))
        if (-1e6 > d.c || 1e6 < d.d || 4095 < d.d - d.c) k++;
        else if ((l++, 0 != d.c && mo(a, d), (e = d.d | 0), 1 != e)) {
          g = 2;
          for (c = 4; e >= c; ) g++, (c += c);
          p += g;
          b = ho(a, function (a, b) {
            var c,
              d,
              e = a.Pa[b.F];
            c = 1;
            for (d = 2; c < b.i; c++, d += d) e += d * a.Pa[b.C + (c - 1)];
            a.Pa[b.F] = e;
            return 0;
          });
          b.F = d.C;
          b.C = 0;
          b.i = g;
          e < c - 1 ? ((c = eo(a)), m++, (c.c = -s), (c.d = e)) : (c = null);
          d.d = 1;
          null != c && go(c, d, 1);
          h = 1;
          for (c = 2; h < g; h++, c += c)
            for (
              e = fo(a),
                e.Ra = 1,
                e.c = 0,
                e.d = 1,
                e.u = c * d.u,
                0 == b.C && (b.C = e.C),
                f = d.k;
              null != f;
              f = f.I
            )
              go(f.n, e, c * f.j);
        }
    0 < l &&
      x(l + " integer variable(s) were replaced by " + p + " binary ones");
    0 < m && x(m + " row(s) were added due to binarization");
    0 < k && x("Binarization failed for " + k + " integer variable(s)");
  }
  function Do(a, b) {
    var c, d, e;
    d = null;
    for (c = a.k; null != c; c = c.B)
      (e = {}), (e.fa = b * c.j), (e.jc = c.f), (e.e = d), (d = e);
    return d;
  }
  function Eo(a, b, c) {
    var d, e, f;
    for (d = a; null != d; d = d.e);
    e = 0;
    for (d = a; null != d; d = d.e)
      if (1 != d.fa)
        if (-1 == d.fa) e++;
        else break;
    if (null == d && b == 1 - e) return 1;
    for (d = a; null != d; d = d.e) 0 > d.fa && (b -= d.fa);
    for (d = a; null != d; d = d.e) if (Math.abs(d.fa) > b) return 0;
    e = null;
    for (d = a; null != d; d = d.e)
      if (null == e || Math.abs(e.fa) > Math.abs(d.fa)) e = d;
    f = null;
    for (d = a; null != d; d = d.e)
      d != e && (null == f || Math.abs(f.fa) > Math.abs(d.fa)) && (f = d);
    if (Math.abs(e.fa) + Math.abs(f.fa) <= b + (0.001 + 1e-6 * Math.abs(b)))
      return 0;
    b = 1;
    for (d = a; null != d; d = d.e)
      0 < d.fa ? (d.fa = 1) : ((d.fa = -1), (b -= 1));
    c(b);
    return 2;
  }
  function Fo(a, b) {
    var c,
      d,
      e,
      f,
      g = 0,
      h;
    for (f = 0; 1 >= f; f++) {
      if (0 == f) {
        if (b.d == +s) continue;
        e = Do(b, 1);
        h = +b.d;
      } else {
        if (b.c == -s) continue;
        e = Do(b, -1);
        h = -b.c;
      }
      c = Eo(e, h, function (a) {
        h = a;
      });
      if ((1 == f && 1 == c) || 2 == c) {
        g++;
        if (b.c == -s || b.d == +s) c = null;
        else
          for (
            c = eo(a),
              0 == f ? ((c.c = b.c), (c.d = +s)) : ((c.c = -s), (c.d = b.d)),
              d = b.k;
            null != d;
            d = d.B
          )
            go(c, d.f, d.j);
        io(b);
        b.c = -s;
        for (b.d = h; null != e; e = e.e) go(b, e.jc, e.fa);
        null != c && (b = c);
      }
    }
    return g;
  }
  function Go(a, b, c) {
    var d, e;
    for (d = a; null != d; d = d.e);
    e = 0;
    for (d = a; null != d; d = d.e)
      if (1 != d.fa)
        if (-1 == d.fa) e++;
        else break;
    if (null == d && b == 1 - e) return 1;
    for (d = a; null != d; d = d.e) 0 > d.fa && (b -= d.fa);
    if (0.001 > b) return 0;
    e = 1e-9 + 1e-12 * Math.abs(b);
    for (d = a; null != d; d = d.e) if (Math.abs(d.fa) < b - e) return 0;
    b = 1;
    for (d = a; null != d; d = d.e)
      0 < d.fa ? (d.fa = 1) : ((d.fa = -1), (b -= 1));
    c(b);
    return 2;
  }
  function Ho(a, b) {
    var c,
      d,
      e,
      f,
      g = 0,
      h;
    for (f = 0; 1 >= f; f++) {
      if (0 == f) {
        if (b.c == -s) continue;
        e = Do(b, 1);
        h = +b.c;
      } else {
        if (b.d == +s) continue;
        e = Do(b, -1);
        h = -b.d;
      }
      c = Go(e, h, function (a) {
        h = a;
      });
      if ((1 == f && 1 == c) || 2 == c) {
        g++;
        if (b.c == -s || b.d == +s) c = null;
        else
          for (
            c = eo(a),
              0 == f ? ((c.c = -s), (c.d = b.d)) : ((c.c = b.c), (c.d = +s)),
              d = b.k;
            null != d;
            d = d.B
          )
            go(c, d.f, d.j);
        io(b);
        b.c = h;
        for (b.d = +s; null != e; e = e.e) go(b, e.jc, e.fa);
        null != c && (b = c);
      }
    }
    return g;
  }
  function Io(a, b, c) {
    var d,
      e = 0,
      f,
      g;
    f = 0;
    for (d = a; null != d; d = d.e)
      if (0 < d.fa) {
        if (d.jc.c == -s) return e;
        f += d.fa * d.jc.c;
      } else {
        if (d.jc.d == +s) return e;
        f += d.fa * d.jc.d;
      }
    for (d = a; null != d; d = d.e)
      d.jc.Ra &&
        0 == d.jc.c &&
        1 == d.jc.d &&
        (0 < d.fa
          ? ((a = f),
            b - d.fa < a &&
              a < b &&
              ((g = b - a),
              0.001 <= g && d.fa - g >= 0.01 * (1 + d.fa) && ((d.fa = g), e++)))
          : ((a = f - d.fa),
            b < a &&
              a < b - d.fa &&
              ((g = d.fa + (a - b)),
              -0.001 >= g &&
                g - d.fa >= 0.01 * (1 - d.fa) &&
                ((d.fa = g), (f += a - b), (b = a), e++))));
    c(b);
    return e;
  }
  function Jo(a, b) {
    var c,
      d,
      e,
      f,
      g = Array(2),
      h;
    for (f = g[0] = g[1] = 0; 1 >= f; f++) {
      if (0 == f) {
        if (b.c == -s) continue;
        e = Do(b, 1);
        h = +b.c;
      } else {
        if (b.d == +s) continue;
        e = Do(b, -1);
        h = -b.d;
      }
      g[f] = Io(e, h, function (a) {
        h = a;
      });
      if (0 < g[f]) {
        if (b.c == -s || b.d == +s) c = null;
        else
          for (
            c = eo(a),
              0 == f ? ((c.c = -s), (c.d = b.d)) : ((c.c = b.c), (c.d = +s)),
              d = b.k;
            null != d;
            d = d.B
          )
            go(c, d.f, d.j);
        io(b);
        b.c = h;
        b.d = +s;
        for (d = e; null != d; d = d.e) go(b, d.jc, d.fa);
        null != c && (b = c);
      }
    }
    return g[0] + g[1];
  }
  function Ko(a, b, c) {
    function d() {
      for (f = b.k; null != f; f = g) {
        e = f.f;
        g = f.B;
        for (h = e.k; null != h; h = h.I) $n(a, h.n);
        no(a, e);
      }
      lo(a, b);
      return 0;
    }
    var e, f, g, h, k;
    if (null == b.k) {
      k = qo(a, b);
      if (0 == k) return 0;
      if (1 == k) return ac;
    }
    if (null == b.k.B)
      if (((e = b.k.f), b.c == b.d)) {
        k = to(a, b);
        if (0 == k) {
          for (f = e.k; null != f; f = f.I) $n(a, f.n);
          no(a, e);
          return 0;
        }
        if (1 == k || 2 == k) return ac;
      } else {
        k = wo(a, b);
        if (0 <= k && 3 >= k) {
          co(a, e);
          if (2 <= k) for (f = e.k; null != f; f = f.I) $n(a, f.n);
          3 == k && no(a, e);
          return 0;
        }
        if (4 == k) return ac;
      }
    k = Ao(b);
    if (51 == k) return ac;
    if (0 == (k & 15)) b.c != -s && Bo(a, b, 0);
    else if (1 != (k & 15) && 2 == (k & 15) && 0 == zo(a, b, 0)) return d();
    if (0 == (k & 240)) b.d != +s && Bo(a, b, 1);
    else if (16 != (k & 240) && 32 == (k & 240) && 0 == zo(a, b, 1)) return d();
    if (b.c == -s && b.d == +s) {
      for (f = b.k; null != f; f = f.B) co(a, f.f);
      lo(a, b);
      return 0;
    }
    return a.da == Sc && c && 0 > Lo(a, b, 1) ? ac : 0;
  }
  function Lo(a, b, c) {
    var d,
      e,
      f,
      g,
      h,
      k = 0,
      l;
    h = !1;
    e = 1;
    for (d = b.k; null != d; d = d.B)
      (d.f.qb.qb = -s),
        (d.f.ub.ub = +s),
        e < Math.abs(d.j) && (e = Math.abs(d.j));
    g = 1e-6 * e;
    if (b.c != -s) {
      e = null;
      for (d = b.k; null != d; d = d.B)
        if ((0 < d.j && d.f.d == +s) || (0 > d.j && d.f.c == -s))
          if (null == e) e = d;
          else {
            h = !0;
            break;
          }
      if (!h) {
        h = b.c;
        for (d = b.k; null != d; d = d.B)
          d != e && (h = 0 < d.j ? h - d.j * d.f.d : h - d.j * d.f.c);
        if (null == e)
          for (d = b.k; null != d; d = d.B)
            d.j >= +g
              ? (d.f.qb.qb = d.f.d + h / d.j)
              : d.j <= -g && (d.f.ub.ub = d.f.c + h / d.j);
        else
          e.j >= +g
            ? (e.f.qb.qb = h / e.j)
            : e.j <= -g && (e.f.ub.ub = h / e.j);
      }
    }
    h = !1;
    if (b.d != +s) {
      e = null;
      for (d = b.k; null != d; d = d.B)
        if ((0 < d.j && d.f.c == -s) || (0 > d.j && d.f.d == +s))
          if (null == e) e = d;
          else {
            h = !0;
            break;
          }
      if (!h) {
        h = b.d;
        for (d = b.k; null != d; d = d.B)
          d != e && (h = 0 < d.j ? h - d.j * d.f.c : h - d.j * d.f.d);
        if (null == e)
          for (d = b.k; null != d; d = d.B)
            d.j >= +g
              ? (d.f.ub.ub = d.f.c + h / d.j)
              : d.j <= -g && (d.f.qb.qb = d.f.d + h / d.j);
        else
          e.j >= +g
            ? (e.f.ub.ub = h / e.j)
            : e.j <= -g && (e.f.qb.qb = h / e.j);
      }
    }
    for (e = b.k; null != e; )
      for (d = e.f, e = e.B, g = 0; 1 >= g; g++) {
        f = d.c;
        l = d.d;
        if (0 == g) {
          if (d.qb.qb == -s) continue;
          h = uo(d, d.qb.qb);
        } else {
          if (d.ub.ub == +s) continue;
          h = vo(d, d.ub.ub);
        }
        if (0 == h || 1 == h) (d.c = f), (d.d = l);
        else if (2 == h || 3 == h) {
          k++;
          if (c) for (f = d.k; null != f; f = f.I) f.n != b && $n(a, f.n);
          if (3 == h) {
            no(a, d);
            break;
          }
        } else if (4 == h) return -1;
      }
    return k;
  }
  function Mo(a, b) {
    var c, d, e;
    if (null == b.k) {
      e = ro(a, b);
      if (0 == e) return 0;
      if (1 == e) return bc;
    }
    if (null == b.k.I) {
      c = b.k.n;
      var f = function () {
        xo(a, b);
        if (c.c == -s && c.d == +s) {
          for (d = c.k; null != d; d = d.B) co(a, d.f);
          lo(a, c);
        } else $n(a, c);
        return 0;
      };
      if (c.c == c.d) {
        if (!b.Ra) return f();
      } else if (!b.Ra) {
        e = yo(a, b);
        if (0 == e) return f();
        if (1 != e && 2 == e) return bc;
      }
    }
    return 0;
  }
  function $b(a, b) {
    var c, d, e;
    for (c = a.Db; null != c; c = d)
      (d = c.e), c.c == -s && c.d == +s && lo(a, c);
    for (c = a.Db; null != c; c = d)
      (d = c.e), c.c != -s && c.d != +s && c.c < c.d && oo(a, c);
    for (c = a.Kb; null != c; c = d) (d = c.e), c.c == c.d && no(a, c);
    for (c = a.Kb; null != c; c = d)
      (d = c.e),
        c.c != -s &&
          c.d != +s &&
          c.c < c.d &&
          ((e = po(a, c)), 0 != e && 1 == e && no(a, c));
    for (c = a.Db; null != c; c = c.e) c.ja = 1;
    for (c = a.Kb; null != c; c = c.e) c.ja = 1;
    for (d = 1; d; ) {
      for (d = 0; ; ) {
        c = a.Db;
        if (null == c || !c.ja) break;
        d = a;
        e = c;
        e.ja && ((e.ja = 0), Zn(d, e), Yn(d, e, 1));
        c = Ko(a, c, b);
        if (0 != c) return c;
        d = 1;
      }
      for (;;) {
        c = a.Kb;
        if (null == c || !c.ja) break;
        d = a;
        e = c;
        e.ja && ((e.ja = 0), bo(d, e), ao(d, e, 1));
        c = Mo(a, c);
        if (0 != c) return c;
        d = 1;
      }
    }
    if (a.da == Sc && !b)
      for (c = a.Db; null != c; c = c.e) if (0 > Lo(a, c, 0)) return (c = ac);
    return 0;
  }
  function Tc(a, b) {
    var c, d, e, f, g;
    c = $b(a, 1);
    if (0 != c) return c;
    b.sd && Co(a);
    g = 0;
    for (c = a.Pc; null != c; c = d)
      if (
        ((d = c.ca),
        (c.c != -s || c.d != +s) && c.c != c.d && null != c.k && null != c.k.B)
      ) {
        for (
          f = c.k;
          null != f && ((e = f.f), e.Ra && 0 == e.c && 1 == e.d);
          f = f.B
        );
        null == f && (g += Fo(a, c));
      }
    0 < g && x(g + " hidden packing inequaliti(es) were detected");
    g = 0;
    for (c = a.Pc; null != c; c = d)
      if (
        ((d = c.ca),
        (c.c != -s || c.d != +s) &&
          c.c != c.d &&
          null != c.k &&
          null != c.k.B &&
          null != c.k.B.B)
      ) {
        for (
          f = c.k;
          null != f && ((e = f.f), e.Ra && 0 == e.c && 1 == e.d);
          f = f.B
        );
        null == f && (g += Ho(a, c));
      }
    0 < g && x(g + " hidden covering inequaliti(es) were detected");
    g = 0;
    for (c = a.Pc; null != c; c = d) (d = c.ca), c.c != c.d && (g += Jo(a, c));
    0 < g && x(g + " constraint coefficient(s) were reduced");
    return 0;
  }
  function No(a) {
    var b, c;
    b = 1;
    for (c = 32; 55 >= c; b++, c++) a.nb[b] = (a.nb[b] - a.nb[c]) & 2147483647;
    for (c = 1; 55 >= b; b++, c++) a.nb[b] = (a.nb[b] - a.nb[c]) & 2147483647;
    a.Qf = 54;
    return a.nb[55];
  }
  function sg() {
    var a = {},
      b;
    a.nb = Array(56);
    a.nb[0] = -1;
    for (b = 1; 55 >= b; b++) a.nb[b] = 0;
    a.Qf = 0;
    Md(a, 1);
    return a;
  }
  function Md(a, b) {
    var c,
      d,
      e = 1;
    b = d = (b - 0) & 2147483647;
    a.nb[55] = d;
    for (c = 21; c; c = (c + 21) % 55)
      (a.nb[c] = e),
        (e = (d - e) & 2147483647),
        (b = b & 1 ? 1073741824 + (b >> 1) : b >> 1),
        (e = (e - b) & 2147483647),
        (d = a.nb[c]);
    No(a);
    No(a);
    No(a);
    No(a);
    No(a);
  }
  function dm(a) {
    return 0 <= a.nb[a.Qf] ? a.nb[a.Qf--] : No(a);
  }
  function bn(a) {
    var b;
    do b = dm(a);
    while (2147483648 <= b);
    return b % 16777216;
  }
  function tg(a) {
    a = dm(a) / 2147483647;
    return -0.3 * (1 - a) + 0.7 * a;
  }
  var De = 1,
    Fe = 2,
    We = 1,
    Se = 2,
    Ce = 0,
    Ue = 1e-10;
  function Te(a, b, c) {
    return (b - 1) * a.K + c - (b * (b - 1)) / 2;
  }
  function Oo(a, b, c) {
    var d;
    0 == b
      ? ((b = 1), (d = 0))
      : Math.abs(a) <= Math.abs(b)
      ? ((a = -a / b), (d = 1 / Math.sqrt(1 + a * a)), (b = d * a))
      : ((a = -b / a), (b = 1 / Math.sqrt(1 + a * a)), (d = b * a));
    c(b, d);
  }
  function Ve(a, b, c) {
    for (var d = a.i, e = a.Nb, f = a.v, g, h, k, l, p, m; b < d; b++)
      (l = Te(a, b, b)),
        (h = (b - 1) * a.K + 1),
        (p = (d - 1) * a.K + 1),
        Math.abs(f[l]) < Ue && Math.abs(c[b]) < Ue && (f[l] = c[b] = 0),
        0 != c[b] &&
          Oo(f[l], c[b], function (a, r) {
            g = b;
            for (k = l; g <= d; g++, k++) {
              var n = f[k],
                t = c[g];
              f[k] = a * n - r * t;
              c[g] = r * n + a * t;
            }
            g = 1;
            k = h;
            for (m = p; g <= d; g++, k++, m++)
              (n = e[k]),
                (t = e[m]),
                (e[k] = a * n - r * t),
                (e[m] = r * n + a * t);
          });
    Math.abs(c[d]) < Ue && (c[d] = 0);
    f[Te(a, d, d)] = c[d];
  }
  Ce &&
    (Le = function (a, b) {
      var c = a.i,
        d = a.Nb,
        e = a.v,
        f = a.s,
        g = a.l,
        h,
        k,
        l,
        p,
        m = 0;
      for (h = 1; h <= c; h++)
        for (k = 1; k <= c; k++) {
          p = 0;
          for (l = 1; l <= c; l++)
            p += d[(h - 1) * a.K + l] * g[(l - 1) * a.K + k];
          l = f[k];
          l = h <= l ? e[Te(a, h, l)] : 0;
          p = Math.abs(p - l) / (1 + Math.abs(l));
          m < p && (m = p);
        }
      1e-8 < m && x(b + ": dmax = " + m + "; relative error too large");
    });
  function Ke(a, b, c, d) {
    a.pa < a.i && w("scf_solve_it: singular matrix");
    if (b) {
      b = a.i;
      var e = a.Nb,
        f = a.v,
        g = a.s,
        h = a.ig,
        k,
        l,
        p;
      for (k = 1; k <= b; k++) h[k] = c[g[k] + d];
      for (k = 1; k <= b; k++)
        for (
          l = Te(a, k, k), p = h[k] /= f[l], g = k + 1, l++;
          g <= b;
          g++, l++
        )
          h[g] -= f[l] * p;
      for (g = 1; g <= b; g++) c[g + d] = 0;
      for (k = 1; k <= b; k++)
        for (p = h[k], g = 1, l = (k - 1) * a.K + 1; g <= b; g++, l++)
          c[g + d] += e[l] * p;
    } else {
      b = a.i;
      e = a.Nb;
      f = a.v;
      h = a.s;
      k = a.ig;
      for (var m, g = 1; g <= b; g++) {
        m = 0;
        l = 1;
        for (p = (g - 1) * a.K + 1; l <= b; l++, p++) m += e[p] * c[l + d];
        k[g] = m;
      }
      for (g = b; 1 <= g; g--) {
        m = k[g];
        l = b;
        for (p = Te(a, g, b); l > g; l--, p--) m -= f[p] * k[l];
        k[g] = m / f[p];
      }
      for (g = 1; g <= b; g++) c[h[g] + d] = k[g];
    }
  }
  var gc = (exports.glp_scale_prob = function (a, b) {
    function c(a, b) {
      var c, d, e;
      d = 1;
      for (c = a.n[b].k; null != c; c = c.B)
        if (
          ((e = Math.abs(c.j)),
          (e = e * c.n.ma * c.f.va),
          null == c.ua || d > e)
        )
          d = e;
      return d;
    }
    function d(a, b) {
      var c, d, e;
      d = 1;
      for (c = a.n[b].k; null != c; c = c.B)
        if (
          ((e = Math.abs(c.j)),
          (e = e * c.n.ma * c.f.va),
          null == c.ua || d < e)
        )
          d = e;
      return d;
    }
    function e(a, b) {
      var c, d, e;
      d = 1;
      for (c = a.f[b].k; null != c; c = c.I)
        if (
          ((e = Math.abs(c.j)),
          (e = e * c.n.ma * c.f.va),
          null == c.ra || d > e)
        )
          d = e;
      return d;
    }
    function f(a, b) {
      var c, d, e;
      d = 1;
      for (c = a.f[b].k; null != c; c = c.I)
        if (
          ((e = Math.abs(c.j)),
          (e = e * c.n.ma * c.f.va),
          null == c.ra || d < e)
        )
          d = e;
      return d;
    }
    function g(a) {
      var b, d, e;
      for (b = d = 1; b <= a.g; b++)
        if (((e = c(a, b)), 1 == b || d > e)) d = e;
      return d;
    }
    function h(a) {
      var b, c, e;
      for (b = c = 1; b <= a.g; b++)
        if (((e = d(a, b)), 1 == b || c < e)) c = e;
      return c;
    }
    function k(a, b) {
      var c, e, g;
      for (e = 0; 1 >= e; e++)
        if (e == b)
          for (c = 1; c <= a.g; c++) (g = d(a, c)), zb(a, c, Cb(a, c) / g);
        else for (c = 1; c <= a.i; c++) (g = f(a, c)), Ab(a, c, Db(a, c) / g);
    }
    function l(a) {
      var b, e, f;
      for (b = e = 1; b <= a.g; b++)
        if (((f = d(a, b) / c(a, b)), 1 == b || e < f)) e = f;
      return e;
    }
    function p(a) {
      var b, c, d;
      for (b = c = 1; b <= a.i; b++)
        if (((d = f(a, b) / e(a, b)), 1 == b || c < d)) c = d;
      return c;
    }
    function m(a) {
      var b,
        k,
        m = 0,
        y;
      k = l(a) > p(a);
      for (b = 1; 15 >= b; b++) {
        y = m;
        m = h(a) / g(a);
        if (1 < b && m > 0.9 * y) break;
        y = a;
        for (
          var E = k, C = void 0, D = (C = void 0), H = void 0, D = 0;
          1 >= D;
          D++
        )
          if (D == E)
            for (C = 1; C <= y.g; C++)
              (H = c(y, C) * d(y, C)), zb(y, C, Cb(y, C) / Math.sqrt(H));
          else
            for (C = 1; C <= y.i; C++)
              (H = e(y, C) * f(y, C)), Ab(y, C, Db(y, C) / Math.sqrt(H));
      }
    }
    b & ~(Uc | Vc | Wc | Xc | hc) &&
      w("glp_scale_prob: flags = " + b + "; invalid scaling options");
    b & hc && (b = Uc | Vc | Xc);
    (function (a, b) {
      function c(a, b, d, e) {
        return (
          a + ": min|aij| = " + b + "  max|aij| = " + d + "  ratio = " + e + ""
        );
      }
      var d, e;
      x("Scaling...");
      Eb(a);
      d = g(a);
      e = h(a);
      x(c(" A", d, e, e / d));
      if (
        0.1 <= d &&
        10 >= e &&
        (x("Problem data seem to be well scaled"), b & Xc)
      )
        return;
      b & Uc && (m(a), (d = g(a)), (e = h(a)), x(c("GM", d, e, e / d)));
      b & Vc &&
        (k(a, l(a) > p(a)), (d = g(a)), (e = h(a)), x(c("EQ", d, e, e / d)));
      if (b & Wc) {
        for (d = 1; d <= a.g; d++) zb(a, d, ug(Cb(a, d)));
        for (d = 1; d <= a.i; d++) Ab(a, d, ug(Db(a, d)));
        d = g(a);
        e = h(a);
        x(c("2N", d, e, e / d));
      }
    })(a, b);
  });
  function Pb(a, b) {
    function c(a, b, c, d) {
      var e = a.g,
        f = a.Ca,
        g = a.Ba,
        h = a.Ia;
      a = a.head[b];
      if (a <= e) (e = 1), (c[1] = a), (d[1] = 1);
      else
        for (
          b = f[a - e],
            e = f[a - e + 1] - b,
            ga(c, 1, g, b, e),
            ga(d, 1, h, b, e),
            c = 1;
          c <= e;
          c++
        )
          d[c] = -d[c];
      return e;
    }
    function d(a) {
      var b = od(a.U, a.g, c, a);
      a.valid = 0 == b;
      return b;
    }
    function e(a, b, c) {
      var d = a.g,
        e;
      if (c <= d) {
        var f = Array(2);
        e = Array(2);
        f[1] = c;
        e[1] = 1;
        b = Ne(a.U, b, 1, f, 0, e);
      } else {
        var g = a.Ca,
          f = a.Ba,
          h = a.Ia;
        e = a.jb;
        var k;
        k = g[c - d];
        c = g[c - d + 1];
        g = 0;
        for (d = k; d < c; d++) e[++g] = -h[d];
        b = Ne(a.U, b, g, f, k - 1, e);
      }
      a.valid = 0 == b;
      return b;
    }
    function f(a, b, c) {
      var d = a.g,
        e = a.jb,
        f = a.jb,
        g = a.g,
        h = a.Ca,
        k = a.Ba,
        l = a.Ia,
        m = a.head,
        n,
        p,
        q;
      ga(f, 1, b, 1, g);
      for (b = 1; b <= g; b++)
        if (((q = c[b]), 0 != q))
          if (((n = m[b]), n <= g)) f[n] -= q;
          else
            for (p = h[n - g], n = h[n - g + 1]; p < n; p++)
              f[k[p]] += l[p] * q;
      wd(a.U, e);
      for (a = 1; a <= d; a++) c[a] += e[a];
    }
    function g(a, b, c) {
      var d = a.g,
        e = a.jb,
        f = a.jb,
        g = a.g,
        h = a.Ca,
        k = a.Ba,
        l = a.Ia,
        m = a.head,
        n,
        p,
        q,
        K;
      for (n = 1; n <= g; n++) {
        p = m[n];
        K = b[n];
        if (p <= g) K -= c[p];
        else
          for (q = h[p - g], p = h[p - g + 1]; q < p; q++) K += l[q] * c[k[q]];
        f[n] = K;
      }
      yd(a.U, e);
      for (a = 1; a <= d; a++) c[a] += e[a];
    }
    function h(a, b, c) {
      var d = a.g,
        e = a.De,
        f = a.Kd,
        g = a.Ce,
        h = a.Ee,
        k;
      if (c <= d) (k = e[c] + f[c]++), (g[k] = b), (h[k] = 1);
      else {
        k = a.Ca;
        var l = a.Ba;
        a = a.Ia;
        var m;
        m = k[c - d];
        c = k[c - d + 1];
        for (d = m; d < c; d++)
          (k = l[d]), (k = e[k] + f[k]++), (g[k] = b), (h[k] = -a[d]);
      }
    }
    function k(a, b, c) {
      var d = a.g,
        e = a.De,
        f = a.Kd,
        g = a.Ce,
        h = a.Ee,
        k;
      if (c <= d) {
        for (d = k = e[c]; g[d] != b; d++);
        k += --f[c];
        g[d] = g[k];
        h[d] = h[k];
      } else {
        k = a.Ca;
        a = a.Ba;
        var l, m;
        m = k[c - d];
        for (c = k[c - d + 1]; m < c; m++) {
          l = a[m];
          for (d = k = e[l]; g[d] != b; d++);
          k += --f[l];
          g[d] = g[k];
          h[d] = h[k];
        }
      }
    }
    function l(a, b) {
      var c = a.c,
        d = a.d,
        e = a.m,
        f,
        g;
      f = a.head[a.g + b];
      switch (e[b]) {
        case G:
          g = c[f];
          break;
        case Ua:
          g = d[f];
          break;
        case Ra:
          g = 0;
          break;
        case Na:
          g = c[f];
      }
      return g;
    }
    function p(a, b) {
      var c = a.g,
        d = a.i,
        e = a.Ca,
        g = a.Ba,
        h = a.Ia,
        k = a.head,
        m = a.Gc,
        n,
        p,
        q,
        K;
      for (n = 1; n <= c; n++) m[n] = 0;
      for (n = 1; n <= d; n++)
        if (((p = k[c + n]), (K = l(a, n)), 0 != K))
          if (p <= c) m[p] -= K;
          else
            for (q = e[p - c], p = e[p - c + 1]; q < p; q++)
              m[g[q]] += K * h[q];
      ga(b, 1, m, 1, c);
      wd(a.U, b);
      f(a, m, b);
    }
    function m(a) {
      var b = a.i,
        c = a.Qa,
        d = a.Hc,
        e;
      e = a.g;
      var f = a.u,
        h = a.head,
        k = a.Gc,
        l;
      for (l = 1; l <= e; l++) k[l] = f[h[l]];
      ga(d, 1, k, 1, e);
      yd(a.U, d);
      g(a, k, d);
      for (e = 1; e <= b; e++) {
        f = a.g;
        l = a.u;
        k = h = void 0;
        h = a.head[f + e];
        k = l[h];
        if (h <= f) k -= d[h];
        else {
          l = a.Ca;
          for (
            var m = a.Ba,
              n = a.Ia,
              p = void 0,
              q = void 0,
              p = void 0,
              p = l[h - f],
              q = l[h - f + 1];
            p < q;
            p++
          )
            k += n[p] * d[m[p]];
        }
        c[e] = k;
      }
    }
    function q(a) {
      var b = a.g,
        c = a.i,
        d = a.head,
        e = a.hd,
        f = a.gamma,
        g;
      a.Rb = 1e3;
      ha(e, 1, 0, b + c);
      for (a = 1; a <= c; a++) (g = d[b + a]), (e[g] = 1), (f[a] = 1);
    }
    function r(a, b) {
      var c = a.i,
        d = a.m,
        e = a.Qa,
        f = a.gamma,
        g,
        h,
        k,
        l;
      l = h = 0;
      for (g = 1; g <= c; g++) {
        k = e[g];
        switch (d[g]) {
          case G:
            if (k >= -b) continue;
            break;
          case Ua:
            if (k <= +b) continue;
            break;
          case Ra:
            if (-b <= k && k <= +b) continue;
            break;
          case Na:
            continue;
        }
        k = (k * k) / f[g];
        l < k && ((h = g), (l = k));
      }
      a.F = h;
    }
    function n(a) {
      var b = a.g,
        c = a.yb,
        d = a.Ua,
        e = a.Ua,
        f,
        g;
      g = a.head[b + a.F];
      for (f = 1; f <= b; f++) e[f] = 0;
      if (g <= b) e[g] = -1;
      else {
        var h = a.Ca;
        f = a.Ba;
        var k = a.Ia,
          l;
        l = h[g - b];
        for (g = h[g - b + 1]; l < g; l++) e[f[l]] = k[l];
      }
      wd(a.U, d);
      e = 0;
      for (f = 1; f <= b; f++) 0 != d[f] && (c[++e] = f);
      a.Tb = e;
    }
    function t(a) {
      var b = a.g,
        c = a.yb,
        d = a.Ua,
        e = a.Hc,
        g,
        h;
      h = a.head[b + a.F];
      for (g = 1; g <= b; g++) e[g] = 0;
      if (h <= b) e[h] = -1;
      else {
        var k = a.Ca;
        g = a.Ba;
        var l = a.Ia,
          m;
        m = k[h - b];
        for (h = k[h - b + 1]; m < h; m++) e[g[m]] = l[m];
      }
      f(a, e, d);
      e = 0;
      for (g = 1; g <= b; g++) 0 != d[g] && (c[++e] = g);
      a.Tb = e;
    }
    function y(a, b) {
      var c = a.Tb,
        d = a.yb,
        e = a.Ua,
        f,
        g,
        h;
      g = 0;
      for (f = 1; f <= c; f++) (h = Math.abs(e[d[f]])), g < h && (g = h);
      a.rh = g;
      h = b * (1 + 0.01 * g);
      for (g = 0; g < c; )
        (f = d[c]), Math.abs(e[f]) < h ? c-- : (g++, (d[c] = d[g]), (d[g] = f));
      a.sh = g;
    }
    function E(a, b) {
      var c = a.g,
        d = a.type,
        e = a.c,
        f = a.d,
        g = a.u,
        h = a.head,
        k = a.D,
        l = a.Ha,
        m = a.F,
        n = a.yb,
        p = a.Ua,
        q = a.sh,
        K,
        r,
        v,
        u,
        t,
        y,
        oa,
        z,
        C;
      oa = 0 < a.Qa[m] ? -1 : 1;
      K = h[c + m];
      d[K] == I
        ? ((m = -1), (r = 0), (C = f[K] - e[K]), (t = 1))
        : ((r = m = 0), (C = s), (t = 0));
      for (v = 1; v <= q; v++) {
        c = n[v];
        K = h[c];
        u = oa * p[c];
        if (0 < u)
          if (1 == k && 0 > g[K])
            (y = b * (1 + L * Math.abs(e[K]))),
              (z = (e[K] + y - l[c]) / u),
              (K = G);
          else if (1 == k && 0 < g[K]) continue;
          else if (d[K] == Ta || d[K] == I || d[K] == B)
            (y = b * (1 + L * Math.abs(f[K]))),
              (z = (f[K] + y - l[c]) / u),
              (K = Ua);
          else continue;
        else if (1 == k && 0 < g[K])
          (y = b * (1 + L * Math.abs(f[K]))),
            (z = (f[K] - y - l[c]) / u),
            (K = Ua);
        else if (1 == k && 0 > g[K]) continue;
        else if (d[K] == Sa || d[K] == I || d[K] == B)
          (y = b * (1 + L * Math.abs(e[K]))),
            (z = (e[K] - y - l[c]) / u),
            (K = G);
        else continue;
        0 > z && (z = 0);
        if (C > z || (C == z && t < Math.abs(u)))
          (m = c), (r = K), (C = z), (t = Math.abs(u));
      }
      if (!(0 == b || 0 >= m || 0 == C))
        for (y = C, r = m = 0, C = s, t = 0, v = 1; v <= q; v++) {
          c = n[v];
          K = h[c];
          u = oa * p[c];
          if (0 < u)
            if (1 == k && 0 > g[K]) (z = (e[K] - l[c]) / u), (K = G);
            else if (1 == k && 0 < g[K]) continue;
            else if (d[K] == Ta || d[K] == I || d[K] == B)
              (z = (f[K] - l[c]) / u), (K = Ua);
            else continue;
          else if (1 == k && 0 < g[K]) (z = (f[K] - l[c]) / u), (K = Ua);
          else if (1 == k && 0 > g[K]) continue;
          else if (d[K] == Sa || d[K] == I || d[K] == B)
            (z = (e[K] - l[c]) / u), (K = G);
          else continue;
          0 > z && (z = 0);
          z <= y &&
            t < Math.abs(u) &&
            ((m = c), (r = K), (C = z), (t = Math.abs(u)));
        }
      a.s = m;
      a.nf = 0 < m && d[h[m]] == B ? Na : r;
      a.th = oa * C;
    }
    function C(a, b) {
      var c = a.g,
        d = a.s,
        e;
      for (e = 1; e <= c; e++) b[e] = 0;
      b[d] = 1;
      yd(a.U, b);
    }
    function D(a, b) {
      var c = a.g,
        d = a.s,
        e = a.Hc,
        f;
      for (f = 1; f <= c; f++) e[f] = 0;
      e[d] = 1;
      g(a, e, b);
    }
    function H(a, b) {
      var c = a.g,
        d = a.i,
        e = a.De,
        f = a.Kd,
        g = a.Ce,
        h = a.Ee,
        k = a.Vb,
        l = a.mb,
        m,
        n,
        p,
        q;
      for (m = 1; m <= d; m++) l[m] = 0;
      for (m = 1; m <= c; m++)
        if (((q = b[m]), 0 != q))
          for (n = e[m], p = n + f[m]; n < p; n++) l[g[n]] -= q * h[n];
      c = 0;
      for (m = 1; m <= d; m++) 0 != l[m] && (k[++c] = m);
      a.Bc = c;
    }
    function R(a) {
      var b = a.Ha,
        c = a.F,
        d = a.Tb,
        e = a.yb,
        f = a.Ua,
        g = a.s,
        h = a.th;
      0 < g && (b[g] = l(a, c) + h);
      if (0 != h)
        for (c = 1; c <= d; c++) (a = e[c]), a != g && (b[a] += f[a] * h);
    }
    function V(a) {
      var b = a.u,
        c = a.head,
        d = a.Tb,
        e = a.yb,
        f = a.Ua,
        g,
        h;
      h = b[c[a.g + a.F]];
      for (g = 1; g <= d; g++) (a = e[g]), (h += b[c[a]] * f[a]);
      return h;
    }
    function O(a) {
      var b = a.Qa,
        c = a.F,
        d = a.Bc,
        e = a.Vb;
      a = a.mb;
      var f, g, h;
      h = b[c] /= a[c];
      for (g = 1; g <= d; g++) (f = e[g]), f != c && (b[f] -= a[f] * h);
    }
    function Q(a) {
      var b = a.g,
        c = a.type,
        d = a.Ca,
        e = a.Ba,
        f = a.Ia,
        g = a.head,
        h = a.hd,
        k = a.gamma,
        l = a.F,
        m = a.Tb,
        n = a.yb,
        p = a.Ua,
        q = a.s,
        K = a.Bc,
        r = a.Vb,
        v = a.mb,
        u = a.Hc,
        t,
        y,
        oa,
        z,
        C,
        D;
      a.Rb--;
      z = C = h[g[b + l]] ? 1 : 0;
      for (t = 1; t <= b; t++) u[t] = 0;
      for (y = 1; y <= m; y++)
        (t = n[y]), h[g[t]] ? ((u[t] = t = p[t]), (z += t * t)) : (u[t] = 0);
      yd(a.U, u);
      m = v[l];
      for (y = 1; y <= K; y++)
        if (((a = r[y]), a != l)) {
          t = v[a] / m;
          n = g[b + a];
          if (n <= b) D = u[n];
          else
            for (D = 0, oa = d[n - b], p = d[n - b + 1]; oa < p; oa++)
              D -= f[oa] * u[e[oa]];
          p = k[a] + t * t * z + 2 * t * D;
          t = (h[n] ? 1 : 0) + C * t * t;
          k[a] = p >= t ? p : t;
          2.220446049250313e-16 > k[a] && (k[a] = 2.220446049250313e-16);
        }
      c[g[q]] == B
        ? (k[l] = 1)
        : ((k[l] = z / (m * m)),
          2.220446049250313e-16 > k[l] && (k[l] = 2.220446049250313e-16));
    }
    function F(a) {
      var b = a.g,
        c = a.head,
        d = a.m,
        e = a.F,
        f = a.s;
      a = a.nf;
      var g;
      if (0 > f)
        switch (d[e]) {
          case G:
            d[e] = Ua;
            break;
          case Ua:
            d[e] = G;
        }
      else (g = c[f]), (c[f] = c[b + e]), (c[b + e] = g), (d[e] = a);
    }
    function W(a, b) {
      var c = a.g,
        d = a.i,
        e = a.type,
        f = a.c,
        g = a.d,
        h = a.u,
        k = a.head,
        l = a.Ha,
        m,
        n = 0,
        p;
      b *= 0.9;
      for (m = 1; m <= c + d; m++) h[m] = 0;
      for (d = 1; d <= c; d++) {
        m = k[d];
        if (e[m] == Sa || e[m] == I || e[m] == B)
          (p = b * (1 + L * Math.abs(f[m]))),
            l[d] < f[m] - p && ((h[m] = -1), n++);
        if (e[m] == Ta || e[m] == I || e[m] == B)
          (p = b * (1 + L * Math.abs(g[m]))),
            l[d] > g[m] + p && ((h[m] = 1), n++);
      }
      return n;
    }
    function X(a) {
      var b = a.g,
        c = a.i,
        d = a.u,
        e = a.eb;
      a = a.$a;
      var f;
      for (f = 1; f <= b; f++) d[f] = 0;
      for (f = 1; f <= c; f++) d[b + f] = a * e[f];
    }
    function ca(a, b) {
      var c = a.g,
        d = a.type,
        e = a.c,
        f = a.d,
        g = a.u,
        h = a.head,
        k = a.D,
        l = a.Ha,
        m,
        n,
        p;
      for (m = 1; m <= c; m++)
        if (((n = h[m]), 1 == k && 0 > g[n])) {
          if (((p = b * (1 + L * Math.abs(e[n]))), l[m] > e[n] + p)) return 1;
        } else if (1 == k && 0 < g[n]) {
          if (((p = b * (1 + L * Math.abs(f[n]))), l[m] < f[n] - p)) return 1;
        } else {
          if (d[n] == Sa || d[n] == I || d[n] == B)
            if (((p = b * (1 + L * Math.abs(e[n]))), l[m] < e[n] - p)) return 1;
          if (d[n] == Ta || d[n] == I || d[n] == B)
            if (((p = b * (1 + L * Math.abs(f[n]))), l[m] > f[n] + p)) return 1;
        }
      return 0;
    }
    function ka(a, b) {
      var c = a.g,
        d = a.c,
        e = a.d,
        f = a.u,
        g = a.head,
        h = a.Ha,
        k,
        l,
        m;
      for (k = 1; k <= c; k++)
        if (((l = g[k]), 0 > f[l])) {
          if (((m = b * (1 + L * Math.abs(d[l]))), h[k] < d[l] - m)) return 1;
        } else if (
          0 < f[l] &&
          ((m = b * (1 + L * Math.abs(e[l]))), h[k] > e[l] + m)
        )
          return 1;
      return 0;
    }
    function P(a) {
      var b = a.g,
        c = a.i,
        d = a.eb,
        e = a.head,
        f = a.Ha,
        g,
        h,
        k;
      k = d[0];
      for (g = 1; g <= b; g++) (h = e[g]), h > b && (k += d[h - b] * f[g]);
      for (f = 1; f <= c; f++)
        (h = e[b + f]), h > b && (k += d[h - b] * l(a, f));
      return k;
    }
    function u(a, b, c) {
      var d = a.g,
        e = a.type,
        f = a.c,
        g = a.d,
        h = a.D,
        k = a.head,
        l = a.Ha,
        m,
        n;
      if (
        !(
          b.o < fc ||
          (0 < b.fb && 1e3 * la(a.hc) < b.fb) ||
          a.$ == a.be ||
          (!c && 0 != a.$ % b.bc)
        )
      ) {
        m = n = 0;
        for (b = 1; b <= d; b++)
          (c = k[b]),
            (e[c] == Sa || e[c] == I || e[c] == B) &&
              l[b] < f[c] &&
              (n += f[c] - l[b]),
            (e[c] == Ta || e[c] == I || e[c] == B) &&
              l[b] > g[c] &&
              (n += l[b] - g[c]),
            e[c] == B && m++;
        x(
          (1 == h ? " " : "*") +
            a.$ +
            ": obj = " +
            P(a) +
            "  infeas = " +
            n +
            " (" +
            m +
            ")"
        );
        a.be = a.$;
      }
    }
    function z(a, b, c, d, e) {
      var f = a.g,
        g = a.i,
        h = a.$a,
        k = a.head,
        l = a.m,
        m = a.Ha,
        n = a.Qa;
      b.valid = 1;
      a.valid = 0;
      b.U = a.U;
      a.U = null;
      ga(b.head, 1, k, 1, f);
      b.na = c;
      b.sa = d;
      b.aa = P(a);
      b.$ = a.$;
      b.some = e;
      for (a = 1; a <= f; a++)
        (c = k[a]),
          c <= f
            ? ((c = b.n[c]), (c.m = A), (c.bind = a), (c.r = m[a] / c.ma))
            : ((c = b.f[c - f]), (c.m = A), (c.bind = a), (c.r = m[a] * c.va)),
          (c.J = 0);
      for (m = 1; m <= g; m++)
        if (((c = k[f + m]), c <= f)) {
          c = b.n[c];
          c.m = l[m];
          c.bind = 0;
          switch (l[m]) {
            case G:
              c.r = c.c;
              break;
            case Ua:
              c.r = c.d;
              break;
            case Ra:
              c.r = 0;
              break;
            case Na:
              c.r = c.c;
          }
          c.J = (n[m] * c.ma) / h;
        } else {
          c = b.f[c - f];
          c.m = l[m];
          c.bind = 0;
          switch (l[m]) {
            case G:
              c.r = c.c;
              break;
            case Ua:
              c.r = c.d;
              break;
            case Ra:
              c.r = 0;
              break;
            case Na:
              c.r = c.c;
          }
          c.J = n[m] / c.va / h;
        }
    }
    var L = 0.1,
      v,
      S = 2,
      M = 0,
      ba = 0,
      J = 0,
      ea,
      fa,
      sa;
    v = (function (a) {
      var b = a.g,
        c = a.i;
      a = a.L;
      var d = {};
      d.g = b;
      d.i = c;
      d.type = new Int8Array(1 + b + c);
      d.c = new Float64Array(1 + b + c);
      d.d = new Float64Array(1 + b + c);
      d.u = new Float64Array(1 + b + c);
      d.eb = new Float64Array(1 + c);
      d.Ca = new Int32Array(1 + c + 1);
      d.Ba = new Int32Array(1 + a);
      d.Ia = new Float64Array(1 + a);
      d.head = new Int32Array(1 + b + c);
      d.m = new Int8Array(1 + c);
      d.De = new Int32Array(1 + b + 1);
      d.Kd = new Int32Array(1 + b);
      d.Ce = null;
      d.Ee = null;
      d.Ha = new Float64Array(1 + b);
      d.Qa = new Float64Array(1 + c);
      d.hd = new Int8Array(1 + b + c);
      d.gamma = new Float64Array(1 + c);
      d.yb = new Int32Array(1 + b);
      d.Ua = new Float64Array(1 + b);
      d.Vb = new Int32Array(1 + c);
      d.mb = new Float64Array(1 + c);
      d.jb = new Float64Array(1 + b);
      d.Gc = new Float64Array(1 + b);
      d.Hc = new Float64Array(1 + b);
      d.jg = new Float64Array(1 + b);
      return d;
    })(a);
    (function (a, b) {
      var c = a.g,
        d = a.i,
        e = a.type,
        f = a.c,
        g = a.d,
        k = a.u,
        l = a.eb,
        m = a.Ca,
        n = a.Ba,
        p = a.Ia,
        q = a.head,
        K = a.m,
        r = a.hd,
        v = a.gamma,
        t,
        u;
      for (t = 1; t <= c; t++)
        (u = b.n[t]),
          (e[t] = u.type),
          (f[t] = u.c * u.ma),
          (g[t] = u.d * u.ma),
          (k[t] = 0);
      for (t = 1; t <= d; t++)
        (u = b.f[t]),
          (e[c + t] = u.type),
          (f[c + t] = u.c / u.va),
          (g[c + t] = u.d / u.va),
          (k[c + t] = u.u * u.va);
      l[0] = b.ha;
      ga(l, 1, k, c + 1, d);
      e = 0;
      for (t = 1; t <= d; t++) e < Math.abs(l[t]) && (e = Math.abs(l[t]));
      0 == e && (e = 1);
      switch (b.dir) {
        case za:
          a.$a = 1 / e;
          break;
        case Ea:
          a.$a = -1 / e;
      }
      1 > Math.abs(a.$a) && (a.$a *= 1e3);
      for (t = l = 1; t <= d; t++)
        for (m[t] = l, e = b.f[t].k; null != e; e = e.I)
          (n[l] = e.n.ea), (p[l] = e.n.ma * e.j * e.f.va), l++;
      m[d + 1] = l;
      ga(q, 1, b.head, 1, c);
      m = 0;
      for (t = 1; t <= c; t++)
        (u = b.n[t]), u.m != A && (m++, (q[c + m] = t), (K[m] = u.m));
      for (t = 1; t <= d; t++)
        (u = b.f[t]), u.m != A && (m++, (q[c + m] = c + t), (K[m] = u.m));
      a.valid = 1;
      b.valid = 0;
      a.U = b.U;
      b.U = null;
      q = a.g;
      K = a.i;
      m = a.Ca;
      n = a.Ba;
      p = a.De;
      t = a.Kd;
      for (l = 1; l <= q; l++) t[l] = 1;
      for (l = 1; l <= K; l++)
        for (f = m[l], e = m[l + 1]; f < e; f++) t[n[f]]++;
      for (l = p[1] = 1; l <= q; l++)
        t[l] > K && (t[l] = K), (p[l + 1] = p[l] + t[l]);
      a.Ce = new Int32Array(p[q + 1]);
      a.Ee = new Float64Array(p[q + 1]);
      q = a.g;
      K = a.i;
      m = a.head;
      n = a.m;
      ha(a.Kd, 1, 0, q);
      for (p = 1; p <= K; p++) n[p] != Na && ((t = m[q + p]), h(a, p, t));
      a.D = 0;
      a.hc = ja();
      a.Sf = a.$ = b.$;
      a.be = -1;
      a.Rb = 0;
      ha(r, 1, 0, c + d);
      for (t = 1; t <= d; t++) v[t] = 1;
    })(v, a);
    for (b.o >= mc && x("Objective scale factor = " + v.$a + ""); ; ) {
      if (0 == S) {
        sa = d(v);
        if (0 != sa)
          return (
            b.o >= Lb &&
              (x("Error: unable to factorize the basis matrix (" + sa + ")"),
              x("Sorry, basis recovery procedure not implemented yet")),
            (a.U = v.U),
            (v.U = null),
            (a.na = a.sa = Aa),
            (a.aa = 0),
            (a.$ = v.$),
            (a.some = 0),
            (sa = Sb)
          );
        S = v.valid = 1;
        M = ba = 0;
      }
      if (
        0 == M &&
        (p(v, v.Ha),
        (M = 1),
        0 == v.D &&
          (0 < W(v, b.Gb) ? (v.D = 1) : (X(v), (v.D = 2)),
          (ba = 0),
          u(v, b, 1)),
        ca(v, b.Gb))
      ) {
        b.o >= Lb &&
          x(
            "Warning: numerical instability (primal simplex, phase " +
              (1 == v.D ? "I" : "II") +
              ")"
          );
        S = v.D = 0;
        J = 5;
        continue;
      }
      1 != v.D || ka(v, b.Gb) || ((v.D = 2), X(v), (ba = 0), u(v, b, 1));
      0 == ba && (m(v), (ba = 1));
      switch (b.fd) {
        case oc:
          0 == v.Rb && q(v);
      }
      if (2147483647 > b.oc && v.$ - v.Sf >= b.oc) {
        if (1 != M || (2 == v.D && 1 != ba)) {
          1 != M && (M = 0);
          2 == v.D && 1 != ba && (ba = 0);
          continue;
        }
        u(v, b, 1);
        b.o >= Wb && x("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED");
        switch (v.D) {
          case 1:
            ea = Ad;
            X(v);
            m(v);
            break;
          case 2:
            ea = dc;
        }
        r(v, b.tb);
        fa = 0 == v.F ? dc : Ad;
        z(v, a, ea, fa, 0);
        return (sa = rg);
      }
      if (2147483647 > b.sb && 1e3 * la(v.hc) >= b.sb) {
        if (1 != M || (2 == v.D && 1 != ba)) {
          1 != M && (M = 0);
          2 == v.D && 1 != ba && (ba = 0);
          continue;
        }
        u(v, b, 1);
        b.o >= Wb && x("TIME LIMIT EXCEEDED; SEARCH TERMINATED");
        switch (v.D) {
          case 1:
            ea = Ad;
            X(v);
            m(v);
            break;
          case 2:
            ea = dc;
        }
        r(v, b.tb);
        fa = 0 == v.F ? dc : Ad;
        z(v, a, ea, fa, 0);
        return (sa = Qc);
      }
      u(v, b, 0);
      r(v, b.tb);
      if (0 == v.F) {
        if (1 != M || 1 != ba) {
          1 != M && (M = 0);
          1 != ba && (ba = 0);
          continue;
        }
        u(v, b, 1);
        switch (v.D) {
          case 1:
            b.o >= Wb && x("PROBLEM HAS NO FEASIBLE SOLUTION");
            ea = jc;
            X(v);
            m(v);
            r(v, b.tb);
            fa = 0 == v.F ? dc : Ad;
            break;
          case 2:
            b.o >= Wb && x("OPTIMAL SOLUTION FOUND"), (ea = fa = dc);
        }
        z(v, a, ea, fa, 0);
        return (sa = 0);
      }
      n(v);
      J && t(v);
      y(v, b.xe);
      var K = v.Qa[v.F],
        oa = V(v);
      if (
        Math.abs(K - oa) > 1e-5 * (1 + Math.abs(oa)) ||
        !((0 > K && 0 > oa) || (0 < K && 0 < oa))
      )
        if (
          (b.o >= mc && x("d1 = " + K + "; d2 = " + oa + ""), 1 != ba || !J)
        ) {
          1 != ba && (ba = 0);
          J = 5;
          continue;
        }
      v.Qa[v.F] =
        0 < K
          ? 0 < oa
            ? oa
            : 2.220446049250313e-16
          : 0 > oa
          ? oa
          : -2.220446049250313e-16;
      switch (b.ne) {
        case pc:
          E(v, 0);
          break;
        case qc:
          E(v, 0.3 * b.Gb);
      }
      if (0 == v.s) {
        if (1 != M || 1 != ba || !J) {
          1 != M && (M = 0);
          1 != ba && (ba = 0);
          J = 1;
          continue;
        }
        u(v, b, 1);
        switch (v.D) {
          case 1:
            b.o >= Lb && x("Error: unable to choose basic variable on phase I");
            a.U = v.U;
            v.U = null;
            a.na = a.sa = Aa;
            a.aa = 0;
            a.$ = v.$;
            a.some = 0;
            sa = Sb;
            break;
          case 2:
            b.o >= Wb && x("PROBLEM HAS UNBOUNDED SOLUTION"),
              z(v, a, dc, jc, v.head[v.g + v.F]),
              (sa = 0);
        }
        return sa;
      }
      if (
        0 < v.s &&
        ((K = v.Ua[v.s]),
        (oa = 1e-5 * (1 + 0.01 * v.rh)),
        Math.abs(K) < oa &&
          (b.o >= mc && x("piv = " + K + "; eps = " + oa + ""), !J))
      ) {
        J = 5;
        continue;
      }
      0 < v.s && ((K = v.jg), C(v, K), J && D(v, K), H(v, K));
      if (
        0 < v.s &&
        ((K = v.Ua[v.s]),
        (oa = v.mb[v.F]),
        Math.abs(K - oa) > 1e-8 * (1 + Math.abs(K)) ||
          !((0 < K && 0 < oa) || (0 > K && 0 > oa)))
      ) {
        b.o >= mc && x("piv1 = " + K + "; piv2 = " + oa + "");
        if (1 != S || !J) {
          1 != S && (S = 0);
          J = 5;
          continue;
        }
        0 == v.mb[v.F] && (v.Bc++, (v.Vb[v.Bc] = v.F));
        v.mb[v.F] = K;
      }
      R(v);
      M = 2;
      0 < v.s &&
        (O(v),
        (ba = 2),
        1 == v.D && ((K = v.head[v.s]), (v.Qa[v.F] -= v.u[K]), (v.u[K] = 0)));
      if (0 < v.s)
        switch (b.fd) {
          case oc:
            0 < v.Rb && Q(v);
        }
      0 < v.s &&
        ((sa = e(v, v.s, v.head[v.g + v.F])),
        (S = 0 == sa ? 2 : (v.valid = 0)));
      0 < v.s &&
        (k(v, v.F, v.head[v.g + v.F]),
        v.type[v.head[v.s]] != B && h(v, v.F, v.head[v.s]));
      F(v);
      v.$++;
      0 < J && J--;
    }
  }
  function Rb(a, b) {
    function c(a, b, c, d) {
      var e = a.g,
        f = a.Ca,
        g = a.Ba,
        h = a.Ia;
      a = a.head[b];
      if (a <= e) (e = 1), (c[1] = a), (d[1] = 1);
      else
        for (
          b = f[a - e],
            e = f[a - e + 1] - b,
            ga(c, 1, g, b, e),
            ga(d, 1, h, b, e),
            c = 1;
          c <= e;
          c++
        )
          d[c] = -d[c];
      return e;
    }
    function d(a) {
      var b = od(a.U, a.g, c, a);
      a.valid = 0 == b;
      return b;
    }
    function e(a, b, c) {
      var d = a.g,
        e;
      if (c <= d) {
        var f = Array(2);
        e = Array(2);
        f[1] = c;
        e[1] = 1;
        b = Ne(a.U, b, 1, f, 0, e);
      } else {
        var g = a.Ca,
          f = a.Ba,
          h = a.Ia;
        e = a.jb;
        var k;
        k = g[c - d];
        c = g[c - d + 1];
        g = 0;
        for (d = k; d < c; d++) e[++g] = -h[d];
        b = Ne(a.U, b, g, f, k - 1, e);
      }
      a.valid = 0 == b;
      return b;
    }
    function f(a, b, c) {
      var d = a.g,
        e = a.jb,
        f = a.jb,
        g = a.g,
        h = a.Ca,
        k = a.Ba,
        l = a.Ia,
        m = a.head,
        n,
        p,
        q;
      ga(f, 1, b, 1, g);
      for (b = 1; b <= g; b++)
        if (((q = c[b]), 0 != q))
          if (((n = m[b]), n <= g)) f[n] -= q;
          else
            for (p = h[n - g], n = h[n - g + 1]; p < n; p++)
              f[k[p]] += l[p] * q;
      wd(a.U, e);
      for (a = 1; a <= d; a++) c[a] += e[a];
    }
    function g(a, b, c) {
      var d = a.g,
        e = a.jb,
        f = a.jb,
        g = a.g,
        h = a.Ca,
        k = a.Ba,
        l = a.Ia,
        m = a.head,
        n,
        p,
        q,
        r;
      for (n = 1; n <= g; n++) {
        p = m[n];
        r = b[n];
        if (p <= g) r -= c[p];
        else
          for (q = h[p - g], p = h[p - g + 1]; q < p; q++) r += l[q] * c[k[q]];
        f[n] = r;
      }
      yd(a.U, e);
      for (a = 1; a <= d; a++) c[a] += e[a];
    }
    function h(a, b) {
      var c = a.c,
        d = a.d,
        e = a.m,
        f,
        g;
      f = a.head[a.g + b];
      switch (e[b]) {
        case G:
          g = c[f];
          break;
        case Ua:
          g = d[f];
          break;
        case Ra:
          g = 0;
          break;
        case Na:
          g = c[f];
      }
      return g;
    }
    function k(a, b) {
      var c = a.g,
        d = a.i,
        e = a.Ca,
        g = a.Ba,
        k = a.Ia,
        l = a.head,
        m = a.Gc,
        n,
        p,
        q,
        r;
      for (n = 1; n <= c; n++) m[n] = 0;
      for (n = 1; n <= d; n++)
        if (((p = l[c + n]), (r = h(a, n)), 0 != r))
          if (p <= c) m[p] -= r;
          else
            for (q = e[p - c], p = e[p - c + 1]; q < p; q++)
              m[g[q]] += r * k[q];
      ga(b, 1, m, 1, c);
      wd(a.U, b);
      f(a, m, b);
    }
    function l(a) {
      var b = a.i,
        c = a.Qa,
        d = a.Hc,
        e;
      e = a.g;
      var f = a.u,
        h = a.head,
        k = a.Gc,
        l;
      for (l = 1; l <= e; l++) k[l] = f[h[l]];
      ga(d, 1, k, 1, e);
      yd(a.U, d);
      g(a, k, d);
      for (e = 1; e <= b; e++) {
        f = a.g;
        l = a.u;
        k = h = void 0;
        h = a.head[f + e];
        k = l[h];
        if (h <= f) k -= d[h];
        else {
          l = a.Ca;
          for (
            var m = a.Ba,
              n = a.Ia,
              p = void 0,
              q = void 0,
              p = void 0,
              p = l[h - f],
              q = l[h - f + 1];
            p < q;
            p++
          )
            k += n[p] * d[m[p]];
        }
        c[e] = k;
      }
    }
    function p(a) {
      var b = a.g,
        c = a.i,
        d = a.head,
        e = a.hd,
        f = a.gamma;
      a.Rb = 1e3;
      ha(e, 1, 0, b + c);
      for (a = 1; a <= b; a++) (c = d[a]), (e[c] = 1), (f[a] = 1);
    }
    function m(a, b) {
      var c = a.g,
        d = a.type,
        e = a.c,
        f = a.d,
        g = a.head,
        h = a.Ha,
        k = a.gamma,
        l,
        m,
        n,
        p,
        q,
        r,
        t;
      q = p = n = 0;
      for (l = 1; l <= c; l++) {
        m = g[l];
        t = 0;
        if (d[m] == Sa || d[m] == I || d[m] == B)
          (r = b * (1 + P * Math.abs(e[m]))),
            h[l] < e[m] - r && (t = e[m] - h[l]);
        if (d[m] == Ta || d[m] == I || d[m] == B)
          (r = b * (1 + P * Math.abs(f[m]))),
            h[l] > f[m] + r && (t = f[m] - h[l]);
        0 != t &&
          ((m = k[l]),
          2.220446049250313e-16 > m && (m = 2.220446049250313e-16),
          (m = (t * t) / m),
          q < m && ((n = l), (p = t), (q = m)));
      }
      a.s = n;
      a.Se = p;
    }
    function q(a, b) {
      var c = a.g,
        d = a.s,
        e;
      for (e = 1; e <= c; e++) b[e] = 0;
      b[d] = 1;
      yd(a.U, ea);
    }
    function r(a, b) {
      var c = a.g,
        d = a.s,
        e = a.Hc,
        f;
      for (f = 1; f <= c; f++) e[f] = 0;
      e[d] = 1;
      g(a, e, b);
    }
    function n(a, b) {
      var c = a.g,
        d,
        e;
      e = 0;
      for (d = 1; d <= c; d++) 0 != b[d] && e++;
      if (0.2 <= e / c) {
        c = a.g;
        d = a.i;
        e = a.Ca;
        var f = a.Ba,
          g = a.Ia,
          h = a.head,
          k = a.m,
          l = a.Vb,
          m = a.mb,
          n,
          p,
          q,
          r,
          t;
        r = 0;
        for (n = 1; n <= d; n++)
          if (k[n] == Na) m[n] = 0;
          else {
            p = h[c + n];
            if (p <= c) t = -b[p];
            else
              for (q = e[p - c], p = e[p - c + 1], t = 0; q < p; q++)
                t += b[f[q]] * g[q];
            0 != t && (l[++r] = n);
            m[n] = t;
          }
        a.Bc = r;
      } else {
        f = a.g;
        c = a.i;
        g = a.lg;
        h = a.kg;
        k = a.mg;
        l = a.bind;
        m = a.m;
        d = a.Vb;
        e = a.mb;
        for (r = 1; r <= c; r++) e[r] = 0;
        for (n = 1; n <= f; n++)
          if (((p = b[n]), 0 != p))
            for (
              r = l[n] - f,
                1 <= r && m[r] != Na && (e[r] -= p),
                r = g[n],
                q = g[n + 1],
                t = r;
              t < q;
              t++
            )
              (r = l[f + h[t]] - f), 1 <= r && m[r] != Na && (e[r] += p * k[t]);
        f = 0;
        for (r = 1; r <= c; r++) 0 != e[r] && (d[++f] = r);
        a.Bc = f;
      }
    }
    function t(a, b) {
      var c = a.Bc,
        d = a.Vb,
        e = a.mb,
        f,
        g,
        h;
      g = 0;
      for (f = 1; f <= c; f++) (h = Math.abs(e[d[f]])), g < h && (g = h);
      a.uh = g;
      h = b * (1 + 0.01 * g);
      for (g = 0; g < c; )
        (f = d[c]), Math.abs(e[f]) < h ? c-- : (g++, (d[c] = d[g]), (d[g] = f));
      a.vh = g;
    }
    function y(a, b) {
      var c = a.m,
        d = a.Qa,
        e = a.Vb,
        f = a.mb,
        g = a.vh,
        h,
        k,
        l,
        m,
        n,
        p,
        q,
        r,
        t;
      p = 0 < a.Se ? 1 : -1;
      l = 0;
      r = s;
      n = 0;
      for (k = 1; k <= g; k++) {
        h = e[k];
        m = p * f[h];
        if (0 < m)
          if (c[h] == G || c[h] == Ra) q = (d[h] + b) / m;
          else continue;
        else if (c[h] == Ua || c[h] == Ra) q = (d[h] - b) / m;
        else continue;
        0 > q && (q = 0);
        if (r > q || (r == q && n < Math.abs(m)))
          (l = h), (r = q), (n = Math.abs(m));
      }
      if (0 != b && 0 != l && 0 != r)
        for (t = r, l = 0, r = s, n = 0, k = 1; k <= g; k++) {
          h = e[k];
          m = p * f[h];
          if (0 < m)
            if (c[h] == G || c[h] == Ra) q = d[h] / m;
            else continue;
          else if (c[h] == Ua || c[h] == Ra) q = d[h] / m;
          else continue;
          0 > q && (q = 0);
          q <= t && n < Math.abs(m) && ((l = h), (r = q), (n = Math.abs(m)));
        }
      a.F = l;
      a.bh = p * r;
    }
    function E(a) {
      var b = a.g,
        c = a.yb,
        d = a.Ua,
        e = a.Ua,
        f,
        g;
      g = a.head[b + a.F];
      for (f = 1; f <= b; f++) e[f] = 0;
      if (g <= b) e[g] = -1;
      else {
        var h = a.Ca;
        f = a.Ba;
        var k = a.Ia,
          l;
        l = h[g - b];
        for (g = h[g - b + 1]; l < g; l++) e[f[l]] = k[l];
      }
      wd(a.U, d);
      e = 0;
      for (f = 1; f <= b; f++) 0 != d[f] && (c[++e] = f);
      a.Tb = e;
    }
    function C(a) {
      var b = a.g,
        c = a.yb,
        d = a.Ua,
        e = a.Hc,
        g,
        h;
      h = a.head[b + a.F];
      for (g = 1; g <= b; g++) e[g] = 0;
      if (h <= b) e[h] = -1;
      else {
        var k = a.Ca;
        g = a.Ba;
        var l = a.Ia,
          m;
        m = k[h - b];
        for (h = k[h - b + 1]; m < h; m++) e[g[m]] = l[m];
      }
      f(a, e, d);
      e = 0;
      for (g = 1; g <= b; g++) 0 != d[g] && (c[++e] = g);
      a.Tb = e;
    }
    function D(a) {
      var b = a.Qa,
        c = a.Bc,
        d = a.Vb,
        e = a.mb,
        f = a.F;
      a = a.bh;
      var g, h;
      b[f] = a;
      if (0 != a)
        for (h = 1; h <= c; h++) (g = d[h]), g != f && (b[g] -= e[g] * a);
    }
    function H(a) {
      var b = a.Ha,
        c = a.s,
        d = a.F,
        e = a.Tb,
        f = a.yb,
        g = a.Ua,
        k;
      k = a.Se / g[c];
      b[c] = h(a, d) + k;
      if (0 != k)
        for (d = 1; d <= e; d++) (a = f[d]), a != c && (b[a] += g[a] * k);
    }
    function R(a) {
      var b = a.g,
        c = a.type,
        d = a.head,
        e = a.hd,
        f = a.gamma,
        g = a.s,
        h = a.Bc,
        k = a.Vb,
        l = a.mb,
        m = a.F,
        n = a.Tb,
        p = a.yb,
        q = a.Ua,
        r = a.Hc,
        t,
        u,
        v,
        y,
        z,
        C;
      a.Rb--;
      z = C = e[d[g]] ? 1 : 0;
      for (t = 1; t <= b; t++) r[t] = 0;
      for (y = 1; y <= h; y++)
        if (((u = k[y]), (v = d[b + u]), e[v]))
          if (((u = l[u]), (z += u * u), v <= b)) r[v] += u;
          else {
            var D = a.Ca;
            t = a.Ba;
            var F = a.Ia,
              E;
            E = D[v - b];
            for (v = D[v - b + 1]; E < v; E++) r[t[E]] -= u * F[E];
          }
      wd(a.U, r);
      a = q[g];
      for (y = 1; y <= n; y++)
        (t = p[y]),
          (v = d[t]),
          t != g &&
            c[d[t]] != Ka &&
            ((u = q[t] / a),
            (h = f[t] + u * u * z + 2 * u * r[t]),
            (u = (e[v] ? 1 : 0) + C * u * u),
            (f[t] = h >= u ? h : u),
            2.220446049250313e-16 > f[t] && (f[t] = 2.220446049250313e-16));
      c[d[b + m]] == Ka
        ? (f[g] = 1)
        : ((f[g] = z / (a * a)),
          2.220446049250313e-16 > f[g] && (f[g] = 2.220446049250313e-16));
      v = d[g];
      if (c[v] == B && e[v])
        for (e[v] = 0, y = 1; y <= n; y++) {
          t = p[y];
          if (t == g) {
            if (c[d[b + m]] == Ka) continue;
            u = 1 / q[g];
          } else {
            if (c[d[t]] == Ka) continue;
            u = q[t] / q[g];
          }
          f[t] -= u * u;
          2.220446049250313e-16 > f[t] && (f[t] = 2.220446049250313e-16);
        }
    }
    function V(a) {
      var b = a.g,
        c = a.type,
        d = a.head,
        e = a.bind,
        f = a.m,
        g = a.s,
        h = a.Se;
      a = a.F;
      var k;
      k = d[g];
      d[g] = d[b + a];
      d[b + a] = k;
      e[d[g]] = g;
      e[d[b + a]] = b + a;
      f[a] = c[k] == B ? Na : 0 < h ? G : Ua;
    }
    function O(a, b) {
      var c = a.g,
        d = a.i,
        e = a.ac,
        f = a.head,
        g = a.Qa,
        h,
        k;
      for (h = 1; h <= d; h++)
        if (
          ((k = f[c + h]),
          (g[h] < -b && (e[k] == Sa || e[k] == Ka)) ||
            (g[h] > +b && (e[k] == Ta || e[k] == Ka)))
        )
          return 1;
      return 0;
    }
    function Q(a) {
      var b = a.g,
        c = a.i,
        d = a.type,
        e = a.c,
        f = a.d,
        g = a.ac,
        h = a.head,
        k = a.m;
      a = a.Qa;
      var l;
      for (l = 1; l <= b + c; l++)
        switch (g[l]) {
          case Ka:
            d[l] = I;
            e[l] = -1e3;
            f[l] = 1e3;
            break;
          case Sa:
            d[l] = I;
            e[l] = 0;
            f[l] = 1;
            break;
          case Ta:
            d[l] = I;
            e[l] = -1;
            f[l] = 0;
            break;
          case I:
          case B:
            (d[l] = B), (e[l] = f[l] = 0);
        }
      for (e = 1; e <= c; e++)
        (l = h[b + e]), (k[e] = d[l] == B ? Na : 0 <= a[e] ? G : Ua);
    }
    function F(a) {
      var b = a.g,
        c = a.i,
        d = a.type,
        e = a.c,
        f = a.d,
        g = a.bd,
        h = a.cd,
        k = a.head,
        l = a.m,
        m = a.Qa;
      ga(d, 1, a.ac, 1, b + c);
      ga(e, 1, g, 1, b + c);
      ga(f, 1, h, 1, b + c);
      for (a = 1; a <= c; a++)
        switch (((g = k[b + a]), d[g])) {
          case Ka:
            l[a] = Ra;
            break;
          case Sa:
            l[a] = G;
            break;
          case Ta:
            l[a] = Ua;
            break;
          case I:
            l[a] =
              2.220446049250313e-16 <= m[a]
                ? G
                : -2.220446049250313e-16 >= m[a]
                ? Ua
                : Math.abs(e[g]) <= Math.abs(f[g])
                ? G
                : Ua;
            break;
          case B:
            l[a] = Na;
        }
    }
    function W(a, b) {
      var c = a.i,
        d = a.m,
        e = a.Qa,
        f;
      for (f = 1; f <= c; f++)
        if (
          (e[f] < -b && (d[f] == G || d[f] == Ra)) ||
          (e[f] > +b && (d[f] == Ua || d[f] == Ra))
        )
          return 1;
      return 0;
    }
    function X(a) {
      var b = a.g,
        c = a.i,
        d = a.eb,
        e = a.head,
        f = a.Ha,
        g,
        k,
        l;
      l = d[0];
      for (g = 1; g <= b; g++) (k = e[g]), k > b && (l += d[k - b] * f[g]);
      for (f = 1; f <= c; f++)
        (k = e[b + f]), k > b && (l += d[k - b] * h(a, f));
      return l;
    }
    function ca(a, b, c) {
      var d = a.g,
        e = a.i,
        f = a.u,
        g = a.ac,
        k = a.head,
        l = a.m,
        m = a.D,
        n = a.Ha,
        p = a.Qa;
      if (
        !(
          b.o < fc ||
          (0 < b.fb && 1e3 * la(a.hc) < b.fb) ||
          a.$ == a.be ||
          (!c && 0 != a.$ % b.bc)
        )
      ) {
        b = 0;
        if (1 == m) {
          for (l = 1; l <= d; l++) b -= f[k[l]] * n[l];
          for (n = 1; n <= e; n++) b -= f[k[d + n]] * h(a, n);
        } else
          for (n = 1; n <= e; n++)
            0 > p[n] && (l[n] == G || l[n] == Ra) && (b -= p[n]),
              0 < p[n] && (l[n] == Ua || l[n] == Ra) && (b += p[n]);
        e = 0;
        for (l = 1; l <= d; l++) g[k[l]] == B && e++;
        1 == a.D
          ? x(" " + a.$ + ":  infeas = " + b + " (" + e + ")")
          : x(
              "|" + a.$ + ": obj = " + X(a) + "  infeas = " + b + " (" + e + ")"
            );
        a.be = a.$;
      }
    }
    function ka(a, b, c, d, e) {
      var f = a.g,
        g = a.i,
        h = a.$a,
        k = a.head,
        l = a.m,
        m = a.Ha,
        n = a.Qa;
      b.valid = 1;
      a.valid = 0;
      b.U = a.U;
      a.U = null;
      ga(b.head, 1, k, 1, f);
      b.na = c;
      b.sa = d;
      b.aa = X(a);
      b.$ = a.$;
      b.some = e;
      for (a = 1; a <= f; a++)
        (c = k[a]),
          c <= f
            ? ((c = b.n[c]), (c.m = A), (c.bind = a), (c.r = m[a] / c.ma))
            : ((c = b.f[c - f]), (c.m = A), (c.bind = a), (c.r = m[a] * c.va)),
          (c.J = 0);
      for (m = 1; m <= g; m++)
        if (((c = k[f + m]), c <= f)) {
          c = b.n[c];
          c.m = l[m];
          c.bind = 0;
          switch (l[m]) {
            case G:
              c.r = c.c;
              break;
            case Ua:
              c.r = c.d;
              break;
            case Ra:
              c.r = 0;
              break;
            case Na:
              c.r = c.c;
          }
          c.J = (n[m] * c.ma) / h;
        } else {
          c = b.f[c - f];
          c.m = l[m];
          c.bind = 0;
          switch (l[m]) {
            case G:
              c.r = c.c;
              break;
            case Ua:
              c.r = c.d;
              break;
            case Ra:
              c.r = 0;
              break;
            case Na:
              c.r = c.c;
          }
          c.J = n[m] / c.va / h;
        }
    }
    var P = 0.1,
      u,
      z = 2,
      L = 0,
      v = 0,
      S = 0,
      M,
      ba,
      J;
    u = (function (a) {
      var b = a.g,
        c = a.i;
      a = a.L;
      var d = {};
      d.g = b;
      d.i = c;
      d.type = new Int8Array(1 + b + c);
      d.c = new Float64Array(1 + b + c);
      d.d = new Float64Array(1 + b + c);
      d.u = new Float64Array(1 + b + c);
      d.ac = new Int8Array(1 + b + c);
      d.bd = new Float64Array(1 + b + c);
      d.cd = new Float64Array(1 + b + c);
      d.eb = new Float64Array(1 + c);
      d.Ca = new Int32Array(1 + c + 1);
      d.Ba = new Int32Array(1 + a);
      d.Ia = new Float64Array(1 + a);
      d.lg = new Int32Array(1 + b + 1);
      d.kg = new Int32Array(1 + a);
      d.mg = new Float64Array(1 + a);
      d.head = new Int32Array(1 + b + c);
      d.bind = new Int32Array(1 + b + c);
      d.m = new Int8Array(1 + c);
      d.Ha = new Float64Array(1 + b);
      d.Qa = new Float64Array(1 + c);
      d.hd = new Int8Array(1 + b + c);
      d.gamma = new Float64Array(1 + b);
      d.Vb = new Int32Array(1 + c);
      d.mb = new Float64Array(1 + c);
      d.yb = new Int32Array(1 + b);
      d.Ua = new Float64Array(1 + b);
      d.jb = new Float64Array(1 + b);
      d.Gc = new Float64Array(1 + b);
      d.Hc = new Float64Array(1 + b);
      d.jg = new Float64Array(1 + b);
      return d;
    })(a);
    (function (a, b) {
      var c = a.g,
        d = a.i,
        e = a.type,
        f = a.c,
        g = a.d,
        h = a.u,
        k = a.ac,
        l = a.bd,
        m = a.cd,
        n = a.eb,
        p = a.Ca,
        q = a.Ba,
        r = a.Ia,
        t = a.lg,
        u = a.kg,
        v = a.mg,
        y = a.head,
        z = a.bind,
        C = a.m,
        D = a.hd,
        F = a.gamma,
        E,
        H;
      for (E = 1; E <= c; E++)
        (H = b.n[E]),
          (e[E] = H.type),
          (f[E] = H.c * H.ma),
          (g[E] = H.d * H.ma),
          (h[E] = 0);
      for (E = 1; E <= d; E++)
        (H = b.f[E]),
          (e[c + E] = H.type),
          (f[c + E] = H.c / H.va),
          (g[c + E] = H.d / H.va),
          (h[c + E] = H.u * H.va);
      ga(k, 1, e, 1, c + d);
      ga(l, 1, f, 1, c + d);
      ga(m, 1, g, 1, c + d);
      n[0] = b.ha;
      ga(n, 1, h, c + 1, d);
      e = 0;
      for (E = 1; E <= d; E++) e < Math.abs(n[E]) && (e = Math.abs(n[E]));
      0 == e && (e = 1);
      switch (b.dir) {
        case za:
          a.$a = 1 / e;
          break;
        case Ea:
          a.$a = -1 / e;
      }
      1 > Math.abs(a.$a) && (a.$a *= 1e3);
      for (E = 1; E <= d; E++) h[c + E] *= a.$a;
      for (E = h = 1; E <= d; E++)
        for (p[E] = h, n = b.f[E].k; null != n; n = n.I)
          (q[h] = n.n.ea), (r[h] = n.n.ma * n.j * n.f.va), h++;
      p[d + 1] = h;
      for (E = h = 1; E <= c; E++)
        for (t[E] = h, n = b.n[E].k; null != n; n = n.B)
          (u[h] = n.f.C), (v[h] = n.n.ma * n.j * n.f.va), h++;
      t[c + 1] = h;
      ga(y, 1, b.head, 1, c);
      p = 0;
      for (E = 1; E <= c; E++)
        (H = b.n[E]), H.m != A && (p++, (y[c + p] = E), (C[p] = H.m));
      for (E = 1; E <= d; E++)
        (H = b.f[E]), H.m != A && (p++, (y[c + p] = c + E), (C[p] = H.m));
      for (p = 1; p <= c + d; p++) z[y[p]] = p;
      a.valid = 1;
      b.valid = 0;
      a.U = b.U;
      b.U = null;
      a.D = 0;
      a.hc = ja();
      a.Sf = a.$ = b.$;
      a.be = -1;
      a.Rb = 0;
      ha(D, 1, 0, c + d);
      for (E = 1; E <= c; E++) F[E] = 1;
    })(u, a);
    for (b.o >= mc && x("Objective scale factor = " + u.$a + ""); ; ) {
      if (0 == z) {
        J = d(u);
        if (0 != J)
          return (
            b.o >= Lb &&
              (x("Error: unable to factorize the basis matrix (" + J + ")"),
              x("Sorry, basis recovery procedure not implemented yet")),
            (a.U = u.U),
            (u.U = null),
            (a.na = a.sa = Aa),
            (a.aa = 0),
            (a.$ = u.$),
            (a.some = 0),
            (J = Sb)
          );
        z = u.valid = 1;
        L = v = 0;
      }
      if (
        0 == v &&
        (l(u),
        (v = 1),
        0 == u.D &&
          (0 != O(u, 0.9 * b.tb) ? ((u.D = 1), Q(u)) : ((u.D = 2), F(u)),
          (L = u.Rb = 0)),
        0 != W(u, b.tb))
      ) {
        b.o >= Lb &&
          x(
            "Warning: numerical instability (dual simplex, phase " +
              (1 == u.D ? "I" : "II") +
              ")"
          );
        if (b.cb == Qb) return ka(u, a, Aa, Aa, 0), (J = Sb);
        z = u.D = 0;
        S = 5;
        continue;
      }
      1 == u.D &&
        0 == O(u, b.tb) &&
        (ca(u, b, 1),
        (u.D = 2),
        1 != v && (l(u), (v = 1)),
        F(u),
        (L = u.Rb = 0));
      0 == L && (k(u, u.Ha), 2 == u.D && (u.Ha[0] = X(u)), (L = 1));
      switch (b.fd) {
        case oc:
          0 == u.Rb && p(u);
      }
      if (2 == u.D && 0 > u.$a && b.hf > -s && u.Ha[0] <= b.hf) {
        if (1 != L || 1 != v) {
          1 != L && (L = 0);
          1 != v && (v = 0);
          continue;
        }
        ca(u, b, 1);
        b.o >= Wb && x("OBJECTIVE LOWER LIMIT REACHED; SEARCH TERMINATED");
        ka(u, a, Ad, dc, 0);
        return (J = Vf);
      }
      if (2 == u.D && 0 < u.$a && b.jf < +s && u.Ha[0] >= b.jf) {
        if (1 != L || 1 != v) {
          1 != L && (L = 0);
          1 != v && (v = 0);
          continue;
        }
        ca(u, b, 1);
        b.o >= Wb && x("OBJECTIVE UPPER LIMIT REACHED; SEARCH TERMINATED");
        ka(u, a, Ad, dc, 0);
        return (J = Wf);
      }
      if (2147483647 > b.oc && u.$ - u.Sf >= b.oc) {
        if ((2 == u.D && 1 != L) || 1 != v) {
          2 == u.D && 1 != L && (L = 0);
          1 != v && (v = 0);
          continue;
        }
        ca(u, b, 1);
        b.o >= Wb && x("ITERATION LIMIT EXCEEDED; SEARCH TERMINATED");
        switch (u.D) {
          case 1:
            ba = Ad;
            F(u);
            k(u, u.Ha);
            break;
          case 2:
            ba = dc;
        }
        ka(u, a, Ad, ba, 0);
        return (J = rg);
      }
      if (2147483647 > b.sb && 1e3 * la(u.hc) >= b.sb) {
        if ((2 == u.D && 1 != L) || 1 != v) {
          2 == u.D && 1 != L && (L = 0);
          1 != v && (v = 0);
          continue;
        }
        ca(u, b, 1);
        b.o >= Wb && x("TIME LIMIT EXCEEDED; SEARCH TERMINATED");
        switch (u.D) {
          case 1:
            ba = Ad;
            F(u);
            k(u, u.Ha);
            break;
          case 2:
            ba = dc;
        }
        ka(u, a, Ad, ba, 0);
        return (J = Qc);
      }
      ca(u, b, 0);
      m(u, b.Gb);
      if (0 == u.s) {
        if (1 != L || 1 != v) {
          1 != L && (L = 0);
          1 != v && (v = 0);
          continue;
        }
        ca(u, b, 1);
        switch (u.D) {
          case 1:
            b.o >= Wb && x("PROBLEM HAS NO DUAL FEASIBLE SOLUTION");
            F(u);
            k(u, u.Ha);
            M = Ad;
            ba = jc;
            break;
          case 2:
            b.o >= Wb && x("OPTIMAL SOLUTION FOUND"), (M = ba = dc);
        }
        ka(u, a, M, ba, 0);
        return (J = 0);
      }
      var ea = u.jg;
      q(u, ea);
      S && r(u, ea);
      n(u, ea);
      t(u, b.Gb);
      switch (b.ne) {
        case pc:
          y(u, 0);
          break;
        case qc:
          y(u, 0.3 * b.tb);
      }
      if (0 == u.F) {
        if (1 != L || 1 != v || !S) {
          1 != L && (L = 0);
          1 != v && (v = 0);
          S = 1;
          continue;
        }
        ca(u, b, 1);
        switch (u.D) {
          case 1:
            b.o >= Lb && x("Error: unable to choose basic variable on phase I");
            a.U = u.U;
            u.U = null;
            a.na = a.sa = Aa;
            a.aa = 0;
            a.$ = u.$;
            a.some = 0;
            J = Sb;
            break;
          case 2:
            b.o >= Wb && x("PROBLEM HAS NO FEASIBLE SOLUTION"),
              ka(u, a, jc, dc, u.head[u.s]),
              (J = 0);
        }
        return J;
      }
      var fa = u.mb[u.F],
        sa = 1e-5 * (1 + 0.01 * u.uh);
      if (
        Math.abs(fa) < sa &&
        (b.o >= mc && x("piv = " + fa + "; eps = " + sa + ""), !S)
      ) {
        S = 5;
        continue;
      }
      E(u);
      S && C(u);
      fa = u.Ua[u.s];
      sa = u.mb[u.F];
      if (
        Math.abs(fa - sa) > 1e-8 * (1 + Math.abs(fa)) ||
        !((0 < fa && 0 < sa) || (0 > fa && 0 > sa))
      ) {
        b.o >= mc && x("piv1 = " + fa + "; piv2 = " + sa + "");
        if (1 != z || !S) {
          1 != z && (z = 0);
          S = 5;
          continue;
        }
        0 == u.Ua[u.s] && (u.Tb++, (u.yb[u.Tb] = u.s));
        u.Ua[u.s] = sa;
      }
      H(u);
      2 == u.D && (u.Ha[0] += (u.Qa[u.F] / u.$a) * (u.Se / u.Ua[u.s]));
      L = 2;
      D(u);
      v = 2;
      switch (b.fd) {
        case oc:
          0 < u.Rb && R(u);
      }
      J = e(u, u.s, u.head[u.g + u.F]);
      z = 0 == J ? 2 : (u.valid = 0);
      V(u);
      u.$++;
      0 < S && S--;
    }
  }
})((typeof exports === "object" && exports) || this);
