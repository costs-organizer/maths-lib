const {
  glp_create_prob,
  glp_set_obj_dir,
  GLP_MIN,
  glp_add_cols,
  glp_add_rows,
  glp_set_col_bnds,
  GLP_LO,
  glp_set_obj_coef,
  glp_set_col_kind,
  glp_set_row_bnds,
  GLP_FX,
  glp_set_mat_row,
  SMCP,
  GLP_ON,
  glp_scale_prob,
  GLP_SF_AUTO,
  glp_simplex,
  glp_get_col_prim,
  glp_get_obj_val,
  GLP_OPT,
  GLP_FEAS,
  GLP_INFEAS,
  GLP_NOFEAS,
  GLP_UNBND,
  GLP_UNDEF,
  glp_write_lp,
  glp_get_status,
} = require("./glpk");

function vectorCopy(a) {
  return new Float64Array(a);
}

const EPS = 2.2205e-16;

module.exports.mat = (elems, rowwise) => {
  var k;
  var concatWithNumbers = false;
  var elemtypes = new Array(elems.length);
  for (k = 0; k < elems.length; k++) {
    elemtypes[k] = type(elems[k]);
    if (elemtypes[k] == "number") concatWithNumbers = true;
  }
  if (typeof rowwise == "undefined") {
    if (type(elems) == "vector") return new Float64Array(elems);
    var rowwise = true;
    for (k = 0; k < elems.length; k++) {
      if (!Array.isArray(elems[k]) || elemtypes[k] == "vector") {
        rowwise = false;
        if (elemtypes[k] == "string") return elems;
      }
    }
  }
  if (elems.length == 0) {
    return [];
  }
  var m = 0;
  var n = 0;
  var i;
  var j;
  if (rowwise) {
    var res = new Array();
    for (k = 0; k < elems.length; k++) {
      switch (elemtypes[k]) {
        case "matrix":
          res.push(elems[k].val);
          m += elems[k].m;
          n = elems[k].n;
          break;
        case "vector":
          if (concatWithNumbers) {
            for (var l = 0; l < elems[k].length; l++) res.push(elems[k][l]);
            n = 1;
            m += elems[k].length;
          } else {
            res.push(elems[k]);
            m += 1;
            n = elems[k].length;
          }
          break;
        case "number":
          res.push(elems[k]);
          m += 1;
          n = 1;
          break;
        case "spvector":
          return spmat(elems);
        default:
          return elems;
      }
    }
    if (n == 1) {
      var M = new Float64Array(res);
      return M;
    }
    var M = new Matrix(m, n);
    var p = 0;
    for (k = 0; k < res.length; k++) {
      if (res[k].buffer) {
        M.val.set(res[k], p);
        p += res[k].length;
      } else {
        for (j = 0; j < res[k].length; j++) M.val[p + j] = res[k][j];
        p += res[k].length;
      }
    }
    return M;
  } else {
    m = size(elems[0], 1);
    for (k = 0; k < elems.length; k++) {
      if (elemtypes[k] == "matrix") n += elems[k].n;
      else n++;
      if (size(elems[k], 1) != m) return "undefined";
    }
    var res = new Matrix(m, n);
    var c;
    for (i = 0; i < m; i++) {
      c = 0;
      for (k = 0; k < elems.length; k++) {
        switch (elemtypes[k]) {
          case "matrix":
            for (j = 0; j < elems[k].n; j++) {
              res.val[i * n + j + c] = elems[k].val[i * elems[k].n + j];
            }
            c += elems[k].n;
            break;
          case "vector":
            res.val[i * n + c] = elems[k][i];
            c++;
            break;
          case "number":
            res.val[i * n + c] = elems[k];
            c++;
            break;
          default:
            break;
        }
      }
    }
    return res;
  }
};

const setVectorScalar = (A, rowsrange, B) => {
  var i;
  for (i = 0; i < rowsrange.length; i++) A[rowsrange[i]] = B;
};
const setVectorVector = (A, rowsrange, B) => {
  var i;
  for (i = 0; i < rowsrange.length; i++) A[rowsrange[i]] = B[i];
};
const setMatrixScalar = (A, rowsrange, colsrange, B) => {
  var i;
  var j;
  var m = rowsrange.length;
  var n = colsrange.length;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) A.val[rowsrange[i] * A.n + colsrange[j]] = B;
};
const setMatrixMatrix = (A, rowsrange, colsrange, B) => {
  var i;
  var j;
  var m = rowsrange.length;
  var n = colsrange.length;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A.val[rowsrange[i] * A.n + colsrange[j]] = B.val[i * B.n + j];
};
const setMatrixColVector = (A, rowsrange, col, B) => {
  var i;
  var m = rowsrange.length;
  for (i = 0; i < m; i++) A.val[rowsrange[i] * A.n + col] = B[i];
};
const setMatrixRowVector = (A, row, colsrange, B) => {
  var j;
  var n = colsrange.length;
  for (j = 0; j < n; j++) A.val[row * A.n + colsrange[j]] = B[j];
};
const setRows = (A, rowsrange, B) => {
  var i;
  var j;
  var m = rowsrange.length;
  var rA;
  switch (type(B)) {
    case "vector":
      for (i = 0; i < m; i++) {
        rA = rowsrange[i] * A.n;
        for (j = 0; j < B.length; j++) A.val[rA + j] = B[j];
      }
      break;
    case "matrix":
      var rB = 0;
      for (i = 0; i < m; i++) {
        rA = rowsrange[i] * A.n;
        for (j = 0; j < B.n; j++) A.val[rA + j] = B.val[rB + j];
        rB += B.n;
      }
      break;
    default:
      for (i = 0; i < m; i++) {
        rA = rowsrange[i] * A.n;
        for (j = 0; j < A.n; j++) A.val[rA + j] = B;
      }
      break;
  }
};
const setCols = (A, colsrange, B) => {
  var i;
  var m = A.m;
  var n = colsrange.length;
  var r = 0;
  switch (type(B)) {
    case "vector":
      for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) A.val[r + colsrange[j]] = B[i];
        r += A.n;
      }
      break;
    case "matrix":
      for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) A.val[r + colsrange[j]] = B.val[i * B.n + j];
        r += A.n;
      }
      break;
    default:
      for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) A.val[r + colsrange[j]] = B;
        r += A.n;
      }
      break;
  }
};

module.exports.set = (A, rowsrange, colsrange, B) => {
  var i;
  var j;
  var k;
  var l;
  var n;
  var typerows = typeof rowsrange;
  var typecols = typeof colsrange;
  if (arguments.length == 1) return undefined;
  var typeA = type(A);
  if (typeA == "vector") {
    B = colsrange;
    if (typerows == "number") {
      A[rowsrange] = B;
      return B;
    } else if (rowsrange.length == 0) rowsrange = range(A.length);
    if (size(B, 1) == 1) {
      setVectorScalar(A, rowsrange, B);
    } else {
      setVectorVector(A, rowsrange, B);
    }
    return B;
  } else if (typeA == "matrix") {
    if (typerows == "number") rowsrange = [rowsrange];
    if (typecols == "number") colsrange = [colsrange];
    if (rowsrange.length == 1 && colsrange.length == 1) {
      A.val[rowsrange[0] * A.n + colsrange[0]] = B;
      return B;
    }
    if (rowsrange.length == 0) {
      setCols(A, colsrange, B);
      return B;
    }
    if (colsrange.length == 0) {
      setRows(A, rowsrange, B);
      return B;
    }
    var sB = size(B);
    var tB = type(B);
    if (sB[0] == 1 && sB[1] == 1) {
      if (tB == "number") setMatrixScalar(A, rowsrange, colsrange, B);
      else if (tB == "vector") setMatrixScalar(A, rowsrange, colsrange, B[0]);
      else setMatrixScalar(A, rowsrange, colsrange, B.val[0]);
    } else {
      if (colsrange.length == 1)
        setMatrixColVector(A, rowsrange, colsrange[0], B);
      else if (rowsrange.length == 1) {
        if (tB == "vector") setMatrixRowVector(A, rowsrange[0], colsrange, B);
        else setMatrixRowVector(A, rowsrange[0], colsrange, B.val);
      } else setMatrixMatrix(A, rowsrange, colsrange, B);
    }
    return B;
  } else if (typeA == "ComplexVector") {
    B = colsrange;
    if (typerows == "number") {
      A.set(rowsrange, B);
      return B;
    } else if (rowsrange.length == 0) rowsrange = range(A.length);
    if (size(B, 1) == 1) {
      A.setVectorScalar(rowsrange, B);
    } else {
      A.setVectorVector(rowsrange, B);
    }
    return B;
  }
};

module.exports.isZero = (x) => {
  return Math.abs(x) < EPS;
};

module.exports.swaprows = (A, i, j) => {
  if (i != j) {
    var ri = i * A.n;
    var rj = j * A.n;
    var tmp = vectorCopy(A.val.subarray(ri, ri + A.n));
    A.val.set(vectorCopy(A.val.subarray(rj, rj + A.n)), ri);
    A.val.set(tmp, rj);
  }
};

module.exports.ones = (rows, cols) => {
  if (arguments.length == 1 || cols == 1) {
    var v = new Float64Array(rows);
    for (var i = 0; i < rows; i++) v[i] = 1;
    return v;
  } else {
    var M = new Matrix(rows, cols);
    const mn = rows * cols;
    for (var i = 0; i < mn; i++) {
      M.val[i] = 1;
    }
    return M;
  }
};
module.exports.zeros = (rows, cols) => {
  if (arguments.length == 1 || cols == 1) {
    return new Float64Array(rows);
  } else {
    return new Matrix(rows, cols);
  }
};
module.exports.eye = (m, n) => {
  if (typeof n == "undefined") var n = m;
  if (m == 1 && n == 1) return 1;
  var I = zeros(m, n);
  const e = m < n ? m : n;
  for (var i = 0; i < e; i++) {
    I.val[i * (n + 1)] = 1;
  }
  return I;
};

module.exports.glp = (
  c,
  A,
  b,
  Aeq,
  beq,
  lb,
  ub,
  integer_variables,
  verbose
) => {
  var prob = glp_create_prob();
  glp_set_obj_dir(prob, GLP_MIN);
  if (typeof Aeq == "undefined") var Aeq = [];
  glp_add_cols(prob, c.length);
  if (A.length + Aeq.length > 0) glp_add_rows(prob, A.length + Aeq.length);
  var i;
  var j;
  var indexes;
  var values;
  var n = c.length;
  if (lb) {
    var lbdense = vectorCopy(lb);
    for (i = 0; i < lbdense.length; i++) {
      if (!isFinite(lbdense[i])) lbdense[i] = NaN;
    }
  } else var lbdense = [];
  if (ub) {
    var ubdense = vectorCopy(ub);
    for (i = 0; i < ubdense.length; i++) {
      if (!isFinite(ubdense[i])) lbdense[i] = NaN;
    }
  } else var ubdense = [];
  for (i = 0; i < c.length; i++) {
    var lbi = NaN;
    var ubi = NaN;
    if (lbdense.length > 0) lbi = lbdense[i];
    if (ubdense.length > 0) ubi = ubdense[i];
    if (!isNaN(lbi) && !isNaN(ubi))
      glp_set_col_bnds(prob, i + 1, GLP_DB, lbi, ubi);
    else if (!isNaN(lbi)) glp_set_col_bnds(prob, i + 1, GLP_LO, lbi);
    else if (!isNaN(ubi)) glp_set_col_bnds(prob, i + 1, GLP_UP, 0, ubi);
    else glp_set_col_bnds(prob, i + 1, GLP_FR);
    glp_set_obj_coef(prob, i + 1, c[i]);
  }
  if (integer_variables) {
    for (i = 0; i < integer_variables.length; i++)
      glp_set_col_kind(prob, integer_variables[i] + 1, GLP_IV);
  }
  if (A.length == 1 && typeof b == "number") b = [b];
  for (i = 0; i < A.length; i++) {
    glp_set_row_bnds(prob, i + 1, GLP_UP, 0, b[i]);
    indexes = new Array();
    values = new Array();
    indexes.push(0);
    values.push(0);
    for (j = 0; j < n; j++) {
      if (!isZero(A.val[i * n + j])) {
        indexes.push(j + 1);
        values.push(A.val[i * n + j]);
      }
    }
    glp_set_mat_row(prob, i + 1, indexes.length - 1, indexes, values);
  }
  if (Aeq.length == 1 && typeof beq == "number") beq = [beq];
  for (i = 0; i < Aeq.length; i++) {
    glp_set_row_bnds(prob, A.length + i + 1, GLP_FX, beq[i]);
    indexes = new Array();
    values = new Array();
    indexes.push(0);
    values.push(0);
    for (j = 0; j < n; j++) {
      if (!isZero(Aeq.val[i * n + j])) {
        indexes.push(j + 1);
        values.push(Aeq.val[i * n + j]);
      }
    }
    glp_set_mat_row(
      prob,
      A.length + i + 1,
      indexes.length - 1,
      indexes,
      values
    );
  }
  var rc;
  if (integer_variables && integer_variables.length > 0) {
    var iocp = new IOCP({ presolve: GLP_ON });
    glp_scale_prob(prob, GLP_SF_AUTO);
    rc = glp_intopt(prob, iocp);
    if (rc == 0) {
      var sol = zeros(n);
      for (i = 0; i < n; i++) {
        sol[i] = glp_mip_col_val(prob, i + 1);
      }
      if (verbose) {
        var obj = glp_mip_obj_val(prob);
        console.log("Status : " + glp_mip_status(prob));
        console.log("Obj : " + obj);
      }
      return sol;
    } else return "Status : " + glp_get_prim_stat(prob);
  } else {
    var smcp = new SMCP({ presolve: GLP_ON });
    glp_scale_prob(prob, GLP_SF_AUTO);
    rc = glp_simplex(prob, smcp);
    if (rc == 0) {
      var sol = zeros(n);
      for (i = 0; i < n; i++) {
        sol[i] = glp_get_col_prim(prob, i + 1);
      }
      if (verbose) {
        var obj = glp_get_obj_val(prob);
        console.log(
          "Status : " +
            glp_get_status(prob) +
            "(OPT=" +
            GLP_OPT +
            ",FEAS=" +
            GLP_FEAS +
            ",INFEAS=" +
            GLP_INFEAS +
            ",NOFEAS=" +
            GLP_NOFEAS +
            ",UNBND=" +
            GLP_UNBND +
            ",UNDEF=" +
            GLP_UNDEF +
            ")"
        );
        console.log("Obj : " + obj);
      }
      return sol;
    } else {
      GLPLASTLP = "";
      glp_write_lp(prob, undefined, function (str) {
        GLPLASTLP += str + "<br>";
      });
      return (
        "RC=" +
        rc +
        " ; Status : " +
        glp_get_status(prob) +
        "(OPT=" +
        GLP_OPT +
        ",FEAS=" +
        GLP_FEAS +
        ",INFEAS=" +
        GLP_INFEAS +
        ",NOFEAS=" +
        GLP_NOFEAS +
        ",UNBND=" +
        GLP_UNBND +
        ",UNDEF=" +
        GLP_UNDEF +
        ")"
      );
    }
  }
};

function Matrix(m, n, values) {
  this.length = m;
  this.m = m;
  this.n = n;
  this.size = [m, n];
  this.type = "matrix";
  if (arguments.length == 2) this.val = new Float64Array(m * n);
  else if (arguments.length == 3) this.val = new Float64Array(values);
  else if (arguments.length == 4) this.val = values;
}
Matrix.prototype.get = function (i, j) {
  return this.val[i * this.n + j];
};
Matrix.prototype.set = function (i, j, v) {
  this.val[i * this.n + j] = v;
};
Matrix.prototype.row = function (i) {
  return this.val.subarray(i * this.n, (i + 1) * this.n);
};
Matrix.prototype.toArray = function () {
  var A = new Array(this.m);
  var ri = 0;
  for (var i = 0; i < this.m; i++) {
    A[i] = new Array(this.n);
    for (var j = 0; j < this.n; j++) A[i][j] = this.val[ri + j];
    ri += this.n;
  }
  return A;
};
Matrix.prototype.toArrayOfFloat64Array = function () {
  var A = new Array(this.m);
  for (var i = 0; i < this.m; i++)
    A[i] = this.val.subarray(i * this.n, (i + 1) * this.n);
  return A;
};

module.exports.Matrix = Matrix;
