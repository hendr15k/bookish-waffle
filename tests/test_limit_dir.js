const { CAS } = require('../js/cas.js');
const { Num, Sym, Div } = require('../js/expression.js');

const cas = new CAS();
const x = new Sym('x');
const zero = new Num(0);
const expr = new Div(new Num(1), x);

// Check limit(1/x, x, 0)
// With dir = 0 (default)
const lim0 = cas._limit(expr, x, zero, 0, 0);
console.log("limit(1/x, x, 0) = " + lim0.toString());

// Check limit(1/x, x, 0, 1) -> Infinity
const limPlus = cas._limit(expr, x, zero, 0, 1);
console.log("limit(1/x, x, 0, 1) = " + limPlus.toString());

// Check limit(1/x, x, 0, -1) -> -Infinity
const limMinus = cas._limit(expr, x, zero, 0, -1);
console.log("limit(1/x, x, 0, -1) = " + limMinus.toString());
