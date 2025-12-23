const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Div } = require('../js/expression.js');

// Mock globalThis for tests if needed
globalThis.HELP_DATA = {};

const cas = new CAS();
const parser = { parse: (s) => { throw new Error("Parser not available in this test"); } }; // Mock parser

// Construct limit(1/x, x, 0, 1) -> Infinity
// Construct limit(1/x, x, 0, -1) -> -Infinity

// Manually construct expressions
const x = new Sym('x');
const expr = new Div(new Num(1), x);
const point = new Num(0);

// Test Right Limit
// limit(1/x, x, 0, 1)
const limitRight = cas._limit(expr, x, point, 0, 1);
console.log("Limit 1/x -> 0+ : " + limitRight.toString());

// Test Left Limit
// limit(1/x, x, 0, -1)
const limitLeft = cas._limit(expr, x, point, 0, -1);
console.log("Limit 1/x -> 0- : " + limitLeft.toString());

// Test Default Limit (should be unsigned Infinity for 1/x?)
// limit(1/x, x, 0, 0)
const limitDefault = cas._limit(expr, x, point, 0, 0);
console.log("Limit 1/x -> 0  : " + limitDefault.toString());
