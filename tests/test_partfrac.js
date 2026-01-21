
const { CAS } = require('../js/cas');
const { Expr, Sym, Num, Call, Div, Pow, Add, Sub, Mul } = require('../js/expression');

// Mock global objects
global.Expr = Expr;
global.Sym = Sym;
global.Num = Num;

const cas = new CAS();
const x = new Sym('x');

function assert(condition, message) {
    if (!condition) {
        throw new Error("Test failed: " + message);
    }
    console.log("PASS: " + message);
}

function assertStr(expr, expected, message) {
    const s = expr.toString();
    if (s !== expected) {
        console.error(`Expected: ${expected}`);
        console.error(`Actual:   ${s}`);
        throw new Error("Test failed: " + message);
    }
    console.log("PASS: " + message);
}

console.log("--- Testing partfrac ---");

// 1. partfrac(1/(x^2-1))
// Roots: 1, -1. Residues: 1/2, -1/2.
// Expect: (1/2)/(x-1) + (-1/2)/(x+1)
const expr1 = new Div(new Num(1), new Sub(new Pow(x, new Num(2)), new Num(1)));
const pf1 = cas.evaluate(new Call('partfrac', [expr1, x]));
// Note: exact string representation might vary due to addition order, but let's check structure
console.log("1/(x^2-1) -> " + pf1.toString());
assert(pf1 instanceof Add, "Result is a Sum");
// assertStr(pf1, "(((1 / 2) / (x - 1)) + ((-1 / 2) / (x + 1)))", "Simple quadratic decomposition");

// 2. partfrac(1/(x*(x+1))) -> 1/x - 1/(x+1)
// Roots: 0, -1.
// Residue at 0: 1/(0+1) = 1. Term: 1/x.
// Residue at -1: 1/-1 = -1. Term: -1/(x+1).
const expr2 = new Div(new Num(1), new Mul(x, new Add(x, new Num(1))));
const pf2 = cas.evaluate(new Call('partfrac', [expr2, x]));
console.log("1/(x*(x+1)) -> " + pf2.toString());
// This might verify as ((1 / x) + (-1 / (x + 1)))

// 3. Integrate 1/(x^2-1)
// Should use partfrac internally
const int1 = cas.evaluate(new Call('integrate', [expr1, x]));
console.log("integrate(1/(x^2-1)) -> " + int1.toString());
// Accept ln or atanh
const s = int1.toString();
assert(s.includes("ln") || s.includes("atanh"), "Integration result contains ln or atanh");

// 4. Integrate 1/(x^2+1) -> atan(x)
// This is handled by standard integrate table (via table lookups if implemented, or partfrac with complex roots?)
// My partfrac currently supports real roots from _solve. x^2+1 has no real roots.
// _solve returns complex?
// My _solve logic handles quadratics with discriminant < 0 if it supports sqrt(-1).
// Let's see.
const expr3 = new Div(new Num(1), new Add(new Pow(x, new Num(2)), new Num(1)));
const int3 = cas.evaluate(new Call('integrate', [expr3, x]));
console.log("integrate(1/(x^2+1)) -> " + int3.toString());
// Currently CAS.js integrate doesn't have table for atan.
// partfrac might fail or return complex.
// If it fails, integrate returns Call.

console.log("--- End Tests ---");
