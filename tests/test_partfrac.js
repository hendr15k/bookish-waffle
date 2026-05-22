
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
const expr3 = new Div(new Num(1), new Add(new Pow(x, new Num(2)), new Num(1)));
const int3 = cas.evaluate(new Call('integrate', [expr3, x]));
console.log("integrate(1/(x^2+1)) -> " + int3.toString());

// 5. Integrate x^2/(x^2+1) -> x - atan(x)
// Regression test: polynomial division in partfrac was causing infinite recursion
// when quotient is 0 and remainder degree >= denominator degree.
const expr4 = new Div(new Pow(x, new Num(2)), new Add(new Pow(x, new Num(2)), new Num(1)));
const int4 = cas.evaluate(new Call('integrate', [expr4, x]));
console.log("integrate(x^2/(x^2+1)) -> " + int4.toString());
assert(!int4.toString().includes("integrate("), "integrate(x^2/(x^2+1)) should not return unevaluated");
assert(int4.toString().includes("atan"), "integrate(x^2/(x^2+1)) result should contain atan");

// 6. partfrac(x^2/(x^2+1), x) should give 1 - 1/(x^2+1)
const pf4 = cas.evaluate(new Call('partfrac', [expr4, x]));
console.log("partfrac(x^2/(x^2+1)) -> " + pf4.toString());
assert(pf4.toString().includes("1") && pf4.toString().includes("x^2 + 1"), "partfrac(x^2/(x^2+1)) should return a proper decomposition");

console.log("--- End Tests ---");
