
const { CAS } = require('../js/cas');
const { Expr, Sym, Num, Call, Div, Pow, Add, Sub, Mul } = require('../js/expression');

// Mock global objects
global.Expr = Expr;
global.Sym = Sym;
global.Num = Num;
global.Call = Call;
global.Div = Div;
global.Pow = Pow;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;

const cas = new CAS();
const x = new Sym('x');
const a = new Sym('a');

function assert(condition, message) {
    if (!condition) {
        console.error("FAIL: " + message);
        process.exit(1);
    }
    console.log("PASS: " + message);
}

console.log("--- Testing Symbolic Partial Fraction Decomposition ---");

// Test 1: partfrac(1/(x-a)^2, x)
const expr1 = new Div(new Num(1), new Pow(new Sub(x, a), new Num(2)));
console.log("Expr1: " + expr1.toString());
const pf1 = cas.evaluate(new Call('partfrac', [expr1, x]));
console.log("Result1: " + pf1.toString());
// Should be 1/(x-a)^2 or similar structure
assert(pf1.toString().includes("(x - a))^2"), "Result1 should contain (x-a)^2 denominator");
assert(!pf1.toString().includes("NaN"), "Result1 should not be NaN");
assert(pf1.toString() !== "0", "Result1 should not be 0");

// Test 2: partfrac(1/((x-1)*(x-a)), x)
const expr2 = new Div(new Num(1), new Mul(new Sub(x, new Num(1)), new Sub(x, a)));
console.log("Expr2: " + expr2.toString());
const pf2 = cas.evaluate(new Call('partfrac', [expr2, x]));
console.log("Result2: " + pf2.toString());
// Expect sum of two terms
assert(pf2 instanceof Add, "Result2 should be an Add expression");
// Check denominators
const s2 = pf2.toString();
assert(s2.includes("x - 1") && s2.includes("x - a"), "Result2 should contain x-1 and x-a denominators");

// Test 3: Integration using partfrac for symbolic roots
// integrate(1/((x-1)(x-a)), x)
const intExpr = cas.evaluate(new Call('integrate', [expr2, x]));
console.log("Integrate Result: " + intExpr.toString());
const s3 = intExpr.toString();
assert(s3.includes("ln"), "Integration should result in logarithms");

console.log("--- End Tests ---");
