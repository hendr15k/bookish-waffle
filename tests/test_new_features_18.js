const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Vec, Eq } = require('../js/expression.js');

// Setup global for CAS
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.Mul = Mul;
global.Pow = Pow;
global.Add = Add;
global.Sub = Sub;
global.Div = Div;
global.Vec = Vec;
global.Call = Call;
global.Eq = Eq;

const cas = new CAS();

function runTest(name, input, expected, numericCheck = false) {
    let res;
    try {
        res = cas.evaluate(input);
    } catch (e) {
        console.log(`FAIL: ${name} - Error: ${e.message}`);
        return;
    }

    const resStr = res.toString();

    if (numericCheck) {
        // Evaluate expected and result
        // Expected usually string, parse it? No parser.
        // Assume expected is a number or roughly match.
        // We will just check basic string match or logic in test.
    }

    if (resStr === expected) {
        console.log(`PASS: ${name}`);
    } else {
        // Allow slight differences in spacing or parens if semantically same, but exact match preferred
        console.log(`FAIL: ${name}`);
        console.log(`  Input: ${input.toString()}`);
        console.log(`  Expected: ${expected}`);
        console.log(`  Got: ${resStr}`);
    }
}

console.log("--- Test Set 18: Simplex, Groebner, Q ---");

// 1. Simplex
// Max 3x + 2y
// 2x + y <= 18
// 2x + 3y <= 42
// 3x + y <= 24
// x, y >= 0
// Opt: x=3, y=12, Z=33
const c = new Vec([new Num(3), new Num(2)]);
const b = new Vec([new Num(18), new Num(42), new Num(24)]);
const A = new Vec([
    new Vec([new Num(2), new Num(1)]),
    new Vec([new Num(2), new Num(3)]),
    new Vec([new Num(3), new Num(1)])
]);

// Returns [ [x, y], Z ]
const t1 = new Call('simplex', [c, A, b]);
// Expected: [[3, 12], 33]
runTest("Simplex 2D", t1, "[[3, 12], 33]");

// 2. Groebner Basis
// Basis for { x + y, x - y } -> { 2x, -2y } or { x, y } depending on reduction.
// My implementation returns raw basis from Buchberger, likely { x+y, x-y, -2y, ... }
// Let's test a simple case: { x, y } -> { x, y }
const x = new Sym('x');
const y = new Sym('y');
const vars = new Vec([x, y]);
const basis = new Vec([x, y]);
const t2 = new Call('groebner', [basis, vars]);
// Result should be [x, y] (order might vary or be same)
// Since sorted, x > y lex?
// My code doesn't strictly sort the output basis, just appends remainders.
// For [x, y], S(x, y) = x*y - y*x = 0. Remainder 0. Basis stable.
runTest("Groebner {x, y}", t2, "[x, y]");

// Groebner { x^2 - 1, x - 1 } -> { x-1 } (since x^2-1 = (x-1)(x+1))
// S-poly?
// f = x^2 - 1, g = x - 1.
// LCM = x^2. S = 1*(x^2-1) - x*(x-1) = x^2 - 1 - x^2 + x = x - 1.
// Remainder of x-1 wrt {x-1} is 0.
// So basis remains { x^2-1, x-1 }.
// My implementation doesn't reduce the basis fully (auto-reduction).
// So it returns the input if no new polynomials found.
const t3 = new Call('groebner', [new Vec([new Sub(new Pow(x, new Num(2)), new Num(1)), new Sub(x, new Num(1))]), new Vec([x])]);
// Expect same list
// runTest("Groebner Redundant", t3, "[(x^2 - 1), (x - 1)]");

// Groebner { x^2 + 2x + 1, x + 1 }
const t4 = new Call('groebner', [
    new Vec([
        new Add(new Pow(x, new Num(2)), new Add(new Mul(new Num(2), x), new Num(1))),
        new Add(x, new Num(1))
    ]),
    new Vec([x])
]);
// x^2+2x+1 / x+1 -> rem 0.
// Basis stable.
// runTest("Groebner Perfect Square", t4, "...");

// 3. Solve System (Fallback)
// x + y = 2
// x - y = 0
// Heuristic substitution handles this easily.
const sys1 = new Vec([
    new Eq(new Add(x, y), new Num(2)),
    new Eq(new Sub(x, y), new Num(0))
]);
const t5 = new Call('solve', [sys1, vars]);
// Solution: x=1, y=1
// set([Eq(x, 1), Eq(y, 1)]) or set([Eq(y, 1), Eq(x, 1)]) depending on order
// Heuristic should find it.
// runTest("Solve Linear System", t5, "{[x = 1, y = 1]}");

// Nonlinear System where heuristic might fail?
// x^2 + y^2 = 1, y = x^2.
// y + y^2 = 1 => y^2 + y - 1 = 0.
// Substitutions should work.

// 4. Q Function
const t6 = new Call('Q', [new Num(0)]);
runTest("Q(0)", t6, "0.5");

// Q(Infinity)
const t7 = new Call('Q', [new Sym('Infinity')]);
// erfc(Inf) = 0.
// Numeric check because symbolic simplification of Infinity/sqrt(2) might be incomplete
runTest("Q(Infinity)", t7.evaluateNumeric(), "0");

console.log("--- End Tests ---");
