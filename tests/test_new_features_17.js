const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Vec } = require('../js/expression.js');

// Setup global for CAS to see classes
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.Mul = Mul;
global.Pow = Pow;
global.Add = Add;
global.Div = Div;
global.Vec = Vec;
global.Call = Call;

const cas = new CAS();

function runTest(name, input, expected) {
    let res;
    try {
        if (typeof input === 'string') {
            // Parsing not available here, assume manual construction or evaluate symbolic
            throw new Error("String input not supported in this test harness without Parser");
        }
        res = cas.evaluate(input);
    } catch (e) {
        console.log(`FAIL: ${name} - Error: ${e.message}`);
        return;
    }

    const resStr = res.toString();
    if (resStr === expected) {
        console.log(`PASS: ${name}`);
    } else {
        console.log(`FAIL: ${name}`);
        console.log(`  Input: ${input.toString()}`);
        console.log(`  Expected: ${expected}`);
        console.log(`  Got: ${resStr}`);
    }
}

console.log("--- Test Set 17: Bug Fixes and New Features ---");

// 1. Pow.simplify fix: (-2x)^2 -> 4x^2
// Input: Pow(Mul(-2, x), 2)
const t1 = new Pow(new Mul(new Num(-2), new Sym('x')), new Num(2));
// Note: recursiveEval calls simplify(). Pow construction doesn't simplify immediately unless called.
runTest("Pow Simplify (-2x)^2", t1, "(4 * x^2)");

// 2. rect(r, theta)
// rect(2, pi) -> 2 * (cos(pi) + i*sin(pi)) -> 2 * (-1 + 0) -> -2
// Input: Call('rect', [2, pi])
const t2 = new Call('rect', [new Num(2), new Sym('pi')]);
runTest("rect(2, pi)", t2, "-2");

// rect(1, pi/2) -> 1 * (0 + i) -> i
const t3 = new Call('rect', [new Num(1), new Div(new Sym('pi'), new Num(2))]);
runTest("rect(1, pi/2)", t3, "i");

// 3. clamp(x, min, max)
// clamp(5, 0, 10) -> 5
const t4 = new Call('clamp', [new Num(5), new Num(0), new Num(10)]);
runTest("clamp(5, 0, 10)", t4, "5");

// clamp(-5, 0, 10) -> 0
const t5 = new Call('clamp', [new Num(-5), new Num(0), new Num(10)]);
runTest("clamp(-5, 0, 10)", t5, "0");

// clamp(15, 0, 10) -> 10
const t6 = new Call('clamp', [new Num(15), new Num(0), new Num(10)]);
runTest("clamp(15, 0, 10)", t6, "10");

// 4. map(list, func)
// map([0, pi/2], sin) -> [sin(0), sin(pi/2)] -> [0, 1]
const list = new Vec([new Num(0), new Div(new Sym('pi'), new Num(2))]);
const t7 = new Call('map', [list, new Sym('sin')]);
runTest("map([0, pi/2], sin)", t7, "[0, 1]");

// Test Eigenvalues Numeric 5x5
// Diagonal 1..5
const rows = [];
for(let i=0; i<5; i++) {
    const row = [];
    for(let j=0; j<5; j++) row.push(new Num(i===j?i+1:0));
    rows.push(new Vec(row));
}
const D5 = new Vec(rows);
const t8 = new Call('eigenvals', [D5]);
runTest("eigenvals(D5)", t8, "[1, 2, 3, 4, 5]");

console.log("--- End Tests ---");
