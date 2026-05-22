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

// map with expression placeholder x
const list2 = new Vec([new Num(1), new Num(2), new Num(3)]);
const t7b = new Call('map', [list2, new Add(new Sym('x'), new Num(1))]);
runTest("map([1, 2, 3], x+1)", t7b, "[2, 3, 4]");

// map with expression placeholder inside a call
const t7c = new Call('map', [list, new Call('sin', [new Sym('x')])]);
runTest("map([0, pi/2], sin(x))", t7c, "[0, 1]");

// map with underscore placeholder
const t7d = new Call('map', [list2, new Mul(new Sym('_'), new Num(2))]);
runTest("map([1, 2, 3], _*2)", t7d, "[2, 4, 6]");

// Test Eigenvalues 3x3 diagonal
// Diagonal 1..3
const rows = [];
for(let i=0; i<3; i++) {
    const row = [];
    for(let j=0; j<3; j++) row.push(new Num(i===j?i+1:0));
    rows.push(new Vec(row));
}
const D3 = new Vec(rows);
const t8 = new Call('eigenvals', [D3]);
runTest("eigenvals(D3)", t8, "[1, 2, 3]");

// 5. Bug fix: Pow.simplify with negative exponent (HEN-56)
// (c*x)^(-1) should give 1/(c*x), not fractional coefficient
// (2*x)^(-1) -> 1/(2*x)
const t9 = new Pow(new Mul(new Num(2), new Sym('x')), new Num(-1));
runTest("Pow Simplify (2*x)^(-1)", t9, "(1 / (2 * x))");

// (-2*x)^(-1) -> -1/(2*x) (simplified from 1/(-2*x))
const t10 = new Pow(new Mul(new Num(-2), new Sym('x')), new Num(-1));
runTest("Pow Simplify (-2*x)^(-1)", t10, "(-1 / (2 * x))");

// (c*x)^(-2) -> 1/(c^2 * x^2)
const t11 = new Pow(new Mul(new Num(2), new Sym('x')), new Num(-2));
runTest("Pow Simplify (2*x)^(-2)", t11, "(1 / (4 * x^2))");

// (c*x)^2 should still distribute (positive exponent)
const t12 = new Pow(new Mul(new Num(2), new Sym('x')), new Num(2));
runTest("Pow Simplify (2*x)^2", t12, "(4 * x^2)");

// HEN-90: Negative integer exponents → fractions
// x^(-1) -> 1/x
runTest("Pow x^(-1) -> 1/x", new Pow(new Sym('x'), new Num(-1)), "(1 / x)");

// x^(-2) -> 1/x^2
runTest("Pow x^(-2) -> 1/x^2", new Pow(new Sym('x'), new Num(-2)), "(1 / x^2)");

// x^(-3) -> 1/x^3
runTest("Pow x^(-3) -> 1/x^3", new Pow(new Sym('x'), new Num(-3)), "(1 / x^3)");

// HEN-90: Power-of-powers (a^x)^y -> a^(x*y)
// (a^x)^y -> a^(x*y)
const pow_base_a = new Sym('a');
const pow_exp_x = new Sym('x');
const pow_exp_y = new Sym('y');
const t_pow1 = new Pow(new Pow(pow_base_a, pow_exp_x), pow_exp_y);
const r_pow1 = cas.evaluate(t_pow1);
if (r_pow1 instanceof Pow && r_pow1.left.toString() === 'a' && r_pow1.right.toString() === '(x * y)') {
    console.log("PASS: Pow (a^x)^y -> a^(x*y)");
} else {
    console.log("FAIL: Pow (a^x)^y -> a^(x*y)");
    console.log("  Got:", r_pow1.toString(), "type:", r_pow1.constructor.name);
}

// (x^2)^3 -> x^6
const t_pow2 = new Pow(new Pow(new Sym('x'), new Num(2)), new Num(3));
runTest("Pow (x^2)^3 -> x^6", t_pow2, "x^6");

// ((a^x)^y)^z -> a^((x*y)*z)
const pow_exp_z = new Sym('z');
const inner = new Pow(new Pow(new Sym('a'), new Sym('x')), new Sym('y'));
const t_pow3 = new Pow(inner, pow_exp_z);
const r_pow3 = cas.evaluate(t_pow3);
if (r_pow3 instanceof Pow && r_pow3.left.toString() === 'a' && r_pow3.right.toString() === '((x * y) * z)') {
    console.log("PASS: Pow ((a^x)^y)^z -> a^((x*y)*z)");
} else {
    console.log("FAIL: Pow ((a^x)^y)^z -> a^((x*y)*z)");
    console.log("  Got:", r_pow3.toString(), "type:", r_pow3.constructor.name);
}

console.log("--- End Tests ---");
