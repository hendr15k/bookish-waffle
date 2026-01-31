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
            throw new Error("String input not supported in this test harness without Parser");
        }
        res = cas.evaluate(input);
    } catch (e) {
        console.log(`FAIL: ${name} - Error: ${e.message}`);
        return;
    }

    // Normalize output by removing spaces for comparison
    const resStr = res.toString().replace(/\s+/g, '');
    const expStr = expected.replace(/\s+/g, '');

    if (resStr === expStr) {
        console.log(`PASS: ${name}`);
    } else {
        console.log(`FAIL: ${name}`);
        console.log(`  Input: ${input.toString()}`);
        console.log(`  Expected: ${expected}`);
        console.log(`  Got: ${res.toString()}`);
    }
}

console.log("--- Test Set 18: Integrals and Inverse ---");

const x = new Sym('x');

// 1. Integrals of Inverse Hyperbolic
// integrate(asinh(x), x) -> x*asinh(x) - sqrt(x^2+1)
const t1 = new Call('integrate', [new Call('asinh', [x]), x]);
// Expected: x*asinh(x) - sqrt(x^2+1)
// Note: Parentheses are strict in toString
runTest("integrate(asinh(x))", t1, "((x * asinh(x)) - sqrt((x^2 + 1)))");

// integrate(acosh(x), x) -> x*acosh(x) - sqrt(x^2-1)
const t2 = new Call('integrate', [new Call('acosh', [x]), x]);
runTest("integrate(acosh(x))", t2, "((x * acosh(x)) - sqrt((x^2 - 1)))");

// integrate(atanh(x), x) -> x*atanh(x) + 0.5*ln(1-x^2)
// Note: 0.5*ln(u) simplifies to ln(u^0.5)
const t3 = new Call('integrate', [new Call('atanh', [x]), x]);
runTest("integrate(atanh(x))", t3, "((x * atanh(x)) + ln(((1 - x^2))^0.5))");

// 2. Integrals of Inverse Trig (Reciprocal)
// integrate(acot(x), x) -> x*acot(x) + 0.5*ln(1+x^2)
const t4 = new Call('integrate', [new Call('acot', [x]), x]);
runTest("integrate(acot(x))", t4, "((x * acot(x)) + ln(((1 + x^2))^0.5))");

// integrate(asec(x), x) -> x*asec(x) - acosh(x)
const t5 = new Call('integrate', [new Call('asec', [x]), x]);
runTest("integrate(asec(x))", t5, "((x * asec(x)) - acosh(x))");

// integrate(acsc(x), x) -> x*acsc(x) + acosh(x)
const t6 = new Call('integrate', [new Call('acsc', [x]), x]);
runTest("integrate(acsc(x))", t6, "((x * acsc(x)) + acosh(x))");


// 3. Symbolic Matrix Inverse (2x2)
// [[a, b], [c, d]]^-1
// Determinant = ad - bc
// Result = 1/(ad-bc) * [[d, -b], [-c, a]]
const a = new Sym('a');
const b = new Sym('b');
const c = new Sym('c');
const d = new Sym('d');
const M = new Vec([
    new Vec([a, b]),
    new Vec([c, d])
]);

const invM = new Call('inv', [M]);
// We expect the explicit formula, simplified.
const resInv = cas.evaluate(invM);

if (resInv instanceof Vec && resInv.elements.length === 2) {
    const el00 = resInv.elements[0].elements[0].toString().replace(/\s+/g, '');
    const el01 = resInv.elements[0].elements[1].toString().replace(/\s+/g, '');

    // Check if it matches expected structure roughly
    // d / (a*d - b*c)
    // Note: simplify might reorder (a*d - b*c) to (a*d - c*b) etc.
    const det = "(a*d-b*c)";
    const expected00 = `(d/${det})`;

    if (el00.includes("/(a*d") || el00.includes("/((a*d")) {
         console.log("PASS: 2x2 Inverse Structure");
    } else {
         console.log("FAIL: 2x2 Inverse Structure");
         console.log("Got 00:", el00);
    }
} else {
    console.log("FAIL: 2x2 Inverse did not return 2x2 matrix");
}

console.log("--- End Tests ---");
