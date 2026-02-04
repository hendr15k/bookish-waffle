
const { CAS } = require('../js/cas');
const { Num, Sym, Call, Add, Mul } = require('../js/expression');

const cas = new CAS();

console.log("--- Testing Logic Operators ---");

const t = new Sym('true');
const f = new Sym('false');

const r1 = cas.evaluate(new Call('nand', [t, f]));
console.log(`nand(true, false) = ${r1.toString()} (Expected: 1)`);
if (r1.value !== 1) throw new Error("nand failed");

const r2 = cas.evaluate(new Call('nor', [f, f]));
console.log(`nor(false, false) = ${r2.toString()} (Expected: 1)`);
if (r2.value !== 1) throw new Error("nor failed");

const r3 = cas.evaluate(new Call('xnor', [t, f]));
console.log(`xnor(true, false) = ${r3.toString()} (Expected: 0)`);
if (r3.value !== 0) throw new Error("xnor failed");

console.log("\n--- Testing Differential Equations (Complex Roots) ---");

const x = new Sym('x');
// Use y(x) instead of y to prevent evaluation to 0 during differentiation
const y = new Call('y', [x]);

// 1. Harmonic Oscillator: y'' + y = 0
// desolve(diff(y(x), x, 2) + y(x), y(x))

const d2y = new Call('diff', [y, x, new Num(2)]);
const ode1 = new Add(d2y, y);
const res1 = cas.evaluate(new Call('desolve', [ode1, y]));
console.log(`desolve(y'' + y = 0) = ${res1.toString()}`);

if (!res1.toString().includes('cos(x)') || !res1.toString().includes('sin(x)')) {
    throw new Error("Harmonic oscillator failed to produce trig functions: " + res1.toString());
}

// 2. Damped Oscillator: y'' + 2y' + 2y = 0
// r = -1 +/- i. Alpha=-1, Beta=1.

const dy = new Call('diff', [y, x]);
const ode2 = new Add(new Add(d2y, new Mul(new Num(2), dy)), new Mul(new Num(2), y));
const res2 = cas.evaluate(new Call('desolve', [ode2, y]));
console.log(`desolve(y'' + 2y' + 2y = 0) = ${res2.toString()}`);

// Expect exp(-x) * ... or similar
// Note: Simplify might produce exp((-1 * x))
if (!res2.toString().includes('exp') || (!res2.toString().includes('-x') && !res2.toString().includes('(-1 * x)'))) {
     throw new Error("Damped oscillator missing exponential decay: " + res2.toString());
}
if (!res2.toString().includes('cos(x)') || !res2.toString().includes('sin(x)')) {
    throw new Error("Damped oscillator missing trig functions: " + res2.toString());
}

console.log("\n--- All Tests Passed ---");
