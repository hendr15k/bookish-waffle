const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Call, Vec } = require('../js/expression');

const cas = new CAS();
const assert = (cond, msg) => {
    if (!cond) {
        console.error("FAILED: " + msg);
        process.exit(1);
    } else {
        console.log("PASSED: " + msg);
    }
};

const assertApprox = (val, expected, tol=1e-6, msg="") => {
    if (Math.abs(val - expected) > tol) {
        console.error(`FAILED: ${msg} (Got ${val}, Expected ${expected})`);
        process.exit(1);
    } else {
        console.log(`PASSED: ${msg}`);
    }
};

console.log("--- Test New Features 8: Distributions & Special Functions ---");

// 1. Sinc
const s0 = cas.evaluate(new Call('sinc', [new Num(0)]));
assert(s0.value === 1, "sinc(0) should be 1");
const s1 = cas.evaluate(new Call('sinc', [new Num(Math.PI/2)])); // sin(pi/2)/(pi/2) = 1/(pi/2) = 2/pi
assertApprox(s1.value, 2/Math.PI, 1e-9, "sinc(pi/2)");

// 2. LambertW
const w_e = cas.evaluate(new Call('LambertW', [new Num(Math.E)]));
assertApprox(w_e.value, 1, 1e-6, "LambertW(e) approx 1");
const w_0 = cas.evaluate(new Call('LambertW', [new Num(0)]));
assertApprox(w_0.value, 0, 1e-9, "LambertW(0) = 0");

// 3. Inverse Erf
const invErf = cas.evaluate(new Call('erfinv', [new Num(0.520499877)])); // erf(0.5) approx 0.5205
assertApprox(invErf.value, 0.5, 1e-4, "erfinv(erf(0.5)) approx 0.5");

// 4. Jacobi Symbol
const j1 = cas.evaluate(new Call('jacobiSymbol', [new Num(2), new Num(7)])); // 3^2 = 9 = 2 mod 7 -> 1
assert(j1.value === 1, "jacobiSymbol(2, 7) = 1");
const j2 = cas.evaluate(new Call('jacobiSymbol', [new Num(3), new Num(7)])); // 3 is non-residue mod 7 -> -1
assert(j2.value === -1, "jacobiSymbol(3, 7) = -1");

// 5. ChebyshevU (Second kind)
// U_0(x) = 1
// U_1(x) = 2x
// U_2(x) = 4x^2 - 1
const u2 = cas.evaluate(new Call('chebyshevU', [new Num(2), new Sym('x')]));
// We expand to ensure canonical form 4x^2 - 1
const u2_expanded = u2.expand().simplify();
const u2_str = u2_expanded.toString().replace(/\s/g, '');
assert(u2_str === "((4*(x^2))-1)" || u2_str === "4*x^2-1" || u2_str === "(((2*x))^2-1)",
       `chebyshevU(2, x) expanded structure. Got ${u2_str}`);

// 6. Bernoulli Polynomial
// B_2(x) = x^2 - x + 1/6
const b2 = cas.evaluate(new Call('bernoulliPoly', [new Num(2), new Sym('x')]));
// We check coefficients or value at x=0, x=1
const b2_0 = b2.substitute(new Sym('x'), new Num(0)).simplify();
assertApprox(b2_0.value, 1/6, 1e-9, "B_2(0) = 1/6");
const b2_1 = b2.substitute(new Sym('x'), new Num(1)).simplify();
assertApprox(b2_1.value, 1/6, 1e-9, "B_2(1) = 1/6");

// 7. Gamma Distribution
// gammaPDF(x, k, theta) = x^(k-1) * exp(-x/theta) / (theta^k * gamma(k))
// k=2, theta=1, x=1 -> 1 * e^-1 / (1 * 1) = 1/e
const gPDF = cas.evaluate(new Call('gammaPDF', [new Num(1), new Num(2), new Num(1)]));
assertApprox(gPDF.value, 1/Math.E, 1e-6, "gammaPDF(1, 2, 1)");

// 8. Uniform Distribution
const uPDF = cas.evaluate(new Call('uniformPDF', [new Num(0.5), new Num(0), new Num(1)]));
assert(uPDF.value === 1, "uniformPDF(0.5, 0, 1) = 1");
const uCDF = cas.evaluate(new Call('uniformCDF', [new Num(0.5), new Num(0), new Num(1)]));
assert(uCDF.value === 0.5, "uniformCDF(0.5, 0, 1) = 0.5");

console.log("All tests passed!");
