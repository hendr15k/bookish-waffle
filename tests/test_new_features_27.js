
const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Lt, Gt } = require('../js/expression');

const cas = new CAS();

function assertApprox(actual, expected, message) {
    if (Math.abs(actual - expected) > 1e-6) {
        console.error(`FAIL: ${message}. Expected ${expected}, got ${actual}`);
        throw new Error(`FAIL: ${message}`);
    } else {
        console.log(`PASS: ${message}`);
    }
}

function assertString(actual, expected, message) {
    if (actual.replace(/\s/g, "") !== expected.replace(/\s/g, "")) {
        console.error(`FAIL: ${message}. Expected "${expected}", got "${actual}"`);
        throw new Error(`FAIL: ${message}`);
    } else {
        console.log(`PASS: ${message}`);
    }
}

console.log("--- Testing New Features 27 (Atan2, Airy, Piecewise Integration) ---");

// 1. Atan2
try {
    const a1 = cas.evaluate(new Call('atan2', [new Num(0), new Num(1)])); // 0
    assertString(a1.toString(), "0", "atan2(0, 1)");

    const a2 = cas.evaluate(new Call('atan2', [new Num(1), new Num(1)])); // pi/4
    assertString(a2.toString(), "(pi / 4)", "atan2(1, 1)");

    // Numeric
    const a3 = cas.evaluate(new Call('atan2', [new Num(1), new Num(2)])).evaluateNumeric();
    assertApprox(a3, Math.atan2(1, 2), "atan2(1, 2) numeric");

    // Derivative
    // d/dx atan2(y, x) = -y/(x^2+y^2)
    const dAtan2 = cas.evaluate(new Call('diff', [new Call('atan2', [new Sym('y'), new Sym('x')]), new Sym('x')]));
    // Expected: (-y / (x^2 + y^2))
    console.log("diff(atan2(y,x), x) -> " + dAtan2.toString());
    const valD = dAtan2.substitute(new Sym('x'), new Num(1)).substitute(new Sym('y'), new Num(1)).evaluateNumeric();
    assertApprox(valD, -0.5, "derivative of atan2 check");

} catch(e) {
    console.error("Atan2 Error:", e);
    throw e;
}

// 2. Airy Functions
try {
    const ai0 = cas.evaluate(new Call('airyAi', [new Num(0)])).evaluateNumeric();
    assertApprox(ai0, 0.35502805, "airyAi(0)");

    const bi0 = cas.evaluate(new Call('airyBi', [new Num(0)])).evaluateNumeric();
    assertApprox(bi0, 0.61492663, "airyBi(0)");

    const aip0 = cas.evaluate(new Call('airyAiPrime', [new Num(0)])).evaluateNumeric();
    assertApprox(aip0, -0.25881940, "airyAiPrime(0)");

    // Symbolic derivative
    const dAi = cas.evaluate(new Call('diff', [new Call('airyAi', [new Sym('x')]), new Sym('x')]));
    // Should return mul(airyAiPrime(x), 1) -> airyAiPrime(x)
    console.log("diff(airyAi) -> " + dAi.toString());
    // Simplify might leave it as product if not optimized, but Mul.simplify with 1 handles it.
    // Check if dAi string contains airyAiPrime
    if (!dAi.toString().includes("airyAiPrime")) throw new Error("Derivative of airyAi failed");

} catch(e) {
    console.error("Airy Error:", e);
    throw e;
}

// 3. Piecewise Integration
try {
    // integrate(piecewise(x < 0, 0, 1), x, -1, 1)
    // -1 to 0: 0. 0 to 1: 1. Total 1.
    // Note: Use Lt class, not Call
    const pw = new Call('piecewise', [new Lt(new Sym('x'), new Num(0)), new Num(0), new Num(1)]);
    const int1 = cas.evaluate(new Call('integrate', [pw, new Sym('x'), new Num(-1), new Num(1)]));
    console.log("Integrate Step: " + int1.toString());
    assertApprox(int1.evaluateNumeric(), 1, "Integrate Step Function");

    // integrate(piecewise(x < 0, x, x^2), x, -1, 1)
    // -1 to 0: x -> [x^2/2] -> 0 - 0.5 = -0.5
    // 0 to 1: x^2 -> [x^3/3] -> 1/3 - 0 = 0.333...
    // Sum = -0.5 + 0.3333 = -0.1666...
    const pw2 = new Call('piecewise', [new Lt(new Sym('x'), new Num(0)), new Sym('x'), new Pow(new Sym('x'), new Num(2))]);
    const int2 = cas.evaluate(new Call('integrate', [pw2, new Sym('x'), new Num(-1), new Num(1)]));
    console.log("Integrate Mixed: " + int2.toString());
    assertApprox(int2.evaluateNumeric(), -1.0/6.0, "Integrate Mixed Piecewise");

} catch(e) {
    console.error("Piecewise Integration Error:", e);
    throw e;
}

console.log("All tests passed");
