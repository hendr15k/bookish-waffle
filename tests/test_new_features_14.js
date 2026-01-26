
const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec, Eq } = require('../js/expression');

// Helper to check equality
function assertEqual(actual, expected, message) {
    if (actual.toString() !== expected.toString()) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${expected.toString()}`);
        console.error(`  Actual:   ${actual.toString()}`);
    } else {
        console.log(`PASS: ${message}`);
    }
}

// Helper for approximate equality
function assertApprox(actual, expected, message) {
    const valA = actual.evaluateNumeric();
    const valB = expected.evaluateNumeric();
    if (Math.abs(valA - valB) > 1e-6) {
        console.error(`FAIL: ${message}`);
        console.error(`  Expected: ${valB}`);
        console.error(`  Actual:   ${valA}`);
    } else {
        console.log(`PASS: ${message}`);
    }
}

const cas = new CAS();

// 1. Lagrange Multipliers
console.log("\n--- Testing Lagrange Multipliers ---");
try {
    // Minimize x+y subject to x^2+y^2=1
    const f = new Add(new Sym('x'), new Sym('y'));
    const constraint = new Eq(new Add(new Pow(new Sym('x'), new Num(2)), new Pow(new Sym('y'), new Num(2))), new Num(1));
    const res = cas.evaluate(new Call('lagrange', [f, new Vec([constraint]), new Vec([new Sym('x'), new Sym('y')])]));
    console.log("Lagrange Result:", res.toString());

    // We expect solutions where x=y (or -y).
    // Just check if it returns a set of solutions.
    if (res instanceof Call && res.funcName === 'set') {
        console.log("PASS: Lagrange returned solution set");
    } else {
        console.error("FAIL: Lagrange did not return a set");
    }
} catch (e) {
    console.error("FAIL: Lagrange Error", e);
}

// 2. Mixed Partial Derivatives
console.log("\n--- Testing Mixed Partial Derivatives ---");
try {
    // diff(x*y^2*z, x, y, z)
    // d/dx -> y^2*z
    // d/dy -> 2*y*z
    // d/dz -> 2*y
    const expr = new Mul(new Sym('x'), new Mul(new Pow(new Sym('y'), new Num(2)), new Sym('z')));
    const diffMixed = cas.evaluate(new Call('diff', [expr, new Sym('x'), new Sym('y'), new Sym('z')]));

    const expected = new Mul(new Num(2), new Sym('y'));
    assertEqual(diffMixed.simplify(), expected, "diff(x*y^2*z, x, y, z)");

    // Test order argument preservation: diff(x^3, x, 2)
    const expr2 = new Pow(new Sym('x'), new Num(3));
    const diffOrder = cas.evaluate(new Call('diff', [expr2, new Sym('x'), new Num(2)]));
    const expectedOrder = new Mul(new Num(6), new Sym('x'));
    assertEqual(diffOrder.simplify(), expectedOrder, "diff(x^3, x, 2)");

} catch (e) {
    console.error("FAIL: Diff Error", e);
}

// 3. Multiple Integrals
console.log("\n--- Testing Multiple Integrals ---");
try {
    // integrate(x*y, x, 0, 1, y, 0, 1)
    const expr = new Mul(new Sym('x'), new Sym('y'));
    const intRes = cas.evaluate(new Call('integrate', [
        expr,
        new Sym('x'), new Num(0), new Num(1),
        new Sym('y'), new Num(0), new Num(1)
    ]));

    const expected = new Num(0.25);
    assertApprox(intRes, expected, "integrate(x*y, x, 0, 1, y, 0, 1)");

} catch (e) {
    console.error("FAIL: Integrate Error", e);
}

// 4. Unit Vector
console.log("\n--- Testing Unit Vector ---");
try {
    const v = new Vec([new Num(1), new Num(1)]); // [1, 1]
    const uv = cas.evaluate(new Call('unitVector', [v]));
    // Expected: [1/sqrt(2), 1/sqrt(2)]

    const elem0 = uv.elements[0].evaluateNumeric();
    if (Math.abs(elem0 - 1/Math.sqrt(2)) < 1e-6) {
        console.log("PASS: unitVector element 0 correct");
    } else {
        console.error(`FAIL: unitVector element 0. Got ${elem0}, expected ${1/Math.sqrt(2)}`);
    }

} catch (e) {
    console.error("FAIL: unitVector Error", e);
}

// 5. Trig Simplify
console.log("\n--- Testing Trig Simplify ---");
try {
    const expr = new Add(new Pow(new Call('sin', [new Sym('x')]), new Num(2)), new Pow(new Call('cos', [new Sym('x')]), new Num(2)));
    const simp = cas.evaluate(new Call('trigSimplify', [expr]));

    // Should be 1
    assertEqual(simp, new Num(1), "trigSimplify(sin^2 + cos^2)");

} catch (e) {
    console.error("FAIL: trigSimplify Error", e);
}
