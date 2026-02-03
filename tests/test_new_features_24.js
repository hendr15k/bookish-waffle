
const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow } = require('../js/expression.js');

const cas = new CAS();

console.log("Starting tests...");

// 1. Controllability
// A = [[0, 1], [-2, -3]], B = [[0], [1]]
// C = [B, AB] = [[0, 1], [1, -3]] -> Rank 2 -> Controllable
try {
    const A = new Vec([new Vec([new Num(0), new Num(1)]), new Vec([new Num(-2), new Num(-3)])]);
    const B = new Vec([new Vec([new Num(0)]), new Vec([new Num(1)])]);

    const ctrb = cas.evaluate(new Call('ctrb', [A, B]));
    console.log("ctrb:", ctrb.toString());

    const isC = cas.evaluate(new Call('isControllable', [A, B]));
    console.log("isControllable:", isC.toString());
    if (isC.value !== 1) throw new Error("System should be controllable");
} catch (e) {
    console.error("Controllability test failed:", e);
    process.exit(1);
}

// 2. Observability
// C = [[1, 0]]
try {
    const A = new Vec([new Vec([new Num(0), new Num(1)]), new Vec([new Num(-2), new Num(-3)])]);
    const C_mat = new Vec([new Vec([new Num(1), new Num(0)])]);

    const obsv = cas.evaluate(new Call('obsv', [A, C_mat]));
    console.log("obsv:", obsv.toString());

    const isO = cas.evaluate(new Call('isObservable', [A, C_mat]));
    console.log("isObservable:", isO.toString());
    if (isO.value !== 1) throw new Error("System should be observable");
} catch (e) {
    console.error("Observability test failed:", e);
    process.exit(1);
}

// 3. Ackermann
// Poles at -1, -2. Phi(s) = (s+1)(s+2) = s^2 + 3s + 2.
// A = [[0, 1], [0, 0]], B = [[0], [1]] (Integrator chain)
// K = [2, 3]
try {
    const A = new Vec([new Vec([new Num(0), new Num(1)]), new Vec([new Num(0), new Num(0)])]);
    const B = new Vec([new Vec([new Num(0)]), new Vec([new Num(1)])]);
    const poles = new Vec([new Num(-1), new Num(-2)]);

    const K = cas.evaluate(new Call('ackermann', [A, B, poles]));
    console.log("Ackermann K:", K.toString());

    const k1 = K.elements[0].elements[0].value;
    const k2 = K.elements[0].elements[1].value;
    if (Math.abs(k1 - 2) > 1e-9 || Math.abs(k2 - 3) > 1e-9) throw new Error("Ackermann failed");
} catch (e) {
    console.error("Ackermann test failed:", e);
    process.exit(1);
}

// 4. SS2TF
// sI - A = [[s, -1], [0, s]]
// C = [1, 0], B = [0; 1]
// TF = 1/s^2
try {
    const A = new Vec([new Vec([new Num(0), new Num(1)]), new Vec([new Num(0), new Num(0)])]);
    const B = new Vec([new Vec([new Num(0)]), new Vec([new Num(1)])]);
    const C = new Vec([new Vec([new Num(1), new Num(0)])]);
    const D = new Vec([new Vec([new Num(0)])]); // 1x1 D

    const tf = cas.evaluate(new Call('ss2tf', [A, B, C, D]));
    console.log("TF:", tf.toString());

    const s = new Sym('s');
    const val = tf.substitute(s, new Num(2)).evaluateNumeric();
    if (Math.abs(val - 0.25) > 1e-9) throw new Error("SS2TF failed: " + val);
} catch (e) {
    console.error("SS2TF test failed:", e);
    process.exit(1);
}

// 5. Fourier
// f(t) = t, L=1. Series on [-1, 1].
try {
    const t = new Sym('t');
    const expr = t;
    const n = new Num(1);
    const L = new Num(1);

    const series = cas.evaluate(new Call('fourier', [expr, t, n, L]));
    console.log("Fourier:", series.toString());

    // Check coefficient of sin(pi*t) -> 2/pi
    const valAtHalf = series.substitute(t, new Num(0.5)).evaluateNumeric(); // sin(pi/2)=1
    const expected = 2 / Math.PI;
    if (Math.abs(valAtHalf - expected) > 1e-5) throw new Error("Fourier failed: " + valAtHalf + " vs " + expected);
} catch (e) {
    console.error("Fourier test failed:", e);
    process.exit(1);
}

// 6. Partfrac Repeated Roots
// 1 / ((x-1)^2 * (x+1))
try {
    const x = new Sym('x');
    const num = new Num(1);
    const den = new Mul(new Pow(new Sub(x, new Num(1)), new Num(2)), new Add(x, new Num(1)));
    const expr = new Div(num, den);

    const pf = cas.evaluate(new Call('partfrac', [expr, x]));
    console.log("Partfrac:", pf.toString());

    // Check value at x=0 => 1
    const val = pf.substitute(x, new Num(0)).evaluateNumeric();
    if (Math.abs(val - 1) > 1e-9) throw new Error("Partfrac failed value check: " + val);
} catch (e) {
    console.error("Partfrac test failed:", e);
    process.exit(1);
}

console.log("All tests passed!");
