const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Vec, Eq } = require('../js/expression.js');

// Expose to global for CAS
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.BinaryOp = BinaryOp;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;
global.Call = Call;
global.Vec = Vec;
global.Eq = Eq;

const { CAS } = require('../js/cas.js');
const cas = new CAS();

console.log("Running tests for rk4 and gradient_descent...");

let passed = 0;
let total = 0;

function assertClose(actual, expected, tol=1e-3, msg="") {
    total++;
    let val = actual;
    if (val instanceof Num) val = val.value;

    if (Math.abs(val - expected) < tol) {
        passed++;
        console.log(`PASS: ${msg} (${val} approx ${expected})`);
    } else {
        console.log(`FAIL: ${msg} (Expected ${expected}, got ${val})`);
    }
}

function assertVecClose(actual, expected, tol=1e-3, msg="") {
    total++;
    let passedLocal = true;
    if (!(actual instanceof Vec) || !(expected instanceof Vec)) {
        console.log(`FAIL: ${msg} (Type mismatch)`);
        return;
    }
    if (actual.elements.length !== expected.elements.length) {
        console.log(`FAIL: ${msg} (Length mismatch)`);
        return;
    }

    for(let i=0; i<actual.elements.length; i++) {
        let v1 = actual.elements[i].evaluateNumeric();
        let v2 = expected.elements[i].evaluateNumeric();
        if (Math.abs(v1 - v2) > tol) {
            passedLocal = false;
            console.log(`FAIL: ${msg} element ${i} (Expected ${v2}, got ${v1})`);
        }
    }

    if (passedLocal) {
        passed++;
        console.log(`PASS: ${msg}`);
    }
}

// 1. Test rk4: y' = y, y(0) = 1. y(1) should be e.
// rk4(diffEq, depVar, indepVar, t0, y0, t1, step)
try {
    const y = new Sym('y');
    const t = new Sym('t');
    const diffEq = y; // y' = y
    const res = cas.evaluate(new Call('rk4', [diffEq, y, t, new Num(0), new Num(1), new Num(1), new Num(0.1)]));

    // Last point
    const points = res.elements;
    const lastPt = points[points.length - 1]; // [t, y]
    const lastY = lastPt.elements[1];

    assertClose(lastY, Math.E, 1e-2, "rk4(y'=y, y(0)=1) at t=1");
} catch(e) {
    console.log("FAIL: rk4 error", e);
}

// 2. Test gradient_descent: f(x) = (x-2)^2. Min at x=2.
try {
    const x = new Sym('x');
    const f = new Pow(new Sub(x, new Num(2)), new Num(2));
    // gradient_descent(func, vars, start, alpha, iter)
    // vars must be Vec or Sym. start must be numeric.
    const res = cas.evaluate(new Call('gradient_descent', [f, x, new Num(0)]));

    assertClose(res, 2, 1e-2, "gradient_descent((x-2)^2)");
} catch(e) {
    console.log("FAIL: gradient_descent 1D error", e);
}

// 3. Test gradient_descent 2D: f(x,y) = (x-1)^2 + (y+3)^2. Min at (1, -3).
try {
    const x = new Sym('x');
    const y = new Sym('y');
    const f = new Add(new Pow(new Sub(x, new Num(1)), new Num(2)), new Pow(new Add(y, new Num(3)), new Num(2)));
    const vars = new Vec([x, y]);
    const start = new Vec([new Num(0), new Num(0)]);

    const res = cas.evaluate(new Call('gradient_descent', [f, vars, start, new Num(0.1), new Num(100)]));

    assertVecClose(res, new Vec([new Num(1), new Num(-3)]), 1e-2, "gradient_descent 2D");
} catch(e) {
    console.log("FAIL: gradient_descent 2D error", e);
}

console.log(`\nPassed ${passed} / ${total} tests.`);
if (passed === total) {
    process.exit(0);
} else {
    process.exit(1);
}
