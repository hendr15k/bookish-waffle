const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Vec, Mul, Add, Sub, Div, Pow } = require('../js/expression.js');

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        throw new Error(message || "Assertion failed");
    }
}

function testRK4() {
    console.log("Testing rk4...");
    // y' = -y, y(0)=1. Solution y = exp(-t).
    // ode = -y
    const y = new Sym('y');
    const t = new Sym('t');
    const ode = new Mul(new Num(-1), y);

    const res = cas.evaluate(new Call('rk4', [
        ode,
        y,
        t,
        new Num(1),
        new Num(0),
        new Num(1),
        new Num(0.1)
    ]));

    assert(res instanceof Vec, "Result should be a vector");
    const lastPoint = res.elements[res.elements.length - 1];
    // [t, y]
    const tVal = lastPoint.elements[0].value;
    const yVal = lastPoint.elements[1].value;

    assert(Math.abs(tVal - 1.0) < 1e-9, "Time end correct");
    assert(Math.abs(yVal - Math.exp(-1)) < 1e-4, `Value accuracy scalar: ${yVal} vs ${Math.exp(-1)}`);

    // System: y1' = y2, y2' = -y1 (Harmonic oscillator)
    const y1 = new Sym('y1');
    const y2 = new Sym('y2');
    // ode = [y2, -y1]
    const odeSys = new Vec([y2, new Mul(new Num(-1), y1)]);
    const vars = new Vec([y1, y2]);
    const init = new Vec([new Num(0), new Num(1)]); // y1(0)=0, y2(0)=1 -> y1=sin(t)

    const resSys = cas.evaluate(new Call('rk4', [
        odeSys,
        vars,
        t,
        init,
        new Num(0),
        new Num(Math.PI/2),
        new Num(0.01)
    ]));

    const lastSys = resSys.elements[resSys.elements.length - 1];
    const val1 = lastSys.elements[1].value; // y1 at pi/2 approx 1
    const val2 = lastSys.elements[2].value; // y2 at pi/2 approx 0

    assert(Math.abs(val1 - 1.0) < 1e-3, `System y1: ${val1}`);
    assert(Math.abs(val2 - 0.0) < 1e-3, `System y2: ${val2}`);
    console.log("rk4 passed");
}

function testInv() {
    console.log("Testing inv optimization...");
    // 2x2
    const a = new Sym('a');
    const b = new Sym('b');
    const c = new Sym('c');
    const d = new Sym('d');
    const mat2 = new Vec([new Vec([a, b]), new Vec([c, d])]);
    const inv2 = cas.evaluate(new Call('inv', [mat2]));

    // Check element 0,0: d / (ad-bc)
    const el00 = inv2.elements[0].elements[0];
    console.log("Inv 2x2 [0,0]:", el00.toString());
    // Likely Mul(Div(1, Sub(Mul(a,d), Mul(b,c))), d) or similar structure
    assert(el00.toString().includes("/"), "Should contain division");
    assert(el00.toString().includes("a * d") || el00.toString().includes("d * a"), "Determinant term present");

    // 3x3 Symbolic
    const x = new Sym('x');
    const mat3 = new Vec([
        new Vec([new Num(1), new Num(0), new Num(0)]),
        new Vec([new Num(0), x, new Num(0)]),
        new Vec([new Num(0), new Num(0), new Num(1)])
    ]);
    const inv3 = cas.evaluate(new Call('inv', [mat3]));
    // Should be diag(1, 1/x, 1)
    const mid = inv3.elements[1].elements[1];
    console.log("Inv 3x3 [1,1]:", mid.toString());
    assert(mid.toString() === "(1 / x)", `Middle element should be 1/x, got ${mid.toString()}`);
    console.log("inv passed");
}

function testExp() {
    console.log("Testing exp simplification...");
    const x = new Sym('x');
    const y = new Sym('y');
    // exp(x) * exp(y)
    const mul = new Mul(new Call('exp', [x]), new Call('exp', [y]));
    const simple = mul.simplify();

    console.log("Exp combine:", simple.toString());
    assert(simple.funcName === 'exp', "Result is exp");
    assert(simple.args[0].toString() === "(x + y)" || simple.args[0].toString() === "(y + x)", "Argument is sum");

    // exp(x) * exp(-x)
    const mul2 = new Mul(new Call('exp', [x]), new Call('exp', [new Mul(new Num(-1), x)]));
    const simple2 = mul2.simplify();
    console.log("Exp cancel:", simple2.toString());
    assert(simple2.toString() === "1", "Should simplify to 1");
    console.log("exp passed");
}

try {
    testRK4();
    testInv();
    testExp();
    console.log("All tests passed");
} catch(e) {
    console.error(e);
    process.exit(1);
}
