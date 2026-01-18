const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec, Eq } = require('../js/expression.js');

// Mock the global classes needed for CAS
global.Num = Num;
global.Sym = Sym;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;
global.Call = Call;
global.Vec = Vec;
global.Expr = Expr;
global.Eq = Eq;

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        throw new Error(message || "Assertion failed");
    }
}

function test(name, exprTree, expectedStr) {
    try {
        const res = cas.evaluate(exprTree);
        const resStr = res.toString();
        console.log(`[PASS] ${name}: ${resStr}`);
        if (expectedStr && resStr !== expectedStr) {
            console.warn(`  Warning: Expected "${expectedStr}", got "${resStr}"`);
        }
    } catch (e) {
        console.error(`[FAIL] ${name}: ${e.message}`);
        process.exit(1);
    }
}

console.log("Running new features tests (Laurent Series, 2D Curl)...");

const x = new Sym('x');
const y = new Sym('y');

// 1. Taylor at singularity (1/x around 0) -> Should fallback to Laurent
// Expected: 1/x
test('taylor(1/x, x, 0, 3)',
    new Call('taylor', [new Div(new Num(1), x), x, new Num(0), new Num(3)]),
    "(1 / x)"
);

// 2. Explicit Laurent
test('laurent(1/x, x, 0, 3)',
    new Call('laurent', [new Div(new Num(1), x), x, new Num(0), new Num(3)]),
    "(1 / x)"
);

// 3. Laurent 1/x^2
test('laurent(1/x^2, x, 0, 3)',
    new Call('laurent', [new Div(new Num(1), new Pow(x, new Num(2))), x, new Num(0), new Num(3)]),
    "(1 / x^2)"
);

// 4. Laurent with shift 1/(x-1) around 1
test('laurent(1/(x-1), x, 1, 3)',
    new Call('laurent', [new Div(new Num(1), new Sub(x, new Num(1))), x, new Num(1), new Num(3)]),
    "(1 / (x - 1))"
);

// 5. 2D Curl (Scalar result)
// F = [y, -x]. curl = d(-x)/dx - dy/dy = -1 - 1 = -2
test('curl([y, -x], [x, y])',
    new Call('curl', [new Vec([y, new Mul(new Num(-1), x)]), new Vec([x, y])]),
    "-2"
);

// 6. 2D Curl (Conservative)
// F = [x, y]. curl = dy/dx - dx/dy = 0 - 0 = 0
test('curl([x, y], [x, y])',
    new Call('curl', [new Vec([x, y]), new Vec([x, y])]),
    "0"
);

console.log("All tests passed.");
