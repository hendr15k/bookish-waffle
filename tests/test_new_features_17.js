const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Assignment, Eq, Vec, FunctionDef, Block, toExpr, And, Or, Xor, Implies, Iff, Not, Mod, Neq, Lt, Gt, Le, Ge, At, BooleanEq, If, While, For, Return, Break, Continue } = require('../js/expression.js');
const { CAS } = require('../js/cas.js');

// Mock global environment for CAS if needed (though we import directly)
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.Vec = Vec;
global.Call = Call;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;
global.Eq = Eq;

const cas = new CAS();

console.log("Testing Factorial/Gamma Simplification...");

function checkNumeric(expr, varName, val, expectedVal) {
    const sub = expr.substitute(varName, new Num(val)).evaluateNumeric();
    if (Math.abs(sub - expectedVal) < 1e-9) return true;
    return false;
}

// n! / (n-1)! -> n
const n = new Sym('n');
const factN = new Call('factorial', [n]);
const factNM1 = new Call('factorial', [new Sub(n, new Num(1))]);
const div1 = new Div(factN, factNM1);
const res1 = cas.evaluate(div1);
if (res1.toString() === "n") console.log("PASS: n!/(n-1)!");
else console.log("FAIL: n!/(n-1)! -> " + res1.toString());

// n! / (n-3)! -> n*(n-1)*(n-2)
const factNM3 = new Call('factorial', [new Sub(n, new Num(3))]);
const div2 = new Div(factN, factNM3);
const res2 = cas.evaluate(div2); // Should be n*(n-1)*(n-2)
// Check at n=5. 5!/2! = 120/2 = 60.
// 5*4*3 = 60.
if (checkNumeric(res2, n, 5, 60)) console.log("PASS: n!/(n-3)! (Numeric Check)");
else console.log(`FAIL: n!/(n-3)! -> ${res2}`);

// gamma(n+2)/gamma(n) -> (n+1)n = n^2 + n
const gammaN2 = new Call('gamma', [new Add(n, new Num(2))]);
const gammaN = new Call('gamma', [n]);
const div3 = new Div(gammaN2, gammaN);
const res3 = cas.evaluate(div3);
// Check at n=4. gamma(6)/gamma(4) = 120/6 = 20.
// 4^2 + 4 = 20.
if (checkNumeric(res3, n, 4, 20)) console.log("PASS: gamma(n+2)/gamma(n) (Numeric Check)");
else console.log(`FAIL: gamma(n+2)/gamma(n) -> ${res3}`);


console.log("Testing Control Systems...");

// ss2tf
// A = [[0, 1], [-2, -3]]
// B = [[0], [1]]
// C = [[1, 0]]
// D = [[0]]
// H(s) = 1 / (s^2 + 3s + 2)
const A = new Vec([new Vec([new Num(0), new Num(1)]), new Vec([new Num(-2), new Num(-3)])]);
const B = new Vec([new Vec([new Num(0)]), new Vec([new Num(1)])]);
const C = new Vec([new Vec([new Num(1), new Num(0)])]);
const D = new Vec([new Vec([new Num(0)])]);

const tfCall = new Call('ss2tf', [A, B, C, D]);
const tfRes = cas.evaluate(tfCall);

// Let's check element [0][0]
if (tfRes instanceof Vec) {
    const el = tfRes.elements[0].elements[0];
    const s = new Sym('s');
    // Evaluate at s=1
    // A=[[0,1],[-2,-3]]. sI-A = [[1,-1],[2, 4]]. det = 4+2=6. inv = [[4,1],[-2,1]]/6.
    // C=[1,0]. C*inv = [4/6, 1/6].
    // * B=[0,1]. [4/6, 1/6]*[0,1] = 1/6.
    // + D=0. Result 1/6.
    const valRes = el.substitute(s, new Num(1)).evaluateNumeric();
    if (Math.abs(valRes - 1.0/6.0) < 1e-9) console.log("PASS: ss2tf numeric check at s=1");
    else console.log(`FAIL: ss2tf check. Got ${valRes}, expected ${1.0/6.0}`);
} else {
    console.log("FAIL: ss2tf did not return a matrix");
    console.log(tfRes.toString());
}

// step response
// H(s) = 1/s. Step response is t (ramp).
const sys = new Div(new Num(1), new Sym('s'));
const stepCall = new Call('step', [sys]);
const stepRes = cas.evaluate(stepCall);
if (stepRes.toString() === 't') console.log("PASS: step(1/s) -> t");
else console.log("FAIL: step(1/s) -> " + stepRes.toString());

// impulse response
// H(s) = 1/s. Impulse is 1 (step). ilaplace(1/s) = 1.
const impCall = new Call('impulse', [sys]);
const impRes = cas.evaluate(impCall);
if (impRes.toString() === '1') console.log("PASS: impulse(1/s) -> 1");
else console.log("FAIL: impulse(1/s) -> " + impRes.toString());
