
const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Assignment, Eq, Vec } = require('../js/expression');
const { CAS } = require('../js/cas');

// Make classes global so CAS can use them (if needed by eval logic relying on global)
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
global.Assignment = Assignment;
global.Eq = Eq;
global.Vec = Vec;

const cas = new CAS();

function test(name, exprStr, expectedStr) {
    try {
        const result = cas.evaluate(exprStr);
        const resStr = result.toString();
        if (resStr === expectedStr) {
            console.log(`[PASS] ${name}`);
        } else {
            console.error(`[FAIL] ${name}: Expected ${expectedStr}, got ${resStr}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${name}: Error ${e.message}`);
    }
}

// 1. Sigmoid, ReLU, Softplus
console.log("--- ML Functions ---");
// sigmoid(0) -> 0.5
const sig0 = new Call('sigmoid', [new Num(0)]);
const resSig0 = sig0.evaluateNumeric();
if (Math.abs(resSig0 - 0.5) < 1e-9) console.log("[PASS] sigmoid(0) numeric");
else console.error(`[FAIL] sigmoid(0): ${resSig0}`);

// diff(sigmoid(x), x) -> sigmoid(x)*(1-sigmoid(x))
const x = new Sym('x');
const diffSig = new Call('diff', [new Call('sigmoid', [x]), x]);
const resDiffSig = cas.evaluate(diffSig);
if (resDiffSig.toString().includes("sigmoid(x)")) console.log("[PASS] diff(sigmoid)");
else console.error(`[FAIL] diff(sigmoid): ${resDiffSig.toString()}`);

// relu(-5) -> 0
const reluNeg = new Call('relu', [new Num(-5)]);
if (reluNeg.evaluateNumeric() === 0) console.log("[PASS] relu(-5)");

// integrate(relu(x), x) -> 0.5 * x * relu(x)
const intRelu = new Call('integrate', [new Call('relu', [x]), x]);
const resIntRelu = cas.evaluate(intRelu);
if (resIntRelu.toString().includes("relu(x)")) console.log("[PASS] integrate(relu)");

// softplus(0) -> ln(2) approx 0.693
const sp0 = new Call('softplus', [new Num(0)]);
if (Math.abs(sp0.evaluateNumeric() - Math.log(2)) < 1e-9) console.log("[PASS] softplus(0)");

// integrate(softplus(x), x) -> -Li2(-e^x) = -polylog(2, -exp(x))
const intSP = new Call('integrate', [new Call('softplus', [x]), x]);
const resIntSP = cas.evaluate(intSP);
if (resIntSP.toString().includes("polylog") && resIntSP.toString().includes("exp")) console.log("[PASS] integrate(softplus)");
else console.error(`[FAIL] integrate(softplus): ${resIntSP.toString()}`);

// 2. Calculus Improvements
console.log("--- Calculus Improvements ---");
// diff(sech(x), x)
const diffSech = cas.evaluate(new Call('diff', [new Call('sech', [x]), x]));
if (diffSech.toString().includes("sech(x)") && diffSech.toString().includes("tanh(x)")) console.log("[PASS] diff(sech)");

// integrate(acoth(x), x)
const intAcoth = cas.evaluate(new Call('integrate', [new Call('acoth', [x]), x]));
if (intAcoth.toString().includes("acoth(x)") && intAcoth.toString().includes("ln")) console.log("[PASS] integrate(acoth)");

// 3. Algebra/Stats
console.log("--- Algebra/Stats ---");
// anti_commutator(x, y) -> x*y + y*x
const y = new Sym('y');
const ac = cas.evaluate(new Call('anti_commutator', [x, y]));
if (ac.toString().includes("x * y") && ac.toString().includes("y * x")) console.log("[PASS] anti_commutator");

// entropy([0.5, 0.5]) -> 1
const p = new Vec([new Num(0.5), new Num(0.5)]);
const ent = cas.evaluate(new Call('entropy', [p]));
if (Math.abs(ent.evaluateNumeric() - 1) < 1e-9) console.log("[PASS] entropy([0.5, 0.5])");

// kl_divergence([0.5, 0.5], [0.5, 0.5]) -> 0
const kl = cas.evaluate(new Call('kl_divergence', [p, p]));
if (Math.abs(kl.evaluateNumeric() - 0) < 1e-9) console.log("[PASS] kl_divergence identical");

// 4. Black Scholes
console.log("--- Black Scholes ---");
// S=100, K=100, T=1, r=0.05, sigma=0.2, type='call'
// Using online calc: ~10.45
const bs = new Call('black_scholes', [new Num(100), new Num(100), new Num(1), new Num(0.05), new Num(0.2), new Sym('call')]);
const resBS = cas.evaluate(bs);
const valBS = resBS.evaluateNumeric();
if (valBS > 10 && valBS < 11) console.log(`[PASS] black_scholes call: ${valBS}`);
else console.error(`[FAIL] black_scholes: ${valBS}`);
