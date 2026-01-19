
const assert = require('assert');
require('../js/expression'); // Load Expr classes globally
const { CAS } = require('../js/cas');

console.log("--- Testing Regression & Stats ---");

const cas = new CAS();

function testPowerRegression() {
    console.log("Testing Power Regression y = A x^B...");
    // Data: y = 2 * x^3
    // ln y = ln 2 + 3 ln x
    const data = new Vec([
        new Vec([new Num(1), new Num(2)]),
        new Vec([new Num(2), new Num(16)]),
        new Vec([new Num(3), new Num(54)])
    ]);

    const reg = cas._powerRegression(data);
    console.log("Result:", reg.toString());

    // Check form: 2 * x^3
    // Due to floating point, coefficients might be slightly off
    // 2.000... * x^3.000...

    if (reg instanceof Mul && reg.left instanceof Num && reg.right instanceof Pow) {
        const A = reg.left.value;
        const B = reg.right.right.value;

        if (Math.abs(A - 2) < 1e-5 && Math.abs(B - 3) < 1e-5) {
            console.log("[PASS] Power Regression");
        } else {
            console.error("[FAIL] Power Regression coeffs mismatch", A, B);
        }
    } else {
        console.error("[FAIL] Power Regression struct mismatch", reg.toString());
    }
}

function testLogRegression() {
    console.log("Testing Log Regression y = A + B ln(x)...");
    // Data: y = 1 + 2 ln(x)
    const data = new Vec([
        new Vec([new Num(1), new Num(1)]), // ln(1)=0 -> y=1
        new Vec([new Num(Math.E), new Num(3)]), // ln(e)=1 -> y=1+2=3
        new Vec([new Num(Math.E*Math.E), new Num(5)]) // ln(e^2)=2 -> y=1+4=5
    ]);

    const reg = cas._logRegression(data);
    console.log("Result:", reg.toString());

    // Check form: 1 + 2 * ln(x)
    if (reg instanceof Add) {
        // Order in Add might vary? simplified to Num + Mul
        let A, B;
        if (reg.left instanceof Num) {
            A = reg.left.value;
            // right could be Mul(2, ln(x)) or ln(x^2)
            if (reg.right instanceof Mul && reg.right.left instanceof Num) {
                B = reg.right.left.value;
            } else if (reg.right instanceof Call && reg.right.funcName === 'ln') {
                // Check argument x^2
                const arg = reg.right.args[0];
                if (arg instanceof Pow && arg.right instanceof Num) {
                     B = arg.right.value;
                }
            }
        }

        if (Math.abs(A - 1) < 1e-5 && Math.abs(B - 2) < 1e-5) {
            console.log("[PASS] Log Regression");
        } else {
            console.error("[FAIL] Log Regression coeffs mismatch", A, B);
        }
    } else {
        console.error("[FAIL] Log Regression struct mismatch", reg.toString());
    }
}

function testStudentT() {
    console.log("Testing Student-T PDF...");
    // t-distribution with df=1 (Cauchy)
    // f(0) = 1/pi approx 0.3183
    const val = cas._studentTPDF(new Num(0), new Num(1));
    if (Math.abs(val.value - 1/Math.PI) < 1e-5) {
        console.log("[PASS] Student-T PDF");
    } else {
        console.error("[FAIL] Student-T PDF", val.value);
    }
}

testPowerRegression();
testLogRegression();
testStudentT();
