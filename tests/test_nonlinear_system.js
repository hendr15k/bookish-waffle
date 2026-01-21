
const { CAS } = require('../js/cas');
const expr = require('../js/expression');
Object.assign(global, expr);

const cas = new CAS();

// Test 1: Non-linear system solving
// y = x^2, y = x + 2
// Expected: (2, 4) and (-1, 1)
const eq1 = new Eq(new Sym('y'), new Pow(new Sym('x'), new Num(2)));
const eq2 = new Eq(new Sym('y'), new Add(new Sym('x'), new Num(2)));
const vars = new Vec([new Sym('x'), new Sym('y')]);
const eqs = new Vec([eq1, eq2]);

try {
    const res = cas.evaluate(new Call('solve', [eqs, vars]));
    const resStr = res.toString();
    console.log("Result:", resStr);

    // Check for correct values (order might vary)
    if (!resStr.includes("x = 2") || !resStr.includes("y = 4")) throw new Error("Missing solution (2, 4)");
    if (!resStr.includes("x = -1") || !resStr.includes("y = 1")) throw new Error("Missing solution (-1, 1)");

} catch (e) {
    console.error("Test Failed:", e.message);
    process.exit(1);
}

console.log("Test Passed");
