// Need to load expression first to populate globalThis/module
require('../js/expression.js');
const { CAS } = require('../js/cas.js');

// Now globals like Sym, Num should be available via globalThis if js/expression.js set them,
// but in node require, they are also returned.
// Let's access them from globalThis for safety as CAS relies on them.

const { Sym, Num, Add, Pow, Call, Vec } = globalThis;

const cas = new CAS();

// Helper to check output
const check = (desc, actual, expected) => {
    // Basic approximate check for numbers in string
    if (actual.toString() === expected.toString()) {
        console.log(`PASS: ${desc}`);
    } else {
        console.error(`FAIL: ${desc}`);
        console.error(`  Expected: ${expected}`);
        console.error(`  Actual:   ${actual}`);
    }
};

console.log("--- Testing codegen ---");
// Setup
const x = new Sym('x');
const expr1 = new Add(new Pow(x, new Num(2)), new Call('sin', [x]));

const py = cas._codegen(expr1, 'python').text;
check('codegen python x^2+sin(x)', py, '(x ** 2 + np.sin(x))');

const js = cas._codegen(expr1, 'js').text;
check('codegen js x^2+sin(x)', js, '(x ** 2 + Math.sin(x))');

const expr2 = new Call('ln', [new Sym('e')]);
const c = cas._codegen(expr2, 'c').text;
check('codegen c ln(e)', c, 'log(M_E)');

console.log("\n--- Testing convert ---");
const val1 = cas._convert(new Num(1), 'm', 'cm');
check('convert 1 m to cm', val1, '100');

const val2 = cas._convert(new Num(1), 'kg', 'g');
check('convert 1 kg to g', val2, '1000');

const val3 = cas._convert(new Num(0), 'C', 'K');
check('convert 0 C to K', val3, '273.15');

console.log("\n--- Testing MSE/MAE ---");
const list = new Vec([new Num(1), new Num(2), new Num(3)]);
const mse = cas._mse(list);
// 14/3 approx 4.666666666667
const mseVal = mse.evaluateNumeric();
if (Math.abs(mseVal - 14/3) < 1e-9) console.log("PASS: mse [1,2,3]");
else console.error(`FAIL: mse [1,2,3] got ${mseVal}`);

const list2 = new Vec([new Num(-1), new Num(-2)]);
const mae = cas._mae(list2);
const maeVal = mae.evaluateNumeric();
if (Math.abs(maeVal - 1.5) < 1e-9) console.log("PASS: mae [-1,-2]");
else console.error(`FAIL: mae [-1,-2] got ${maeVal}`);

console.log("\nTests Complete");
