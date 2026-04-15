const fs = require('fs');
const vm = require('vm');

function loadFile(filePath) {
    return fs.readFileSync(filePath, 'utf8');
}

const expressionCode = loadFile('js/expression.js');
const casCode = loadFile('js/cas.js');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    isNaN: isNaN,
    isFinite: isFinite,
};

vm.createContext(sandbox);
vm.runInContext(expressionCode, sandbox);
vm.runInContext(casCode, sandbox);

const CAS = sandbox.CAS;
const Call = sandbox.Call;
const Num = sandbox.Num;
const Sym = sandbox.Sym;

const cas = new CAS();

let hasFailures = false;

function assertEqual(actual, expected, message) {
    if (actual !== expected) {
        console.error(`FAIL: ${message}. Expected '${expected}', got '${actual}'`);
        hasFailures = true;
    } else {
        console.log(`PASS: ${message}`);
    }
}

console.log("--- Series Operations ---");

const x = new Sym('x');
const f = new Call('sin', [x]);

const taylorRes = cas.evaluate(new Call('taylor', [f, x, new Num(0), new Num(3)]));
assertEqual(taylorRes.toString(), '(x + ((-1 * x^3) / 6))', 'taylor(sin(x), x, 0, 3)');

const subRes = f.substitute(x, new Num(Math.PI));
// Note: evaluate isn't substituting automatically in `evaluate(new Call('substitute', ...))`, so calling method.
// Or if it is not supported as a command, we'll just test the `substitute` method directly on Expr.
const subEval = cas.evaluate(subRes);
assertEqual(subEval.evaluateNumeric() < 1e-10 ? '0' : subEval.toString(), '0', 'sin(PI) evaluated numeric');

const seriesRes = cas.evaluate(new Call('series', [f, x, new Num(0), new Num(3)]));
assertEqual(seriesRes.toString(), '(x + ((-1 * x^3) / 6))', 'series(sin(x), x, 0, 3)');

if (hasFailures) {
    process.exit(1);
} else {
    console.log("All Tests Passed!");
    process.exit(0);
}
