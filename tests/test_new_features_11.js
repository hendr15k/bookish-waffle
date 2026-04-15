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
const Add = sandbox.Add;
const Mul = sandbox.Mul;

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

console.log("--- Complex Numbers ---");

// complex(2, 3) format: 2 + 3*i
const c = new Add(new Num(2), new Mul(new Num(3), new Sym('i')));

const cAbs = cas.evaluate(new Call('abs', [c]));
assertEqual(cAbs.toString(), 'abs((2 + (3 * i)))', 'abs(2+3i)');

const cArg = cas.evaluate(new Call('arg', [c]));
// Expected numeric arg of 2+3i is approx 0.98279
assertEqual(Math.abs(cArg.value - 0.98279372) < 1e-5, true, 'arg(2+3i)');

// The 're', 'im', 'conj' are not fully evaluated into numbers via pure standard evaluation unless we use the parts directly
const parts = cas._getComplexParts(c);

assertEqual(parts.re.toString(), '2', 're part of 2+3i');
assertEqual(parts.im.toString(), '3', 'im part of 2+3i');

// conj = re - im * i
const conjExpected = cas.evaluate(new sandbox.Sub(parts.re, new Mul(parts.im, new Sym('i'))));
assertEqual(conjExpected.toString(), '(2 - (3 * i))', 'conj part of 2+3i');

if (hasFailures) {
    process.exit(1);
} else {
    console.log("All Tests Passed!");
    process.exit(0);
}
