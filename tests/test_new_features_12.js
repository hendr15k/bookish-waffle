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
const And = sandbox.And;
const Or = sandbox.Or;
const Not = sandbox.Not;
const Xor = sandbox.Xor;
const Lt = sandbox.Lt;
const Le = sandbox.Le;
const Gt = sandbox.Gt;
const Ge = sandbox.Ge;
const BooleanEq = sandbox.BooleanEq;
const Neq = sandbox.Neq;

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

console.log("--- Boolean Logic ---");

const t = new Num(1);
const f = new Num(0);

// and(1, 0)
assertEqual(cas.evaluate(new And(t, f)).toString(), '0', 'and(1, 0) == 0');
// or(1, 0)
assertEqual(cas.evaluate(new Or(t, f)).toString(), '1', 'or(1, 0) == 1');
// not(1)
assertEqual(cas.evaluate(new Not(t)).toString(), '0', 'not(1) == 0');
// xor(1, 1)
assertEqual(cas.evaluate(new Xor(t, t)).toString(), '0', 'xor(1, 1) == 0');
// xor(1, 0)
assertEqual(cas.evaluate(new Xor(t, f)).toString(), '1', 'xor(1, 0) == 1');

console.log("--- Relational Operators ---");
const num3 = new Num(3);
const num5 = new Num(5);

assertEqual(cas.evaluate(new Lt(num3, num5)).toString(), '1', '3 < 5');
assertEqual(cas.evaluate(new Lt(num5, num3)).toString(), '0', '5 < 3');

assertEqual(cas.evaluate(new Gt(num5, num3)).toString(), '1', '5 > 3');
assertEqual(cas.evaluate(new Gt(num3, num5)).toString(), '0', '3 > 5');

assertEqual(cas.evaluate(new Le(num3, num3)).toString(), '1', '3 <= 3');
assertEqual(cas.evaluate(new Ge(num5, num5)).toString(), '1', '5 >= 5');

assertEqual(cas.evaluate(new BooleanEq(num3, num5)).toString(), '0', '3 == 5');
assertEqual(cas.evaluate(new BooleanEq(num3, num3)).toString(), '1', '3 == 3');

assertEqual(cas.evaluate(new Neq(num3, num5)).toString(), '1', '3 != 5');
assertEqual(cas.evaluate(new Neq(num3, num3)).toString(), '0', '3 != 3');

if (hasFailures) {
    process.exit(1);
} else {
    console.log("All Tests Passed!");
    process.exit(0);
}
