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
const Vec = sandbox.Vec;

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

console.log("--- Matrix Operations ---");

const m1 = new Vec([
    new Vec([new Num(1), new Num(2)]),
    new Vec([new Num(3), new Num(4)])
]);

const m2 = new Vec([
    new Vec([new Num(2), new Num(0)]),
    new Vec([new Num(0), new Num(3)])
]);

const detRes = cas.evaluate(new Call('det', [m1]));
assertEqual(detRes.toString(), '-2', 'det([[1,2],[3,4]])');

const invRes = cas.evaluate(new Call('inv', [m1]));
assertEqual(invRes.toString(), '[[-2, 1], [(3 / 2), (-1 / 2)]]', 'inv([[1,2],[3,4]])');

const rankRes = cas.evaluate(new Call('rank', [m1]));
assertEqual(rankRes.toString(), '2', 'rank([[1,2],[3,4]])');

const transRes = cas.evaluate(new Call('transpose', [m1]));
assertEqual(transRes.toString(), '[[1, 3], [2, 4]]', 'transpose([[1,2],[3,4]])');

const traceRes = cas.evaluate(new Call('trace', [m1]));
assertEqual(traceRes.toString(), '5', 'trace([[1,2],[3,4]])');

const eigenRes = cas.evaluate(new Call('eigenvals', [m2]));
assertEqual(eigenRes.toString(), '[2, 3]', 'eigenvals([[2,0],[0,3]])');

const diagRes = cas.evaluate(new Call('diagonalize', [m2]));
assertEqual(diagRes.toString(), '[[[1, 0], [0, 1]], [[2, 0], [0, 3]]]', 'diagonalize([[2,0],[0,3]])');

if (hasFailures) {
    process.exit(1);
} else {
    console.log("All Tests Passed!");
    process.exit(0);
}
