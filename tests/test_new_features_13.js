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

console.log("--- Set Operations ---");

const s1 = new Vec([new Num(1), new Num(2), new Num(3)]);
const s2 = new Vec([new Num(2), new Num(3), new Num(4)]);

const unionRes = cas.evaluate(new Call('union', [s1, s2]));
assertEqual(unionRes.toString(), '[1, 2, 3, 4]', 'union([1,2,3], [2,3,4])');

const intersectRes = cas.evaluate(new Call('intersect', [s1, s2]));
assertEqual(intersectRes.toString(), '[2, 3]', 'intersect([1,2,3], [2,3,4])');

const diffRes = cas.evaluate(new Call('setdiff', [s1, s2]));
assertEqual(diffRes.toString(), '[1]', 'setdiff([1,2,3], [2,3,4])');

const sizeRes = cas.evaluate(new Call('size', [s1]));
assertEqual(sizeRes.toString(), '3', 'size([1,2,3])');

const isSubsetRes1 = cas.evaluate(new Call('isSubset', [new Vec([new Num(2)]), s1]));
assertEqual(isSubsetRes1.toString(), '1', 'isSubset([2], [1,2,3])');

const isSubsetRes2 = cas.evaluate(new Call('isSubset', [s2, s1]));
assertEqual(isSubsetRes2.toString(), '0', 'isSubset([2,3,4], [1,2,3])');

if (hasFailures) {
    process.exit(1);
} else {
    console.log("All Tests Passed!");
    process.exit(0);
}
