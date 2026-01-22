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
    // Add any other globals
};

vm.createContext(sandbox);
vm.runInContext(expressionCode, sandbox);
vm.runInContext(casCode, sandbox);

// Extract classes from sandbox
const CAS = sandbox.CAS;
const Call = sandbox.Call;
const Num = sandbox.Num;
const Sym = sandbox.Sym;
const Div = sandbox.Div;
const Sub = sandbox.Sub;

const cas = new CAS();

console.log("--- Special Functions Numeric ---");
const x = new Num(2.0);
const si = new Call('Si', [x]).evaluateNumeric();
console.log('Si(2) = ' + si + ' (Expect ~1.605)');
if (Math.abs(si - 1.6054129768) > 1e-5) throw new Error("Si failed");

const ci = new Call('Ci', [x]).evaluateNumeric();
console.log('Ci(2) = ' + ci + ' (Expect ~0.422)');
if (Math.abs(ci - 0.42298082877) > 1e-5) throw new Error("Ci failed");

const airyAi = new Call('airyAi', [new Num(1.0)]).evaluateNumeric();
console.log('Ai(1) = ' + airyAi + ' (Expect ~0.135)');
if (Math.abs(airyAi - 0.1352924163) > 1e-5) throw new Error("Ai failed");

const airyBi = new Call('airyBi', [new Num(1.0)]).evaluateNumeric();
console.log('Bi(1) = ' + airyBi + ' (Expect ~1.207)');
if (Math.abs(airyBi - 1.2074235949) > 1e-5) throw new Error("Bi failed");

console.log("--- Integration Special ---");
const varX = new Sym('x');
// integrate(sin(x)/x, x) -> Si(x)
// Note: CAS.evaluate returns Expr objects
const int1 = cas.evaluate(new Call('integrate', [new Div(new Call('sin', [varX]), varX), varX]));
console.log('int(sin(x)/x) = ' + int1.toString());
if (int1.toString() !== 'Si(x)') throw new Error("int(sin(x)/x) failed");

const int2 = cas.evaluate(new Call('integrate', [new Div(new Call('cos', [varX]), varX), varX]));
console.log('int(cos(x)/x) = ' + int2.toString());
if (int2.toString() !== 'Ci(x)') throw new Error("int(cos(x)/x) failed");

const int3 = cas.evaluate(new Call('integrate', [new Div(new Call('exp', [varX]), varX), varX]));
console.log('int(exp(x)/x) = ' + int3.toString());
if (int3.toString() !== 'Ei(x)') throw new Error("int(exp(x)/x) failed");

const int4 = cas.evaluate(new Call('integrate', [new Div(new Num(1), new Call('ln', [varX])), varX]));
console.log('int(1/ln(x)) = ' + int4.toString());
if (int4.toString() !== 'Li(x)') throw new Error("int(1/ln(x)) failed");

console.log("--- Limits ---");
const lim1 = cas.evaluate(new Call('limit', [new Div(new Call('sin', [varX]), varX), varX, new Num(0)]));
console.log('limit(sin(x)/x, x, 0) = ' + lim1.toString());
if (lim1.toString() !== '1') throw new Error("limit(sin(x)/x) failed");

// (e^x - 1)/x
const expr5 = new Div(new Sub(new Call('exp', [varX]), new Num(1)), varX);
const lim2 = cas.evaluate(new Call('limit', [expr5, varX, new Num(0)]));
console.log('limit((e^x-1)/x, x, 0) = ' + lim2.toString());
if (lim2.toString() !== '1') throw new Error("limit((e^x-1)/x) failed");

console.log("All Tests Passed!");
