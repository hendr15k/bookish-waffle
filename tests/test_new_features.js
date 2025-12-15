const vm = require('vm');
const fs = require('fs');
const path = require('path');

// Load files
const expressionCode = fs.readFileSync(path.join(__dirname, '../js/expression.js'), 'utf8');
const parserCode = fs.readFileSync(path.join(__dirname, '../js/parser.js'), 'utf8');
const casCode = fs.readFileSync(path.join(__dirname, '../js/cas.js'), 'utf8');

// Sandbox
const sandbox = {
    console: console,
    Math: Math,
    Number: Number
};
vm.createContext(sandbox);

vm.runInContext(expressionCode, sandbox);
vm.runInContext(parserCode, sandbox);
vm.runInContext(casCode, sandbox);

// Tests
const tests = [
    // Bernoulli
    { cmd: 'bernoulli(0)', expected: '1' },
    { cmd: 'bernoulli(1)', check: (res) => Math.abs(parseFloat(res.evaluateNumeric()) + 0.5) < 1e-9 },
    { cmd: 'bernoulli(2)', check: (res) => Math.abs(parseFloat(res.evaluateNumeric()) - 1/6) < 1e-9 },

    // Harmonic
    { cmd: 'harmonic(1)', expected: '1' },
    { cmd: 'harmonic(2)', check: (res) => Math.abs(parseFloat(res.evaluateNumeric()) - 1.5) < 1e-9 },
    { cmd: 'harmonic(3)', check: (res) => Math.abs(parseFloat(res.evaluateNumeric()) - 11/6) < 1e-9 },

    // Analytic Geometry
    { cmd: 'lineEquation([0,0], [1,1])', expected: 'y = x' },
    // Circle: (x-0)^2 + (y-0)^2 = 1^2 => x^2 + y^2 = 1. Wait, CAS prints as ((x^2) + (y^2)) = 1 ?
    { cmd: 'circleEquation([0,0], 1)', check: (res) => res.toString().replace(/\s/g,'').replace(/\(/g,'').replace(/\)/g,'') === 'x^2+y^2=1' }
];

console.log('Running tests...');
let failed = 0;

tests.forEach(t => {
    try {
        const lexer = new sandbox.Lexer(t.cmd);
        const parser = new sandbox.Parser(lexer);
        const tree = parser.parse();
        const cas = new sandbox.CAS();
        const result = cas.evaluate(tree);
        const resStr = result.toString();

        if (t.check) {
            if (t.check(result)) {
                console.log(`PASS: ${t.cmd}`);
            } else {
                console.error(`FAIL: ${t.cmd}. Got ${resStr}`);
                failed++;
            }
        } else {
            if (resStr === t.expected) {
                console.log(`PASS: ${t.cmd}`);
            } else {
                console.error(`FAIL: ${t.cmd}. Expected ${t.expected}, got ${resStr}`);
                failed++;
            }
        }
    } catch(e) {
        console.error(`ERROR: ${t.cmd}: ${e.message}`);
        failed++;
    }
});

if (failed > 0) process.exit(1);
console.log('All tests passed');
