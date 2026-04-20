const vm = require('vm');
const fs = require('fs');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    Array: Array,
    Object: Object,
    isNaN: isNaN,
    isFinite: isFinite,
    JSON: JSON,
    String: String,
    Boolean: Boolean,
};

const expr = fs.readFileSync('js/expression.js', 'utf8');
const parser = fs.readFileSync('js/parser.js', 'utf8');
const cas = fs.readFileSync('js/cas.js', 'utf8');
const help = fs.readFileSync('js/help.js', 'utf8');

vm.createContext(sandbox);
vm.runInContext(expr, sandbox);
vm.runInContext(parser, sandbox);
vm.runInContext(help, sandbox);
vm.runInContext(cas, sandbox);

const myCAS = new sandbox.CAS();

const tests = [
    ['1+1', '2'],
    ['nIntegrate(x^2, x, 0, 1)', '1/3 ≈ 0.333'],
    ['fsolve(cos(x)-x, x, 0.5)', '≈ 0.739'],
    ['solve(x^2-4, x)', '{-2, 2}'],
    ['minimize(x^2, x)', '0'],
];

for (const [input, expected] of tests) {
    try {
        const result = myCAS.evaluate(sandbox.Parser.parse(input));
        console.log(`✓ ${input} = ${result.toString()}`);
    } catch(e) {
        console.log(`✗ ${input}: ${e.message}`);
    }
}
