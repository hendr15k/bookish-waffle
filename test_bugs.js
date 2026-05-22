const vm = require('vm');
const fs = require('fs');

const sandbox = {
    console: console, Math: Math, Number: Number, parseFloat: parseFloat, parseInt: parseInt,
    Array: Array, Object: Object, isNaN: isNaN, isFinite: isFinite, JSON: JSON,
    String: String, Boolean: Boolean, window: {},
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
const parse = (s) => new sandbox.Parser(new sandbox.Lexer(s)).parse();

// Test BUG #1: solve(x^5 - 1, x)
console.log('=== BUG #1: solve(x^5 - 1, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(x^5 - 1, x)'));
    console.log('Result:', result.toString());
    console.log('Type:', result.constructor.name);
    if (result.args) console.log('Root count:', result.args.length);
} catch(e) {
    console.log('ERROR:', e.message);
    console.log('Stack:', e.stack.split('\n').slice(0,8).join('\n'));
}

// Test BUG #1b: solve(x^4 + x^2 + 1, x)
console.log('\n=== BUG #1b: solve(x^4 + x^2 + 1, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(x^4 + x^2 + 1, x)'));
    console.log('Result:', result.toString());
    if (result.args) console.log('Root count:', result.args.length);
} catch(e) {
    console.log('ERROR:', e.message);
    console.log('Stack:', e.stack.split('\n').slice(0,8).join('\n'));
}

// Test BUG #2: abs(x) = -3
console.log('\n=== BUG #2: solve(abs(x) = -3, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(abs(x) = -3, x)'));
    console.log('Result:', result.toString());
} catch(e) {
    console.log('ERROR:', e.message);
}

// Test x^7 - 1
console.log('\n=== EXTRA: solve(x^7 - 1, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(x^7 - 1, x)'));
    console.log('Result:', result.toString());
    if (result.args) console.log('Root count:', result.args.length);
} catch(e) {
    console.log('ERROR:', e.message);
}

// Test solve(x^5+1, x)
console.log('\n=== EXTRA: solve(x^5+1, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(x^5+1, x)'));
    console.log('Result:', result.toString());
    if (result.args) console.log('Root count:', result.args.length);
} catch(e) {
    console.log('ERROR:', e.message);
}

// Test basic quadratic still works
console.log('\n=== SANITY: solve(x^2-4, x) ===');
try {
    const result = myCAS.evaluate(parse('solve(x^2-4, x)'));
    console.log('Result:', result.toString());
    if (result.args) console.log('Root count:', result.args.length);
} catch(e) {
    console.log('ERROR:', e.message);
}
