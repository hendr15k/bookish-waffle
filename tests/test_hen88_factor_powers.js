const vm = require('vm');
const fs = require('fs');

const sandbox = { console, Math, Number, parseFloat, parseInt, Array, Object, isNaN, isFinite, JSON, String, Boolean, window: {} };
vm.createContext(sandbox);
vm.runInContext(fs.readFileSync('js/expression.js', 'utf8'), sandbox);
vm.runInContext(fs.readFileSync('js/parser.js', 'utf8'), sandbox);
vm.runInContext(fs.readFileSync('js/help.js', 'utf8'), sandbox);
vm.runInContext(fs.readFileSync('js/cas.js', 'utf8'), sandbox);

const cas = new sandbox.CAS();
const parse = (s) => new sandbox.Parser(new sandbox.Lexer(s)).parse();

let passed = 0, failed = 0;

function testFactor(expr, checkFn, description) {
    try {
        const result = cas.evaluate(parse('factor(' + expr + ')'));
        const str = result.toString();
        if (checkFn(result, str)) {
            console.log('[PASS]', description);
            passed++;
        } else {
            console.error('[FAIL]', description);
            console.error('  Expression:', expr);
            console.error('  Got:', str);
            failed++;
        }
    } catch (e) {
        console.error('[FAIL]', description, ':', e.message);
        failed++;
    }
}

// HEN-88: repeated factors not combined into powers
console.log('=== HEN-88: _factor repeated factors combining ===');

// The original bug case: flat Mul chain instead of power
testFactor('2*x^2+4*x+2',
    (r, s) => s.includes('^2') && s.includes('(x + 1)'),
    'factor(2*x^2+4*x+2) => 2*(x+1)^2');

testFactor('3*x^2+6*x+3',
    (r, s) => s.includes('^2') && s.includes('(x + 1)'),
    'factor(3*x^2+6*x+3) => 3*(x+1)^2');

// No flat Mul chain with repeated factors (the core bug)
testFactor('x^2+2*x+1',
    (r, s) => s.includes('^2') && s.includes('(x + 1)'),
    'factor(x^2+2*x+1) => (x+1)^2');

testFactor('x^2-2*x+1',
    (r, s) => s.includes('^2') && s.includes('(x - 1)'),
    'factor(x^2-2*x+1) => (x-1)^2');

// Cubic with triple root
testFactor('x^3-3*x^2+3*x-1',
    (r, s) => s.includes('^3') && s.includes('(x - 1)'),
    'factor(x^3-3*x^2+3*x-1) => (x-1)^3');

testFactor('x^3+6*x^2+12*x+8',
    (r, s) => s.includes('^3') && s.includes('(x + 2)'),
    'factor(x^3+6*x^2+12*x+8) => (x+2)^3');

// Quartic with quadruple root
testFactor('x^4-4*x^3+6*x^2-4*x+1',
    (r, s) => s.includes('^4') && s.includes('(x - 1)'),
    'factor(x^4-4*x^3+6*x^2-4*x+1) => (x-1)^4');

// Content + repeated: content stays flat, factors combine
testFactor('6*x^2+12*x+6',
    (r, s) => s.includes('^2') && s.includes('6'),
    'factor(6*x^2+12*x+6) => 6*(x+1)^2');

// No flat Mul chain (bug check)
testFactor('2*x^2+4*x+2',
    (r, s) => !s.match(/\(x \+ 1\)\) \* \(x \+ 1\)/),
    'No flat Mul chain for (x+1)*(x+1)');

// Distinct factors should NOT be combined
testFactor('x^2-5*x+6',
    (r, s) => !s.includes('^') && s.includes('x - 2') && s.includes('x - 3'),
    'factor(x^2-5*x+6) => (x-2)(x-3) distinct factors');

// x^4-2x^2+1 => (x-1)^2*(x+1)^2
testFactor('x^4-2*x^2+1',
    (r, s) => s.includes('(x - 1))^2') && s.includes('(x + 1))^2'),
    'factor(x^4-2*x^2+1) => (x-1)^2*(x+1)^2');

console.log('\nResults: ' + passed + ' passed, ' + failed + ' failed');
process.exit(failed > 0 ? 1 : 0);
