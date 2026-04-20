// Test file for new features: bounded optimization + adaptive nIntegrate
var fs = require('fs');
var path = require('path');
var vm = require('vm');

// Use absolute path to the repo
var BASE = '/home/openclaw/.openclaw/workspace/repos/bookish-waffle';

var sandbox = { console: console, Math: Math, Number: Number, parseFloat: parseFloat, parseInt: parseInt };
vm.createContext(sandbox);

var exprCode = fs.readFileSync(path.join(BASE, 'js/expression.js'), 'utf8');
var parserCode = fs.readFileSync(path.join(BASE, 'js/parser.js'), 'utf8');
var casCode = fs.readFileSync(path.join(BASE, 'js/cas.js'), 'utf8');
var exposeCode = 'globalThis.Num=Num;globalThis.Vec=Vec;globalThis.Parser=Parser;globalThis.CAS=CAS;globalThis.Lexer=Lexer;';

var fullCode = [exprCode, parserCode, casCode, exposeCode].join('\n');
vm.runInContext(fullCode, sandbox);

var cas = new sandbox.CAS();
var parser = function(text) { return new sandbox.Parser(new sandbox.Lexer(text)); };

var passed = 0;
var failed = 0;

function getVal(result) {
    if (!result) return result;
    var cname = result.constructor && result.constructor.name;
    if (cname === 'Num') return result.evaluateNumeric();
    if (cname === 'Vec') return result.elements.map(function(e) { return e.evaluateNumeric(); });
    if (typeof result.evaluateNumeric === 'function') return result.evaluateNumeric();
    return result.toString();
}

function test(name, input, check) {
    try {
        var result = cas.evaluate(parser(input));
        var val = getVal(result);
        if (check(val)) {
            console.log('[PASS] ' + name + ': ' + input + ' => ' + JSON.stringify(val));
            passed++;
        } else {
            console.log('[FAIL] ' + name + ': ' + input + ' => ' + JSON.stringify(val) + ' (check failed)');
            failed++;
        }
    } catch(e) {
        console.log('[FAIL] ' + name + ': ' + input + ' => ERROR: ' + e.message);
        failed++;
    }
}

console.log('=== Adaptive nIntegrate ===');
test('nIntegrate x^2 [0,1]', 'nIntegrate(x^2, x, 0, 1)', function(v) {
    return Math.abs(v - 1/3) < 1e-6;
});
test('nIntegrate sin(x) [0,pi]', 'nIntegrate(sin(x), x, 0, pi)', function(v) {
    return Math.abs(v - 2) < 1e-6;
});
test('nIntegrate x^2 n=200', 'nIntegrate(x^2, x, 0, 1, 200)', function(v) {
    return Math.abs(v - 1/3) < 1e-6;
});
test('nIntegrate sin(10x) [0,pi]', 'nIntegrate(sin(10*x), x, 0, pi, 200)', function(v) {
    return Math.abs(v) < 0.01;
});

console.log('\n=== Bounded Golden Section Search ===');
test('fmin x^2 bounded [-2,2]', 'minimize(x^2, x, -2, 2)', function(v) {
    return Array.isArray(v) && v.length === 2 && Math.abs(v[0]) < 0.01 && Math.abs(v[1]) < 0.01;
});
test('fmax sin(x) bounded [0,pi]', 'maximize(sin(x), x, 0, pi)', function(v) {
    return Array.isArray(v) && Math.abs(v[0] - Math.PI/2) < 0.01 && Math.abs(v[1] - 1) < 0.01;
});
test('fmin x^2-4x+7 bounded [0,5]', 'minimize(x^2-4*x+7, x, 0, 5)', function(v) {
    return Array.isArray(v) && Math.abs(v[0] - 2) < 0.01 && Math.abs(v[1] - 3) < 0.01;
});

console.log('\n=== Regression: existing features still work ===');
test('fsolve works', 'fsolve(cos(x)-x, x, 0.5)', function(v) {
    return Math.abs(v - 0.739085) < 0.001;
});
test('symbolic minimize works', 'minimize(x^2, x)', function(v) {
    return Math.abs(v - 0) < 0.001;
});

console.log('\n=== Results: ' + passed + ' passed, ' + failed + ' failed ===');
if (failed > 0) process.exit(1);
