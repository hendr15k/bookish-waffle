const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Eq, Sub, Add, Pow } = require('./js/expression.js');
const cas = new CAS();
const x = new Sym('x');
function cmd(name, ...args) { return cas.evaluate(new Call(name, args)); }

function parseExpr(str) {
    return cas.evaluate(str);
}

let passed = 0;
let failed = 0;

function test(name, fn) {
    try {
        const r = fn();
        const str = r && r.toString ? r.toString() : String(r);
        console.log(`  PASS: ${name} => ${str}`);
        passed++;
    } catch(e) {
        console.log(`  FAIL: ${name} => ${e.message}`);
        failed++;
    }
}

console.log('\n=== REGRESSION TESTS: HEN-75 Bug Fixes ===\n');

// BUG #2: abs(x) = -3 → empty set
console.log('--- BUG #2: abs equation with negative RHS ---');
test('solve(abs(x) = -3) = {}', () => {
    const eq = new Eq(new Call('abs', [x]), new Num(-3));
    const r = cmd('solve', eq, x);
    if (r.toString() !== '{}') throw new Error(`Expected {}, got ${r.toString()}`);
    return r;
});

test('solve(abs(x) = -5) = {}', () => {
    const eq = new Eq(new Call('abs', [x]), new Num(-5));
    const r = cmd('solve', eq, x);
    if (r.toString() !== '{}') throw new Error(`Expected {}, got ${r.toString()}`);
    return r;
});

test('solve(abs(x) = 0) = {0}', () => {
    const eq = new Eq(new Call('abs', [x]), new Num(0));
    const r = cmd('solve', eq, x);
    if (r.toString() !== '{0}') throw new Error(`Expected {0}, got ${r.toString()}`);
    return r;
});

test('solve(abs(x) = 2) = {-2, 2}', () => {
    const eq = new Eq(new Call('abs', [x]), new Num(2));
    const r = cmd('solve', eq, x);
    return r;
});

// BUG #1: High-degree polynomial roots
console.log('\n--- BUG #1: High-degree polynomial roots ---');
test('solve(x^7 - 1) = 7 roots', () => {
    const poly = new Sub(new Pow(x, new Num(7)), new Num(1));
    const sol = cmd('solve', poly, x);
    if (!sol.args || sol.args.length !== 7) throw new Error(`Expected 7 roots, got ${sol.args ? sol.args.length : '?'}`);
    return sol;
});

test('solve(x^5 - 1) = 5 roots', () => {
    const poly = new Sub(new Pow(x, new Num(5)), new Num(1));
    const sol = cmd('solve', poly, x);
    if (!sol.args || sol.args.length !== 5) throw new Error(`Expected 5 roots, got ${sol.args ? sol.args.length : '?'}`);
    return sol;
});

test('solve(x^5 + 1) = 5 roots', () => {
    const poly = new Add(new Pow(x, new Num(5)), new Num(1));
    const sol = cmd('solve', poly, x);
    if (!sol.args || sol.args.length !== 5) throw new Error(`Expected 5 roots, got ${sol.args ? sol.args.length : '?'}`);
    return sol;
});

test('solve(x^4 + x^2 + 1) = 4 roots (no crash)', () => {
    const poly = new Add(new Add(new Pow(x, new Num(4)), new Pow(x, new Num(2))), new Num(1));
    const sol = cmd('solve', poly, x);
    if (!sol.args || sol.args.length !== 4) throw new Error(`Expected 4 roots, got ${sol.args ? sol.args.length : '?'}`);
    return sol;
});

test('solve(x^6 - 1) = 6 roots', () => {
    const poly = new Sub(new Pow(x, new Num(6)), new Num(1));
    const sol = cmd('solve', poly, x);
    if (!sol.args || sol.args.length !== 6) throw new Error(`Expected 6 roots, got ${sol.args ? sol.args.length : '?'}`);
    return sol;
});

console.log('\n--- SANITY: Basic quadratic still works ---');
test('solve(x^2 - 4) = {-2, 2}', () => {
    const poly = new Sub(new Pow(x, new Num(2)), new Num(4));
    return cmd('solve', poly, x);
});

console.log(`\n=== RESULTS: ${passed} passed, ${failed} failed ===\n`);
process.exit(failed > 0 ? 1 : 0);
