const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec, Eq, Assignment, FunctionDef } = require('../js/expression');
const { Lexer } = require('../js/parser');
const { Parser } = require('../js/parser');
const { CAS } = require('../js/cas');

// Mock globalThis.HELP_DATA
globalThis.HELP_DATA = { 'sin': {}, 'cos': {}, 'tan': {} };

const cas = new CAS();

function test(name, input, expectedStr) {
    try {
        const lexer = new Lexer(input);
        const parser = new Parser(lexer);
        const tree = parser.parse();
        const res = cas.evaluate(tree);
        if (res.toString() === expectedStr) {
            console.log(`[PASS] ${name}`);
        } else {
            console.error(`[FAIL] ${name}: Expected "${expectedStr}", got "${res.toString()}"`);
        }
    } catch (e) {
        console.error(`[FAIL] ${name}: Error ${e.message}`);
    }
}

console.log("--- New Features Test ---");

// Test roots(poly, var)
test("roots(x^2 - 1, x)", "roots(x^2 - 1, x)", "[1, -1]");
// Note: _solve for x^2-1 returns set(1, -1) but ordering might vary.
// Let's check logic directly if output varies.
const lexer = new Lexer("roots(x^2 - 4, x)");
const parser = new Parser(lexer);
const tree = parser.parse();
const res = cas.evaluate(tree);
if (res instanceof Vec && res.elements.length === 2) {
    const vals = res.elements.map(e => e.evaluateNumeric()).sort((a,b)=>a-b);
    if (Math.abs(vals[0] - (-2)) < 1e-9 && Math.abs(vals[1] - 2) < 1e-9) {
        console.log("[PASS] roots(x^2-4, x) numeric values");
    } else {
        console.error("[FAIL] roots values mismatch: " + vals);
    }
} else {
    console.error("[FAIL] roots return type: " + res.constructor.name);
}

// Test subs(expr, var=val)
test("subs(x^2, x=2)", "subs(x^2, x=2)", "4");
test("subs(x+y, x=1)", "subs(x+y, x=1)", "(1 + y)");

console.log("--- End New Features Test ---");
