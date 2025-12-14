
const { CAS } = require('../js/cas.js');
const { Parser, Lexer } = require('../js/parser.js');
require('../js/expression.js'); // Globals

const cas = new CAS();
const parser = new Parser(new Lexer(""));

function test(name, input, expected) {
    try {
        const p = new Parser(new Lexer(input));
        const expr = p.parse();
        const result = cas.evaluate(expr);
        const resStr = result ? result.toString() : "null";

        // Normalize strings for comparison (remove spaces)
        const normalize = s => s.replace(/\s+/g, '');

        if (normalize(resStr) === normalize(expected)) {
            console.log(`[PASS] ${name}`);
        } else {
            console.error(`[FAIL] ${name}`);
            console.error(`  Input:    ${input}`);
            console.error(`  Expected: ${expected}`);
            console.error(`  Actual:   ${resStr}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${name} - Exception`);
        console.error(e.stack);
    }
}

console.log("--- Testing Inequality Solving ---");

// Linear Positive Coeff
test("Linear >", "solve(2*x > 4, x)", "x > 2");
test("Linear <", "solve(2*x < 4, x)", "x < 2");
test("Linear >=", "solve(2*x >= 4, x)", "x >= 2");
test("Linear <=", "solve(2*x <= 4, x)", "x <= 2");

// Linear Negative Coeff (Flip)
test("Linear Neg >", "solve(-2*x > 4, x)", "x < -2");
test("Linear Neg <", "solve(-2*x < 4, x)", "x > -2");
test("Linear Neg >=", "solve(-2*x >= 4, x)", "x <= -2");
test("Linear Neg <=", "solve(-2*x <= 4, x)", "x >= -2");

// Quadratic (Two Roots)
// Result: (x < -2 or x > 2)
test("Quadratic >", "solve(x^2 - 4 > 0, x)", "(x < -2 or x > 2)");

// x^2 - 4 < 0
// Result: (x > -2 and x < 2)
test("Quadratic <", "solve(x^2 < 4, x)", "(x > -2 and x < 2)");

// Quadratic with equality
// x^2 - 4 >= 0
// The solver returns OR of intervals and individual roots:
// (((x < -2 or x > 2) or x = -2) or x = 2)
test("Quadratic >=", "solve(x^2 >= 4, x)", "(((x < -2 or x > 2) or x = -2) or x = 2)");

// No Solution
test("No Solution", "solve(x^2 < -1, x)", "false");

// All Reals
test("All Reals", "solve(x^2 > -1, x)", "true");

// Polynomial
// Result: ((x > 0 and x < 1) or x > 2)
test("Cubic >", "solve(x*(x-1)*(x-2) > 0, x)", "((x > 0 and x < 1) or x > 2)");
