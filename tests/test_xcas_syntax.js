const { Lexer, Parser } = require('../js/parser.js');
const { CAS } = require('../js/cas.js');
require('../js/expression.js');

function test(expr, expected) {
    try {
        const lexer = new Lexer(expr);
        const parser = new Parser(lexer);
        const ast = parser.parse();
        const cas = new CAS();
        const result = cas.evaluate(ast);
        const resStr = result.toString();
        if (resStr === expected) {
            console.log(`PASS: ${expr} -> ${resStr}`);
        } else {
            console.error(`FAIL: ${expr} -> ${resStr} (Expected: ${expected})`);
        }
    } catch (e) {
        console.error(`ERROR: ${expr} -> ${e.message}`);
    }
}

console.log("Testing Xcas Syntax Compatibility...");

// Ranges
test("1..5", "[1, 2, 3, 4, 5]");
test("1.5 .. 4.5", "[1.5, 2.5, 3.5, 4.5]");

// Logic
test("true", "1");
test("false", "0");
test("true and false", "0");
test("true or false", "1");
test("not false", "1");
test("1 and 0", "0");
test("1 or 0", "1");

// Combined
test("(1..3) + 2", "[3, 4, 5]");
