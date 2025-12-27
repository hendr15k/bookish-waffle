
const { CAS } = require('../js/cas');
const { Parser, Lexer } = require('../js/parser');
const expr = require('../js/expression');
Object.assign(globalThis, expr);

globalThis.HELP_DATA = {}; // Mock

const cas = new CAS();

function test(name, input, expected) {
    try {
        const lexer = new Lexer(input);
        const parser = new Parser(lexer);
        const expr = parser.parse();
        const res = cas.evaluate(expr);
        if (res.toString() === expected) {
            console.log(`[PASS] ${name}`);
        } else {
            console.error(`[FAIL] ${name}: Expected ${expected}, got ${res.toString()}`);
            process.exit(1);
        }
    } catch (e) {
        console.error(`[FAIL] ${name}: Error ${e.message}`);
        process.exit(1);
    }
}

test("Union", "union([1, 2], [2, 3])", "[1, 2, 3]");
test("Intersect", "intersect([1, 2], [2, 3])", "[2]");
test("SetDiff", "setdiff([1, 2], [2, 3])", "[1]");
test("Root Simplification", "root(27, 3)", "3");
test("Root Simplification 2", "root(8, 3)", "2");
test("Pow Simplification", "8^(1/3)", "2");
test("Pow Simplification 2", "16^(1/4)", "2");
test("Pow Simplification 3", "27^(2/3)", "9");
console.log("All new feature tests passed");
