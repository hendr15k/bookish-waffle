
const { Lexer, Parser } = require('../js/parser');
const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Assignment, Eq, Vec, FunctionDef, Block } = require('../js/expression');

function parse(text) {
    const lexer = new Lexer(text);
    const parser = new Parser(lexer);
    return parser.parse();
}

function verify(text, expected) {
    try {
        const result = parse(text);
        const simplified = result.simplify();
        const resultStr = simplified.toString();
        if (resultStr === expected) {
            console.log(`[PASS] ${text} -> ${resultStr}`);
        } else {
            console.error(`[FAIL] ${text} -> ${resultStr} (expected: ${expected})`);
            process.exit(1);
        }
    } catch (e) {
        console.error(`[ERROR] ${text} -> ${e.message}`);
        process.exit(1);
    }
}

console.log("--- Testing Division Simplification ---");

// Verification of the fix
verify("x / (2 * x)", "(1 / 2)");
verify("x / (x * 2)", "(1 / 2)"); // Commutative check in denominator
verify("(2 * x) / (3 * x)", "(2 / 3)");
verify("(x * 2) / (x * 3)", "(2 / 3)");
verify("(x * 2) / (3 * x)", "(2 / 3)");
verify("(2 * x) / (x * 3)", "(2 / 3)");

// Existing behavior checks (regression testing for this feature)
verify("(2 * x) / x", "2");
verify("(x * 2) / x", "2");
verify("x / x", "1");
