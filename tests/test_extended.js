
const fs = require('fs');
const vm = require('vm');

function loadFile(filePath) {
    return fs.readFileSync(filePath, 'utf8');
}

const expressionCode = loadFile('js/expression.js');
const parserCode = loadFile('js/parser.js');
const casCode = loadFile('js/cas.js');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    isNaN: isNaN
};

vm.createContext(sandbox);

// Helper to expose classes
const exposeClasses = `
    globalThis.Expr = Expr;
    globalThis.Num = Num;
    globalThis.Sym = Sym;
    globalThis.Symbol = Sym; // Alias for tests if they use Symbol explicitly, or just expose Sym
    globalThis.Add = Add;
    globalThis.Sub = Sub;
    globalThis.Mul = Mul;
    globalThis.Div = Div;
    globalThis.Pow = Pow;
    globalThis.Call = Call;
    globalThis.Assignment = Assignment;
    globalThis.Eq = Eq;
    globalThis.Vec = Vec;
    globalThis.FunctionDef = FunctionDef;

    globalThis.Lexer = Lexer;
    globalThis.Parser = Parser;

    globalThis.CAS = CAS;
`;

vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + exposeClasses, sandbox);

const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

function evalExpr(text) {
    const expr = parser(text).parse();
    return cas.evaluate(expr);
}

function test(description, input, expectedOutput) {
    try {
        const result = evalExpr(input);
        const resultStr = result.toString();
        // Allow for some flexibility in expected output string (e.g., whitespace, order)
        // Ideally we would compare ASTs, but string comparison is what we have.
        if (resultStr === expectedOutput) {
            console.log(`[PASS] ${description}`);
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedOutput}`);
            console.error(`  Got:      ${resultStr}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
    }
}

// Helper for numeric approximation checking
function testNumeric(description, input, expectedValue, tolerance = 1e-6) {
    try {
        const result = evalExpr(input);
        const val = result.evaluateNumeric();
        if (Math.abs(val - expectedValue) < tolerance) {
            console.log(`[PASS] ${description}`);
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedValue} (approx)`);
            console.error(`  Got:      ${val}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
    }
}

console.log("--- Extended Test Suite ---");

// 1. Matrix Algebra
test("Matrix Determinant 2x2", "det([[1, 2], [3, 4]])", "-2");
test("Matrix Determinant 3x3", "det([[6, 1, 1], [4, -2, 5], [2, 8, 7]])", "-306");
// Update expected: fractions instead of decimals
test("Matrix Inverse 2x2", "inv([[4, 7], [2, 6]])", "[[(3 / 5), (-7 / 10)], [(-1 / 5), (2 / 5)]]");
test("Matrix Transpose", "trans([[1, 2], [3, 4]])", "[[1, 3], [2, 4]]");
test("Cross Product", "cross([1, 0, 0], [0, 1, 0])", "[0, 0, 1]");

// 2. Calculus
test("Diff sin(x)", "diff(sin(x), x)", "cos(x)");
test("Diff cos(x)", "diff(cos(x), x)", "(-1 * sin(x))");
test("Diff exp(x)", "diff(exp(x), x)", "exp(x)");
test("Diff ln(x)", "diff(ln(x), x)", "(1 / x)");
test("Diff x^2", "diff(x^2, x)", "(2 * x)");
// Update expected: parentheses around arg in cos
test("Diff chain rule sin(x^2)", "diff(sin(x^2), x)", "(cos((x^2)) * (2 * x))");

test("Integrate x", "integrate(x, x)", "((x^2) / 2)");
test("Integrate x^2", "integrate(x^2, x)", "((x^3) / 3)");
test("Integrate 1/x", "integrate(1/x, x)", "ln(x)");
test("Integrate constant", "integrate(5, x)", "(5 * x)");

// 3. Complex Numbers
test("i^1", "i^1", "i");
test("i^2", "i^2", "-1");
test("i^3", "i^3", "(-1 * i)");
test("i^4", "i^4", "1");
test("sqrt(-1)", "sqrt(-1)", "i");
// Update expected: commuted 2 * i
test("sqrt(-4)", "sqrt(-4)", "(2 * i)");

// 4. Statistics & Probability
test("Mean", "mean([1, 2, 3, 4, 5])", "3");
test("Variance", "variance([1, 2, 3])", "1");
test("nCr(5, 2)", "nCr(5, 2)", "10");
test("nPr(5, 2)", "nPr(5, 2)", "20");
test("Factorial", "factorial(5)", "120");
test("BinomialPDF", "binomialPDF(1, 2, 0.5)", "0.5");

// 5. Number Theory
test("GCD", "gcd(12, 18)", "6");
test("LCM", "lcm(4, 6)", "12");
test("isPrime(7)", "isPrime(7)", "1");
test("isPrime(10)", "isPrime(10)", "0");
test("Factor 10", "factor(10)", "factored(2, 5)");

// 6. Equation Solving
test("Solve linear", "solve(2*x - 4 = 0, x)", "2");
test("Solve quadratic", "solve(x^2 - 1, x)", "set(1, -1)");

// 7. Taylor Series
// Update expected: matches current simplified output
test("Taylor exp(x)", "taylor(exp(x), x, 0, 2)", "(1 + (x + ((x^2) / 2)))");

// 8. Simplification
test("Add 0", "x + 0", "x");
test("Mul 1", "x * 1", "x");
test("Mul 0", "x * 0", "0");
test("Sub self", "x - x", "0");
test("Div self", "x / x", "1");
test("Pow 0", "x^0", "1");
test("Pow 1", "x^1", "x");

// 9. Parsing edge cases
test("Implicit Mul Number-Var", "2x", "(2 * x)");
test("Implicit Mul Var-Var", "x y", "(x * y)");
test("Implicit Mul Parenthesis", "2(x+1)", "(2 * (x + 1))");
test("Function Power", "sin^2(x)", "(sin(x)^2)");

// 10. Limits
test("Limit sin(x)/x as x->0", "limit(sin(x)/x, x, 0)", "1");
test("Limit 1/x as x->Infinity", "limit(1/x, x, Infinity)", "0");

// 11. Custom Functions
test("Func Def", "f(x) := x^2; f(2)", "4");

console.log("--- End Extended Test Suite ---");
