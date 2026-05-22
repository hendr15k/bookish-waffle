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
vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + `
    globalThis.Expr = Expr;
    globalThis.Num = Num;
    globalThis.Sym = Sym;
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
`, sandbox);

const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

function evalExpr(text) {
    const expr = parser(text).parse();
    return cas.evaluate(expr);
}

let passed = 0;
let failed = 0;

function test(description, input, expectedOutput) {
    try {
        const result = evalExpr(input);
        const resultStr = result.toString();
        if (resultStr === expectedOutput) {
            console.log(`[PASS] ${description}`);
            passed++;
            return true;
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedOutput}`);
            console.error(`  Got:      ${resultStr}`);
            failed++;
            return false;
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
        failed++;
        return false;
    }
}

console.log('=== Testing _quo polynomial long division ===');
test('quo(10, 3)', 'quo(10, 3)', '3');
test('quo(x^2+1, x-1)', 'quo(x^2+1, x-1)', '(x + 1)');
test('quo(x^3+1, x^2-1)', 'quo(x^3+1, x^2-1)', 'x');
test('quo(x^3, x-1)', 'quo(x^3, x-1)', '(x^2 + (x + 1))');
test('quo(x^4-1, x^2-1)', 'quo(x^4-1, x^2-1)', '(x^2 + 1)');
test('quo(2*x^2+3*x+1, 2*x+1)', 'quo(2*x^2+3*x+1, 2*x+1)', '(x + 1)');

console.log('\n=== Testing _rem polynomial long division ===');
test('rem(10, 3)', 'rem(10, 3)', '1');
test('rem(x^2+1, x-1)', 'rem(x^2+1, x-1)', '2');
test('rem(x^3+1, x^2-1)', 'rem(x^3+1, x^2-1)', '(x + 1)');
test('rem(x^4-1, x^2-1)', 'rem(x^4-1, x^2-1)', '0');
test('rem(x^3, x-1)', 'rem(x^3, x-1)', '1');
test('rem(x^2, x-1)', 'rem(x^2, x-1)', '1');
test('rem(6*x^2+5*x+1, 2*x+1)', 'rem(6*x^2+5*x+1, 2*x+1)', '0');

console.log('\n=== Partfrac integration (uses _quo) ===');
test('partfrac((x^2+1)/(x-1), x)', 'partfrac((x^2+1)/(x-1), x)', '(x + (1 + (2 / (x - 1))))');
test('partfrac((x^3+1)/(x^2-1), x)', 'partfrac((x^3+1)/(x^2-1), x)', '(x + (1 / (x - 1)))');
test('partfrac(1/(x^2-1), x)', 'partfrac(1/(x^2-1), x)', '((-1 / (2 * (x + 1))) + (1 / (2 * (x - 1))))');

console.log('\n=== Integration using _quo ===');
test('integrate(x^2/(x^2+1), x)', 'integrate(x^2/(x^2+1), x)', '(x + (0 - atan(x)))');

console.log(`\n${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
