const fs = require('fs');
const vm = require('vm');

function loadFile(filePath) {
    return fs.readFileSync(filePath, 'utf8');
}

const expressionCode = loadFile('js/expression.js');
const parserCode = loadFile('js/parser.js');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    isNaN: isNaN
};

vm.createContext(sandbox);
vm.runInContext(expressionCode + "\n" + parserCode, sandbox);

const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

function test(description, input, expectedOutput) {
    try {
        const result = parser(input).parse().simplify();
        const resultStr = result.toString();
        if (resultStr === expectedOutput) {
            console.log(`[PASS] ${description}`);
            return true;
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedOutput}`);
            console.error(`  Got:      ${resultStr}`);
            return false;
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
        return false;
    }
}

console.log('=== Testing BUG #8: nested Add partial cancellation ===');

let passed = 0;
let failed = 0;

if (test('(x+y+z)-(x+y) -> z', '(x+y+z)-(x+y)', 'z')) passed++; else failed++;
if (test('(z+x+y)-(y+x) -> z', '(z+x+y)-(y+x)', 'z')) passed++; else failed++;
if (test('(a+b)-(b+a) -> 0', '(a+b)-(b+a)', '0')) passed++; else failed++;
if (test('(x+1)-(1+x) -> 0', '(x+1)-(1+x)', '0')) passed++; else failed++;
if (test('(w+x+y+z)-(w+x) -> y+z', '(w+x+y+z)-(w+x)', '(y + z)')) passed++; else failed++;
if (test('(a+b+c+d)-(c+a) -> b+d', '(a+b+c+d)-(c+a)', '(b + d)')) passed++; else failed++;
if (test('(x+y)-(x) -> y', '(x+y)-(x)', 'y')) passed++; else failed++;
if (test('(x+y)-(y) -> x', '(x+y)-(y)', 'x')) passed++; else failed++;

console.log(`\nResults: ${passed} passed, ${failed} failed`);
process.exit(failed > 0 ? 1 : 0);
