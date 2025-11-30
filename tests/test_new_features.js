
const assert = require('assert');
const vm = require('vm');
const fs = require('fs');
const path = require('path');

// Load the CAS classes
const loadFile = (filePath) => fs.readFileSync(path.join(__dirname, '..', filePath), 'utf8');
const context = {};
vm.createContext(context);
vm.runInContext(loadFile('js/expression.js'), context);
vm.runInContext(loadFile('js/cas.js'), context);
vm.runInContext(loadFile('js/parser.js'), context);

const { CAS, Lexer, Parser, Num, Vec, Call, Mul, Add, Sub } = context;

const cas = new CAS();

function evalExpr(input) {
    const lexer = new Lexer(input);
    const parser = new Parser(lexer);
    const expr = parser.parse();
    return cas.evaluate(expr);
}

function assertApprox(val, expected, tolerance = 1e-9) {
    assert(Math.abs(val - expected) < tolerance, `Expected ${expected}, got ${val}`);
}

// Test new commands

try {

// 1. dot(u, v)
console.log("Testing dot(u, v)...");
let res = evalExpr("dot([1, 2, 3], [4, 5, 6])");
assert.strictEqual(res.toString(), "32"); // 1*4 + 2*5 + 3*6 = 4+10+18 = 32

// 2. norm(v)
console.log("Testing norm(v)...");
res = evalExpr("norm([3, 4])");
assert.strictEqual(res.toString(), "5"); // sqrt(3^2 + 4^2) = 5

// 3. grad(expr, vars)
console.log("Testing grad(f, V)...");
// grad(x^2 + y^2, [x, y]) -> [2x, 2y]
res = evalExpr("grad(x^2 + y^2, [x, y])");
assert.strictEqual(res.toString(), "[(2 * x), (2 * y)]"); // Mul simplification might make it (2*x) or (2 * x)

// 4. curl(F, V)
console.log("Testing curl(F, V)...");
// curl([-y, x, 0], [x, y, z]) -> [0, 0, 1+1=2]
res = evalExpr("curl([-y, x, 0], [x, y, z])");
// c3 = dy(x) - dx(-y) = 1 - (-1) = 2.
assert.strictEqual(res.toString(), "[0, 0, 2]");

// 5. divergence(F, V)
console.log("Testing divergence(F, V)...");
// div([x, y, z], [x, y, z]) -> 1+1+1 = 3
res = evalExpr("divergence([x, y, z], [x, y, z])");
assert.strictEqual(res.toString(), "3");

// 6. rem, quo, mod
console.log("Testing rem, quo, mod...");
res = evalExpr("rem(10, 3)");
assert.strictEqual(res.toString(), "1");
res = evalExpr("quo(10, 3)");
assert.strictEqual(res.toString(), "3");
res = evalExpr("mod(-1, 4)");
assert.strictEqual(res.toString(), "3");

// 7. size(L)
console.log("Testing size(L)...");
res = evalExpr("size([1, 2, 3])");
assert.strictEqual(res.toString(), "3");

// 8. concat(L1, L2)
console.log("Testing concat(L1, L2)...");
res = evalExpr("concat([1, 2], [3, 4])");
assert.strictEqual(res.toString(), "[1, 2, 3, 4]");

// 9. approx(expr)
console.log("Testing approx(expr)...");
res = evalExpr("approx(1/2)");
assert.strictEqual(res.toString(), "0.5");

// 10. arg(z)
console.log("Testing arg(z)...");
res = evalExpr("arg(1)");
assert.strictEqual(res.toString(), "0");
// arg(-1) should be pi.
res = evalExpr("arg(-1)");
// It might be Num(Math.PI) which prints as 3.14...
// Let's check value directly.
assertApprox(res.value, Math.PI);

console.log("All tests passed!");

} catch (e) {
    console.error(e);
    process.exit(1);
}
