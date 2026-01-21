
const assert = require('assert');

// 1. Setup Global Environment (Mimic Browser/Shared Scope)
const expression = require('../js/expression.js');
Object.assign(globalThis, expression);

const parserModule = require('../js/parser.js');
Object.assign(globalThis, parserModule);

const casModule = require('../js/cas.js');
const CAS = casModule.CAS;

const geminiModule = require('../js/gemini_integration.js');
const GeminiAdapter = geminiModule.GeminiAdapter;

console.log("Modules loaded successfully.");

// 2. Initialize
const cas = new CAS();
const adapter = new GeminiAdapter(cas);

// 3. Test getTools() Schema
console.log("Testing getTools()...");
const tools = adapter.getTools();
assert(Array.isArray(tools), "getTools should return an array");
assert(tools.length >= 3, "Should have at least 3 tools");

const calcTool = tools.find(t => t.name === 'calculate');
assert(calcTool, "calculate tool missing");
assert(calcTool.parameters.properties.expression, "calculate tool should have expression param");

const solveTool = tools.find(t => t.name === 'solve_equation');
assert(solveTool, "solve_equation tool missing");

const plotTool = tools.find(t => t.name === 'plot_function');
assert(plotTool, "plot_function tool missing");

console.log("Schema validation passed.");

// 4. Test executeTool()

// Test Calculate
console.log("Testing executeTool('calculate')...");
const calcRes = adapter.executeTool('calculate', { expression: '1 + 1' });
assert.strictEqual(calcRes.status, 'success');
assert.strictEqual(calcRes.result, '2');
console.log("calculate(1+1) ->", calcRes.result);

const calcRes2 = adapter.executeTool('calculate', { expression: 'diff(x^2, x)' });
assert.strictEqual(calcRes2.status, 'success');
// Output might have parentheses
assert(calcRes2.result === '2 * x' || calcRes2.result === '(2 * x)', `Expected 2 * x, got ${calcRes2.result}`);
console.log("calculate(diff(x^2, x)) ->", calcRes2.result);

// Test Solve
console.log("Testing executeTool('solve_equation')...");
const solveRes = adapter.executeTool('solve_equation', { equation: 'x^2 - 4 = 0', variable: 'x' });
assert.strictEqual(solveRes.status, 'success');
// Result might be a Set or list of solutions. String representation depends on CAS implementation.
// Usually "set(-2, 2)" or similar.
console.log("solve(x^2-4=0, x) ->", solveRes.result);
assert(solveRes.result.includes('2'), "Solution should contain 2");
assert(solveRes.result.includes('-2'), "Solution should contain -2");

// Test Plot
console.log("Testing executeTool('plot_function')...");
const plotRes = adapter.executeTool('plot_function', { expression: 'sin(x)', variable: 'x', min: -5, max: 5 });
assert.strictEqual(plotRes.status, 'success');
assert.strictEqual(plotRes.type, 'plot');
assert(plotRes.raw, "Plot result should contain raw object");
assert.strictEqual(plotRes.raw.min, -5);
assert.strictEqual(plotRes.raw.max, 5);
console.log("plot(sin(x)) -> success");

// Test Run Script
console.log("Testing executeTool('run_script')...");
const script = `
a := 5;
b := 10;
a * b
`;
const scriptRes = adapter.executeTool('run_script', { script: script });
if (scriptRes.status === 'error') console.log("Script Error:", scriptRes.error);
assert.strictEqual(scriptRes.status, 'success');
assert.strictEqual(scriptRes.result, '50');
console.log("run_script(a=5; b=10; a*b) ->", scriptRes.result);

// Test Error Handling
console.log("Testing Error Handling...");
const errRes = adapter.executeTool('calculate', { expression: '1 / 0' });
// CAS might handle 1/0 as Infinity or Error?
// let's check parse error
const parseErr = adapter.executeTool('calculate', { expression: '1 + (' }); // Unbalanced parens
assert.strictEqual(parseErr.status, 'error');
console.log("Parse error handled correctly.");

console.log("\nAll Gemini Integration tests passed!");
