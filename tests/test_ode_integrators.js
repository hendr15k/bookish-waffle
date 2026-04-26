const fs = require("fs");
const vm = require("vm");
function loadFile(filePath) { return fs.readFileSync(filePath, "utf8"); }
const expressionCode = loadFile("js/expression.js");
const parserCode = loadFile("js/parser.js");
const casCode = loadFile("js/cas.js");
const helpCode = loadFile("js/help.js");
const sandbox = { console, Math, Number, parseFloat, parseInt, isNaN, window: {} };
vm.createContext(sandbox);
vm.runInContext(expressionCode + "\n" + parserCode + "\n" + helpCode + "\n" + casCode + "\n" + `globalThis.Expr=Expr;globalThis.Num=Num;globalThis.Sym=Sym;globalThis.Add=Add;globalThis.Sub=Sub;globalThis.Mul=Mul;globalThis.Div=Div;globalThis.Pow=Pow;globalThis.Call=Call;globalThis.Assignment=Assignment;globalThis.Eq=Eq;globalThis.Vec=Vec;globalThis.FunctionDef=FunctionDef;globalThis.Lexer=Lexer;globalThis.Parser=Parser;globalThis.CAS=CAS;`, sandbox);
const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));
const evalExpr = (text) => cas.evaluate(parser(text).parse());
function ok(desc, input) {
  try { const r = evalExpr(input); console.log(`[PASS] ${desc}: ${r.toString()}`); }
  catch (e) { console.error(`[FAIL] ${desc}: ${e.message}`); process.exitCode = 1; }
}

const assert = (desc, condition) => { if (!condition) throw new Error(desc); };

const r1 = evalExpr("rk4(y-t, y, t, 1, 0, 1, 0.1)");
assert("rk4 returns a vector", r1 instanceof sandbox.Vec);
assert("rk4 has multiple rows", r1.elements.length > 1);

const r2 = evalExpr("rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, [])");
assert("rk45 returns a vector", r2 instanceof sandbox.Vec);
assert("rk45 includes metadata", r2.elements[r2.elements.length - 1] instanceof sandbox.Vec && r2.elements[r2.elements.length - 1].elements[0] instanceof sandbox.Sym);

ok("odeplot creates plot object", "odeplot(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
ok("ode_plot creates plot object", "ode_plot(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
ok("odetable creates table", "odetable(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
ok("ode_table creates table", "ode_table(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
ok("harmonic oscillator second-order setup", "rk4([y2, -1*y1], [y1, y2], t, [1, 0], 0, 1, 0.1)");
ok("damped oscillation example", "rk45([y2, -0.2*y2 - y1], [y1, y2], t, [1, 0], 0, 10, 0.1, 1e-6, [])");
ok("exponential decay example", "rk45(-0.7*y, y, t, 1, 0, 5, 0.1, 1e-6, [])");
ok("SIR model example", "rk4([-beta*s*i, beta*s*i-gamma*i, gamma*i], [s, i, r], t, [0.99, 0.01, 0], 0, 10, 0.1)");

// ── Bogacki-Shampine 4(5) tests ──────────────────────────────────────────────
const r3 = evalExpr("bs45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, [])");
assert("bs45 returns a vector", r3 instanceof sandbox.Vec);
assert("bs45 has multiple rows", r3.elements.length > 1);
assert("bs45 includes metadata row", r3.elements[r3.elements.length - 1] instanceof sandbox.Vec && r3.elements[r3.elements.length - 1].elements[0] instanceof sandbox.Sym);
assert("bs45 first column is time", r3.elements[0] instanceof sandbox.Vec && r3.elements[0].elements[0] instanceof sandbox.Num);
assert("bs45 metadata __meta__ label", r3.elements[r3.elements.length - 1].elements[0].name === '__meta__');

ok("bs45 exponential decay", "bs45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, [])");
ok("bs45 harmonic oscillator", "bs45([y2, -y1], [y1, y2], t, [0, 1], 0, 10, 0.1, 1e-6, [])");
ok("bs45 damped oscillation", "bs45([y2, -0.2*y2 - y1], [y1, y2], t, [1, 0], 0, 10, 0.1, 1e-6, [])");
ok("bs45 SIR model", "bs45([-beta*s*i, beta*s*i-gamma*i, gamma*i], [s, i, r], t, [0.99, 0.01, 0], 0, 10, 0.1, 1e-4, [])");
ok("odetable works with bs45", "odetable(bs45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
ok("odeplot works with bs45", "odeplot(bs45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))");
