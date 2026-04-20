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

try {
  const parsed = parser("odeplot(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))").parse();
  console.log("odeplot node args length", parsed.args.length);
  console.log("odeplot arg 0 funcName", parsed.args[0].funcName);
  const sol = cas.evaluate(parsed.args[0]);
  console.log("sol evaluated", sol instanceof sandbox.Vec);

  const odeplotRes = cas.evaluate(new sandbox.Call('odeplot', [sol]));
  console.log("odeplot:", odeplotRes);
} catch (e) {
  console.error("Error:", e.stack);
}
