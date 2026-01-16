
const { CAS } = require('../js/cas.js');
const { Lexer, Parser } = require('../js/parser.js');
require('../js/expression.js');

const cas = new CAS();

function eval(str) {
    const lexer = new Lexer(str);
    const parser = new Parser(lexer);
    const ast = parser.parse();
    return cas.evaluate(ast);
}

function check(cmd, expectedPattern, shouldFail = false) {
    try {
        const res = eval(cmd);
        const resStr = res.toString();
        // console.log(`${cmd} -> ${resStr}`);
        const found = resStr.includes(expectedPattern) || new RegExp(expectedPattern).test(resStr);

        if (shouldFail) {
             if (found) console.log(`PASS (Found expected failure mode): ${cmd} -> ${resStr}`);
             else console.error(`FAIL (Unexpected success?): ${cmd} -> ${resStr}`);
        } else {
             if (!found) console.error(`FAIL: ${cmd} -> ${resStr}. Expected ${expectedPattern}`);
             else console.log(`PASS: ${cmd} -> ${resStr}`);
        }
    } catch (e) {
        console.error(`ERROR: ${cmd} -> ${e.message}`);
    }
}

console.log("--- Round 3 ---");

// 1. Solve Cubic (Complex Roots)
// x^3 + x + 1 = 0
// Should return a set of 3 roots (Symbolic mess)
// We just check if it returns a set {...}
check("solve(x^3 + x + 1, x)", "{", false);

// 2. Solve Simple Cubic
// x^3 - 8 = 0 -> 2, ...
check("solve(x^3 - 8, x)", "{", false); // Should have 2, -1+i*sqrt3...

// 3. Factor Multivariate
// x^2 - y^2 -> (x-y)(x+y)
// x^3 - y^3 -> (x-y)(x^2+xy+y^2)
check("factor(x^2 - y^2)", "\\(x - y\\)", false);
check("factor(x^3 - y^3)", "\\(x - y\\)", false);

// 4. Erfi
check("integrate(exp(x^2), x)", "erfi", false);

// 5. Erfi Diff
check("diff(erfi(x), x)", "exp", false);
