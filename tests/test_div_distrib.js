
const fs = require('fs');
const vm = require('vm');

const code1 = fs.readFileSync('js/expression.js', 'utf8');

const sandbox = { console, Math };
vm.createContext(sandbox);
vm.runInContext(code1, sandbox);

vm.runInContext(`
    const x = new Sym('x');
    const y = new Sym('y');
    const z = new Sym('z');
    const two = new Num(2);
    const three = new Num(3);
    const four = new Num(4);

    let failures = 0;
    const assertEq = (desc, expr, expected) => {
        const actual = expr.toString();
        if (actual === expected) {
            console.log("[PASS] " + desc + ": " + actual);
        } else {
            console.log("[FAIL] " + desc + ": expected " + expected + ", got " + actual);
            failures++;
        }
    };

    // Case 1: (2x + 2y) / 2 -> x + y
    const sum1 = new Add(new Mul(two, x), new Mul(two, y));
    const div1 = new Div(sum1, two).simplify();
    assertEq("(2x + 2y) / 2", div1, "(x + y)");

    // Case 2: (4x + 2) / 2 -> 2x + 1
    const sum2 = new Add(new Mul(four, x), two);
    const div2 = new Div(sum2, two).simplify();
    assertEq("(4x + 2) / 2", div2, "((2 * x) + 1)");

    // Case 3: (Ax + Bx) / x -> A + B
    const A = new Sym('A');
    const B = new Sym('B');
    const sum3 = new Add(new Mul(A, x), new Mul(B, x));
    const div3 = new Div(sum3, x).simplify();
    assertEq("(Ax + Bx) / x", div3, "(A + B)");

    // Case 4: (2x + 3y) / 2 -> (2x + 3y) / 2 (Should NOT distribute into partial fractions)
    const sum4 = new Add(new Mul(two, x), new Mul(three, y));
    const div4 = new Div(sum4, two).simplify();
    // expected: stays as division
    assertEq("(2x + 3y) / 2", div4, "(((2 * x) + (3 * y)) / 2)");

    // Case 5: (x + y) / z -> (x + y) / z (Should NOT distribute)
    const sum5 = new Add(x, y);
    const div5 = new Div(sum5, z).simplify();
    assertEq("(x + y) / z", div5, "((x + y) / z)");

    // Case 6: Subtraction (2x - 2y) / 2 -> x - y
    const sub6 = new Sub(new Mul(two, x), new Mul(two, y));
    const div6 = new Div(sub6, two).simplify();
    assertEq("(2x - 2y) / 2", div6, "(x - y)");

    if (failures > 0) {
        console.log("Total Failures: " + failures);
        // process.exit(1); // Do not exit, just log
    } else {
        console.log("All tests passed!");
    }

`, sandbox);
