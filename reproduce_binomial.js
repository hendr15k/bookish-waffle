
const fs = require('fs');
const vm = require('vm');

// Load libraries
const code1 = fs.readFileSync('js/expression.js', 'utf8');
const code2 = fs.readFileSync('js/parser.js', 'utf8');
const code3 = fs.readFileSync('js/cas.js', 'utf8');

const sandbox = { console, Math };
vm.createContext(sandbox);
vm.runInContext(code1, sandbox);
vm.runInContext(code2, sandbox);
vm.runInContext(code3, sandbox);

vm.runInContext(`
    const cas = new CAS();

    // Test Binomial PDF
    // n=10, p=0.5, k=5 -> 10C5 * 0.5^5 * 0.5^5 = 252 * (0.5)^10 = 252 / 1024 approx 0.24609375
    const res = cas._binomialPDF(new Num(5), new Num(10), new Num(0.5));
    console.log("PDF(5, 10, 0.5) =", res.toString());

    if (Math.abs(res.value - 0.24609375) < 1e-9) {
        console.log("PASS");
    } else {
        console.log("FAIL");
    }

    // Test nCr
    const c = cas._nCr(new Num(10), new Num(5));
    console.log("10C5 =", c.toString()); // Should be 252
`, sandbox);
