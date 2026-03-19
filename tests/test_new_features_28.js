const { Parser, Lexer } = require('../js/parser.js');
const { CAS } = require('../js/cas.js');

const cas = new CAS();

function assertApprox(actual, expected, message, tol = 1e-4) {
    if (Math.abs(actual - expected) > tol) {
        console.error(`FAIL: ${message}. Expected ${expected}, got ${actual}`);
    } else {
        console.log(`PASS: ${message}`);
    }
}

function runTests() {
    console.log("=== Testing Expected Value & Variance ===");

    const parse = (expr) => {
        const lexer = new Lexer(expr);
        const parser = new Parser(lexer);
        return parser.parse();
    };

    // Discrete expected value
    let res = cas.evaluate(parse("expectedValue([1, 2, 3], [0.2, 0.5, 0.3])")).evaluateNumeric();
    assertApprox(res, 2.1, "expectedValue([1, 2, 3], [0.2, 0.5, 0.3])");

    // Discrete variance
    res = cas.evaluate(parse("variance([1, 2, 3], [0.2, 0.5, 0.3])")).evaluateNumeric();
    assertApprox(res, 0.49, "variance([1, 2, 3], [0.2, 0.5, 0.3])");

    // Continuous expected value (uniform 0 to 2, PDF = 0.5)
    res = cas.evaluate(parse("expectedValue(0.5, x, 0, 2)")).evaluateNumeric();
    assertApprox(res, 1, "expectedValue(0.5, x, 0, 2)");

    // Continuous variance (uniform 0 to 2, PDF = 0.5)
    res = cas.evaluate(parse("variance(0.5, x, 0, 2)")).evaluateNumeric();
    assertApprox(res, 1/3, "variance(0.5, x, 0, 2)");

    console.log("=== Testing FFT and IFFT ===");
    res = cas.evaluate(parse("fft([1, 0, 1, 0])"));
    let elems = res.elements.map(e => e.evaluateNumeric());
    assertApprox(elems[0], 2, "fft([1, 0, 1, 0]) element 0");
    assertApprox(elems[1], 0, "fft([1, 0, 1, 0]) element 1");
    assertApprox(elems[2], 2, "fft([1, 0, 1, 0]) element 2");
    assertApprox(elems[3], 0, "fft([1, 0, 1, 0]) element 3");

    res = cas.evaluate(parse("ifft([2, 0, 2, 0])"));
    elems = res.elements.map(e => e.evaluateNumeric());
    assertApprox(elems[0], 1, "ifft([2, 0, 2, 0]) element 0");
    assertApprox(elems[1], 0, "ifft([2, 0, 2, 0]) element 1");
    assertApprox(elems[2], 1, "ifft([2, 0, 2, 0]) element 2");
    assertApprox(elems[3], 0, "ifft([2, 0, 2, 0]) element 3");

    console.log("=== Testing RK4 Solver ===");
    res = cas.evaluate(parse("rk4(-2*y, y, t, 1, 0, 0.5, 0.1)"));
    const lastPoint = res.elements[res.elements.length - 1];
    assertApprox(lastPoint.elements[0].evaluateNumeric(), 0.5, "rk4 final t");
    assertApprox(lastPoint.elements[1].evaluateNumeric(), Math.exp(-1), "rk4 final y", 1e-4);

    console.log("=== Testing Probability Distributions ===");
    res = cas.evaluate(parse("normalPDF(0, 0, 1)")).evaluateNumeric();
    assertApprox(res, 1/Math.sqrt(2*Math.PI), "normalPDF(0, 0, 1)");

    res = cas.evaluate(parse("normalCDF(0, 0, 1)")).evaluateNumeric();
    assertApprox(res, 0.5, "normalCDF(0, 0, 1)");

    res = cas.evaluate(parse("exponentialPDF(1, 2)")).evaluateNumeric();
    assertApprox(res, 2 * Math.exp(-2), "exponentialPDF(1, 2)");

    res = cas.evaluate(parse("exponentialCDF(1, 2)")).evaluateNumeric();
    assertApprox(res, 1 - Math.exp(-2), "exponentialCDF(1, 2)");
}

runTests();
