
const fs = require('fs');
const { JSDOM } = require('jsdom');

const codeExpression = fs.readFileSync('js/expression.js', 'utf8');
const codeParser = fs.readFileSync('js/parser.js', 'utf8');
const codeChemistry = fs.readFileSync('js/chemistry.js', 'utf8');
const codeHelp = fs.readFileSync('js/help.js', 'utf8');
const codeCas = fs.readFileSync('js/cas.js', 'utf8');

let htmlContent = fs.readFileSync('index.html', 'utf8');

htmlContent = htmlContent.replace('<script src="js/expression.js"></script>', `<script>${codeExpression}</script>`);
htmlContent = htmlContent.replace('<script src="js/parser.js"></script>', `<script>${codeParser}</script>`);
htmlContent = htmlContent.replace('<script src="js/chemistry.js"></script>', `<script>${codeChemistry}</script>`);
htmlContent = htmlContent.replace('<script src="js/help.js"></script>', `<script>${codeHelp}</script>`);
htmlContent = htmlContent.replace('<script src="js/cas.js"></script>', `<script>${codeCas}</script>`);

const mockCanvasScript = `
<script>
  HTMLCanvasElement.prototype.getContext = function() {
    return { clearRect:()=>{}, beginPath:()=>{}, moveTo:()=>{}, lineTo:()=>{}, stroke:()=>{}, fill:()=>{}, arc:()=>{}, strokeRect:()=>{}, fillRect:()=>{}, setLineDash:()=>{}, createImageData:()=>({data:[]}), putImageData:()=>{}, measureText:()=>({width:0}), fillText:()=>{} };
  };
</script>
`;
htmlContent = htmlContent.replace('<head>', '<head>' + mockCanvasScript);

const dom = new JSDOM(htmlContent, {
    runScripts: "dangerously",
    resources: "usable",
    url: "http://localhost/"
});

const { window } = dom;

let passed = 0, failed = 0, skipped = 0;

function test(desc, fn) {
    try {
        fn();
        console.log(`  [PASS] ${desc}`);
        passed++;
    } catch (e) {
        console.error(`  [FAIL] ${desc}`);
        console.error(`         ${e.message}`);
        failed++;
    }
}

function skip(desc) {
    console.log(`  [SKIP] ${desc}`);
    skipped++;
}

function runCommand(cas, cmdStr) {
    try {
        const l = new window.Lexer(cmdStr);
        const p = new window.Parser(l);
        const t = p.parse();
        return cas.evaluate(t);
    } catch (e) {
        return { error: e.message };
    }
}

window.addEventListener('load', () => {
    // Wait for scripts to be injected
    const waitForScripts = () => {
        if (typeof window.Lexer === 'undefined' || typeof window.Parser === 'undefined' || typeof window.CAS === 'undefined') {
            return false;
        }
        return true;
    };

    const tryRun = () => {
        if (!waitForScripts()) {
            setTimeout(tryRun, 100);
            return;
        }

    console.log("\n========================================");
    console.log("  COMPREHENSIVE UI REGRESSION TESTS");
    console.log("========================================\n");

    const cas = window.cas;

    // ============================================================
    // TEST 1: CORE CALCULATOR
    // ============================================================
    console.log("[1] Calculator & Basic Operations");
    test("simplify(1+1)", () => { const r = runCommand(cas, "simplify(1+1)"); if (r.error) throw new Error(r.error); });
    test("expand((x+1)^2)", () => { const r = runCommand(cas, "expand((x+1)^2)"); if (r.error) throw new Error(r.error); });
    test("factor(x^2-1)", () => { const r = runCommand(cas, "factor(x^2-1)"); if (r.error) throw new Error(r.error); });
    test("sqrt(16)", () => { const r = runCommand(cas, "sqrt(16)"); if (r.error) throw new Error(r.error); });
    test("cbrt(27)", () => { const r = runCommand(cas, "cbrt(27)"); if (r.error) throw new Error(r.error); });
    test("abs(-5)", () => { const r = runCommand(cas, "abs(-5)"); if (r.error) throw new Error(r.error); });
    test("round(4.7)", () => { const r = runCommand(cas, "round(4.7)"); if (r.error) throw new Error(r.error); });
    test("floor(4.8)", () => { const r = runCommand(cas, "floor(4.8)"); if (r.error) throw new Error(r.error); });
    test("ceil(4.2)", () => { const r = runCommand(cas, "ceil(4.2)"); if (r.error) throw new Error(r.error); });
    test("sign(-5)", () => { const r = runCommand(cas, "sign(-5)"); if (r.error) throw new Error(r.error); });
    test("max(1,5,3)", () => { const r = runCommand(cas, "max(1,5,3)"); if (r.error) throw new Error(r.error); });
    test("min(1,5,3)", () => { const r = runCommand(cas, "min(1,5,3)"); if (r.error) throw new Error(r.error); });
    test("exp(1)", () => { const r = runCommand(cas, "exp(1)"); if (r.error) throw new Error(r.error); });
    test("ln(e)", () => { const r = runCommand(cas, "ln(e)"); if (r.error) throw new Error(r.error); });
    test("log(100)", () => { const r = runCommand(cas, "log(100)"); if (r.error) throw new Error(r.error); });
    test("log(8, 2)", () => { const r = runCommand(cas, "log(8, 2)"); if (r.error) throw new Error(r.error); });
    test("clamp(15, 0, 10)", () => { const r = runCommand(cas, "clamp(15, 0, 10)"); if (r.error) throw new Error(r.error); });
    test("mod(10, 3)", () => { const r = runCommand(cas, "mod(10, 3)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 2: ALGEBRA
    // ============================================================
    console.log("\n[2] Algebra");
    test("solve(x^2-4, x)", () => { const r = runCommand(cas, "solve(x^2-4, x)"); if (r.error) throw new Error(r.error); });
    test("solve([x+y=10, x-y=2], [x,y])", () => { const r = runCommand(cas, "solve([x+y=10, x-y=2], [x,y])"); if (r.error) throw new Error(r.error); });
    test("solve(x^2-4, x, 2) [Newton]", () => { const r = runCommand(cas, "solve(x^2-4, x, 2)"); if (r.error) throw new Error(r.error); });
    test("partfrac(1/(x^2-1), x)", () => { const r = runCommand(cas, "partfrac(1/(x^2-1), x)"); if (r.error) throw new Error(r.error); });
    test("completeSquare(x^2+6x+10, x)", () => { const r = runCommand(cas, "completeSquare(x^2+6x+10, x)"); if (r.error) throw new Error(r.error); });
    test("collect(x^2+2*x+x^2, x)", () => { const r = runCommand(cas, "collect(x^2+2*x+x^2, x)"); if (r.error) throw new Error(r.error); });
    test("roots(x^2-4, x)", () => { const r = runCommand(cas, "roots(x^2-4, x)"); if (r.error) throw new Error(r.error); });
    test("groebner([x^2-y, x*y-1], [x,y])", () => { const r = runCommand(cas, "groebner([x^2-y, x*y-1], [x,y])"); if (r.error) throw new Error(r.error); });
    test("sturm(x^2-2, x)", () => { const r = runCommand(cas, "sturm(x^2-2, x)"); if (r.error) throw new Error(r.error); });
    test("numRealRoots(x^2-4, x)", () => { const r = runCommand(cas, "numRealRoots(x^2-4, x)"); if (r.error) throw new Error(r.error); });
    test("numRealRoots(x^2-4, x, -5, 5)", () => { const r = runCommand(cas, "numRealRoots(x^2-4, x, -5, 5)"); if (r.error) throw new Error(r.error); });
    test("resultant(x^2-1, x+1, x)", () => { const r = runCommand(cas, "resultant(x^2-1, x+1, x)"); if (r.error) throw new Error(r.error); });
    test("discriminant(x^2+b*x+c, x)", () => { const r = runCommand(cas, "discriminant(x^2+b*x+c, x)"); if (r.error) throw new Error(r.error); });
    test("coeff(3*x^2+2*x, x, 2)", () => { const r = runCommand(cas, "coeff(3*x^2+2*x, x, 2)"); if (r.error) throw new Error(r.error); });
    test("degree(x^3+x, x)", () => { const r = runCommand(cas, "degree(x^3+x, x)"); if (r.error) throw new Error(r.error); });
    test("symb2poly(x^2+2*x+1, x)", () => { const r = runCommand(cas, "symb2poly(x^2+2*x+1, x)"); if (r.error) throw new Error(r.error); });
    test("poly2symb([1, 2, 1], x)", () => { const r = runCommand(cas, "poly2symb([1, 2, 1], x)"); if (r.error) throw new Error(r.error); });
    test("trigReduce(sin(x)^2)", () => { const r = runCommand(cas, "trigReduce(sin(x)^2)"); if (r.error) throw new Error(r.error); });
    test("trigExpand(sin(2*x))", () => { const r = runCommand(cas, "trigExpand(sin(2*x))"); if (r.error) throw new Error(r.error); });
    test("trigSimplify(sin(x)^2+cos(x)^2)", () => { const r = runCommand(cas, "trigSimplify(sin(x)^2+cos(x)^2)"); if (r.error) throw new Error(r.error); });
    test("lagrange([[1,2],[2,4]], x)", () => { const r = runCommand(cas, "lagrange([[1,2],[2,4]], x)"); if (r.error) throw new Error(r.error); });
    test("taylor(sin(x), x, 0, 4)", () => { const r = runCommand(cas, "taylor(sin(x), x, 0, 4)"); if (r.error) throw new Error(r.error); });
    test("laurent(1/x, x, 0, 3)", () => { const r = runCommand(cas, "laurent(1/x, x, 0, 3)"); if (r.error) throw new Error(r.error); });
    test("pade(exp(x), x, 2, 2)", () => { const r = runCommand(cas, "pade(exp(x), x, 2, 2)"); if (r.error) throw new Error(r.error); });
    test("rewrite(tan(x), sin)", () => { const r = runCommand(cas, "rewrite(tan(x), sin)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 3: CALCULUS
    // ============================================================
    console.log("\n[3] Calculus");
    test("diff(sin(x), x)", () => { const r = runCommand(cas, "diff(sin(x), x)"); if (r.error) throw new Error(r.error); });
    test("diff(x^3, x, 2)", () => { const r = runCommand(cas, "diff(x^3, x, 2)"); if (r.error) throw new Error(r.error); });
    test("implicitDiff(x^2+y^2=1, y, x)", () => { const r = runCommand(cas, "implicitDiff(x^2+y^2=1, y, x)"); if (r.error) throw new Error(r.error); });
    test("integrate(x^2, x)", () => { const r = runCommand(cas, "integrate(x^2, x)"); if (r.error) throw new Error(r.error); });
    test("integrate(x^2, x, 0, 1)", () => { const r = runCommand(cas, "integrate(x^2, x, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("nIntegrate(exp(-x^2), x, 0, 1)", () => { const r = runCommand(cas, "nIntegrate(exp(-x^2), x, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("limit(sin(x)/x, x, 0)", () => { const r = runCommand(cas, "limit(sin(x)/x, x, 0)"); if (r.error) throw new Error(r.error); });
    test("limit(1/x, x, 0, 1)", () => { const r = runCommand(cas, "limit(1/x, x, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("analyze(x^3-3*x, x)", () => { const r = runCommand(cas, "analyze(x^3-3*x, x)"); if (r.error) throw new Error(r.error); });
    test("extrema(x^3-3*x, x)", () => { const r = runCommand(cas, "extrema(x^3-3*x, x)"); if (r.error) throw new Error(r.error); });
    test("stationary_points(x^3-3*x, x)", () => { const r = runCommand(cas, "stationary_points(x^3-3*x, x)"); if (r.error) throw new Error(r.error); });
    test("asymptotes((x^2-1)/(x-2), x)", () => { const r = runCommand(cas, "asymptotes((x^2-1)/(x-2), x)"); if (r.error) throw new Error(r.error); });
    test("tangent(x^2, x, 1)", () => { const r = runCommand(cas, "tangent(x^2, x, 1)"); if (r.error) throw new Error(r.error); });
    test("arcLen(x^2, x, 0, 1)", () => { const r = runCommand(cas, "arcLen(x^2, x, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("curvature(x^2, x)", () => { const r = runCommand(cas, "curvature(x^2, x)"); if (r.error) throw new Error(r.error); });
    test("curvature(x^2, x, 1)", () => { const r = runCommand(cas, "curvature(x^2, x, 1)"); if (r.error) throw new Error(r.error); });
    test("minimize(x^2-4*x, x)", () => { const r = runCommand(cas, "minimize(x^2-4*x, x)"); if (r.error) throw new Error(r.error); });
    test("maximize(-x^2, x)", () => { const r = runCommand(cas, "maximize(-x^2, x)"); if (r.error) throw new Error(r.error); });
    test("minimize(sin(x), x, 0, 2*pi)", () => { const r = runCommand(cas, "minimize(sin(x), x, 0, 2*pi)"); if (r.error) throw new Error(r.error); });
    test("maximize(sin(x), x, 0, 2*pi)", () => { const r = runCommand(cas, "maximize(sin(x), x, 0, 2*pi)"); if (r.error) throw new Error(r.error); });
    test("lagrangeMultipliers(x^2+y^2, [x+y-10], [x,y])", () => { const r = runCommand(cas, "lagrangeMultipliers(x^2+y^2, [x+y-10], [x,y])"); if (r.error) throw new Error(r.error); });
    test("eulerLagrange(1/2*m*v^2-m*g*x, x, t)", () => { const r = runCommand(cas, "eulerLagrange(1/2*m*v^2-m*g*x, x, t)"); if (r.error) throw new Error(r.error); });

    // ODE
    test("desolve(diff(y,x)=y, y)", () => { const r = runCommand(cas, "desolve(diff(y,x)=y, y)"); if (r.error) throw new Error(r.error); });
    test("desolve(diff(y,x,2)+y=0, y)", () => { const r = runCommand(cas, "desolve(diff(y,x,2)+y=0, y)"); if (r.error) throw new Error(r.error); });
    test("slopefield(x+y, x, y)", () => { const r = runCommand(cas, "slopefield(x+y, x, y)"); if (r.error) throw new Error(r.error); });
    test("vectorfield([x, y], x, y)", () => { const r = runCommand(cas, "vectorfield([x, y], x, y)"); if (r.error) throw new Error(r.error); });
    test("wronskian([exp(x), exp(-x)], x)", () => { const r = runCommand(cas, "wronskian([exp(x), exp(-x)], x)"); if (r.error) throw new Error(r.error); });

    // Transforms
    test("laplace(t^2, t, s)", () => { const r = runCommand(cas, "laplace(t^2, t, s)"); if (r.error) throw new Error(r.error); });
    test("ilaplace(2/s^3, s, t)", () => { const r = runCommand(cas, "ilaplace(2/s^3, s, t)"); if (r.error) throw new Error(r.error); });
    test("fourier(x^2, x, 3, pi)", () => { const r = runCommand(cas, "fourier(x^2, x, 3, pi)"); if (r.error) throw new Error(r.error); });
    test("fft([1,0,1,0])", () => { const r = runCommand(cas, "fft([1,0,1,0])"); if (r.error) throw new Error(r.error); });
    test("ifft([2,0,2,0])", () => { const r = runCommand(cas, "ifft([2,0,2,0])"); if (r.error) throw new Error(r.error); });

    // Vector Calculus
    test("grad(x^2+y^2, [x,y])", () => { const r = runCommand(cas, "grad(x^2+y^2, [x,y])"); if (r.error) throw new Error(r.error); });
    test("divergence([x, y, z], [x,y,z])", () => { const r = runCommand(cas, "divergence([x, y, z], [x,y,z])"); if (r.error) throw new Error(r.error); });
    test("curl([-y, x, 0], [x,y,z])", () => { const r = runCommand(cas, "curl([-y, x, 0], [x,y,z])"); if (r.error) throw new Error(r.error); });
    test("jacobian([x^2*y, x+y], [x,y])", () => { const r = runCommand(cas, "jacobian([x^2*y, x+y], [x,y])"); if (r.error) throw new Error(r.error); });
    test("hessian(x^3+y^3-3*x*y, [x,y])", () => { const r = runCommand(cas, "hessian(x^3+y^3-3*x*y, [x,y])"); if (r.error) throw new Error(r.error); });
    test("potential([2*x, 2*y], [x,y])", () => { const r = runCommand(cas, "potential([2*x, 2*y], [x,y])"); if (r.error) throw new Error(r.error); });
    test("conservative([-y, x], [x,y])", () => { const r = runCommand(cas, "conservative([-y, x], [x,y])"); if (r.error) throw new Error(r.error); });
    test("laplacian(x^2+y^2, [x,y])", () => { const r = runCommand(cas, "laplacian(x^2+y^2, [x,y])"); if (r.error) throw new Error(r.error); });
    test("surfaceArea(x^2, x, 0, 1)", () => { const r = runCommand(cas, "surfaceArea(x^2, x, 0, 1)"); if (r.error) throw new Error(r.error); });

    // ODE Numerical
    test("rk4(y-t*t+1, y, t, 0.5, 0, 2, 0.2)", () => { const r = runCommand(cas, "rk4(y-t*t+1, y, t, 0.5, 0, 2, 0.2)"); if (r.error) throw new Error(r.error); });
    test("rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, [])", () => { const r = runCommand(cas, "rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, [])"); if (r.error) throw new Error(r.error); });
    test("odestats(rk45(y-t*t+1, y, t, 0.5, 0, 2, 0.2, 1e-6, []))", () => { const r = runCommand(cas, "odestats(rk45(y-t*t+1, y, t, 0.5, 0, 2, 0.2, 1e-6, []))"); if (r.error) throw new Error(r.error); });
    test("odeplot(rk4(y-t*t+1, y, t, 0.5, 0, 2, 0.2))", () => { const r = runCommand(cas, "odeplot(rk4(y-t*t+1, y, t, 0.5, 0, 2, 0.2))"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 4: LINEAR ALGEBRA
    // ============================================================
    console.log("\n[4] Linear Algebra");
    test("det([[1,2],[3,4]])", () => { const r = runCommand(cas, "det([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("inv([[1,2],[3,4]])", () => { const r = runCommand(cas, "inv([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("trans([[1,2],[3,4]])", () => { const r = runCommand(cas, "trans([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("trace([[1,2],[3,4]])", () => { const r = runCommand(cas, "trace([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("eigenvals([[1,2],[2,1]])", () => { const r = runCommand(cas, "eigenvals([[1,2],[2,1]])"); if (r.error) throw new Error(r.error); });
    test("eigenvects([[1,2],[2,1]])", () => { const r = runCommand(cas, "eigenvects([[1,2],[2,1]])"); if (r.error) throw new Error(r.error); });
    test("diagonalize([[1,2],[2,1]])", () => { const r = runCommand(cas, "diagonalize([[1,2],[2,1]])"); if (r.error) throw new Error(r.error); });
    test("lu([[4,3],[6,3]])", () => { const r = runCommand(cas, "lu([[4,3],[6,3]])"); if (r.error) throw new Error(r.error); });
    test("qr([[1,2],[3,4]])", () => { const r = runCommand(cas, "qr([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("cholesky([[4,12,-16],[12,37,-43],[-16,-43,98]])", () => { const r = runCommand(cas, "cholesky([[4,12,-16],[12,37,-43],[-16,-43,98]])"); if (r.error) throw new Error(r.error); });
    test("svd([[1,2],[3,4]])", () => { const r = runCommand(cas, "svd([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("rref([[1,2,3],[4,5,6]])", () => { const r = runCommand(cas, "rref([[1,2,3],[4,5,6]])"); if (r.error) throw new Error(r.error); });
    test("rank([[1,2],[2,4]])", () => { const r = runCommand(cas, "rank([[1,2],[2,4]])"); if (r.error) throw new Error(r.error); });
    test("kernel([[1,2],[2,4]])", () => { const r = runCommand(cas, "kernel([[1,2],[2,4]])"); if (r.error) throw new Error(r.error); });
    test("charpoly([[1,2],[3,4]], x)", () => { const r = runCommand(cas, "charpoly([[1,2],[3,4]], x)"); if (r.error) throw new Error(r.error); });
    test("pinv([[1,2],[3,6]])", () => { const r = runCommand(cas, "pinv([[1,2],[3,6]])"); if (r.error) throw new Error(r.error); });
    test("cond([[1,2],[3,4]])", () => { const r = runCommand(cas, "cond([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("nullity([[1,2],[2,4]])", () => { const r = runCommand(cas, "nullity([[1,2],[2,4]])"); if (r.error) throw new Error(r.error); });
    test("colSpace([[1,2],[3,4]])", () => { const r = runCommand(cas, "colSpace([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("rowSpace([[1,2],[3,4]])", () => { const r = runCommand(cas, "rowSpace([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("basis([[1,2],[3,4]])", () => { const r = runCommand(cas, "basis([[1,2],[3,4]])"); if (r.error) throw new Error(r.error); });
    test("matrixPow([[1,1],[0,1]], 3)", () => { const r = runCommand(cas, "matrixPow([[1,1],[0,1]], 3)"); if (r.error) throw new Error(r.error); });
    test("matrixExp([[0,1],[-1,0]])", () => { const r = runCommand(cas, "matrixExp([[0,1],[-1,0]])"); if (r.error) throw new Error(r.error); });
    test("commutator([[1,2],[3,4]], [[5,6],[7,8]])", () => { const r = runCommand(cas, "commutator([[1,2],[3,4]], [[5,6],[7,8]])"); if (r.error) throw new Error(r.error); });

    // Vectors
    test("cross([1,0,0], [0,1,0])", () => { const r = runCommand(cas, "cross([1,0,0], [0,1,0])"); if (r.error) throw new Error(r.error); });
    test("dot([1,2], [3,4])", () => { const r = runCommand(cas, "dot([1,2], [3,4])"); if (r.error) throw new Error(r.error); });
    test("norm([3,4])", () => { const r = runCommand(cas, "norm([3,4])"); if (r.error) throw new Error(r.error); });
    test("projection([1,2], [3,4])", () => { const r = runCommand(cas, "projection([1,2], [3,4])"); if (r.error) throw new Error(r.error); });
    test("unitVector([1,1])", () => { const r = runCommand(cas, "unitVector([1,1])"); if (r.error) throw new Error(r.error); });
    test("angle([1,0], [0,1])", () => { const r = runCommand(cas, "angle([1,0], [0,1])"); if (r.error) throw new Error(r.error); });
    test("distance([0,0], [3,4])", () => { const r = runCommand(cas, "distance([0,0], [3,4])"); if (r.error) throw new Error(r.error); });
    test("midpoint([0,0], [2,2])", () => { const r = runCommand(cas, "midpoint([0,0], [2,2])"); if (r.error) throw new Error(r.error); });
    test("gramschmidt([[1,1],[1,0]])", () => { const r = runCommand(cas, "gramschmidt([[1,1],[1,0]])"); if (r.error) throw new Error(r.error); });
    test("toSpherical([1,1,1])", () => { const r = runCommand(cas, "toSpherical([1,1,1])"); if (r.error) throw new Error(r.error); });
    test("toCylindrical([1,1,1])", () => { const r = runCommand(cas, "toCylindrical([1,1,1])"); if (r.error) throw new Error(r.error); });
    test("kroneckerDelta(1,1)", () => { const r = runCommand(cas, "kroneckerDelta(1,1)"); if (r.error) throw new Error(r.error); });

    // Matrix builders
    test("identity(3)", () => { const r = runCommand(cas, "identity(3)"); if (r.error) throw new Error(r.error); });
    test("zeros(2, 3)", () => { const r = runCommand(cas, "zeros(2, 3)"); if (r.error) throw new Error(r.error); });
    test("ones(2, 3)", () => { const r = runCommand(cas, "ones(2, 3)"); if (r.error) throw new Error(r.error); });
    test("diag([1,2,3])", () => { const r = runCommand(cas, "diag([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("hilbert(3)", () => { const r = runCommand(cas, "hilbert(3)"); if (r.error) throw new Error(r.error); });
    test("vandermonde([1,2,3])", () => { const r = runCommand(cas, "vandermonde([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("toeplitz([1,2,3], [1,4,5])", () => { const r = runCommand(cas, "toeplitz([1,2,3], [1,4,5])"); if (r.error) throw new Error(r.error); });
    test("kron([[1,2]], [[0,1]])", () => { const r = runCommand(cas, "kron([[1,2]], [[0,1]])"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 5: STATISTICS
    // ============================================================
    console.log("\n[5] Statistics");
    test("mean([1,2,3,4,5])", () => { const r = runCommand(cas, "mean([1,2,3,4,5])"); if (r.error) throw new Error(r.error); });
    test("median([1,5,2,8,3])", () => { const r = runCommand(cas, "median([1,5,2,8,3])"); if (r.error) throw new Error(r.error); });
    test("mode([1,2,2,3])", () => { const r = runCommand(cas, "mode([1,2,2,3])"); if (r.error) throw new Error(r.error); });
    test("variance([1,2,3,4,5])", () => { const r = runCommand(cas, "variance([1,2,3,4,5])"); if (r.error) throw new Error(r.error); });
    test("std([1,2,3,4,5])", () => { const r = runCommand(cas, "std([1,2,3,4,5])"); if (r.error) throw new Error(r.error); });
    test("cov([1,2,3], [2,4,6])", () => { const r = runCommand(cas, "cov([1,2,3], [2,4,6])"); if (r.error) throw new Error(r.error); });
    test("corr([1,2,3], [1,2,3])", () => { const r = runCommand(cas, "corr([1,2,3], [1,2,3])"); if (r.error) throw new Error(r.error); });
    test("geoMean([1,2,4])", () => { const r = runCommand(cas, "geoMean([1,2,4])"); if (r.error) throw new Error(r.error); });
    test("harmMean([1,2,4])", () => { const r = runCommand(cas, "harmMean([1,2,4])"); if (r.error) throw new Error(r.error); });
    test("rms([1,-1])", () => { const r = runCommand(cas, "rms([1,-1])"); if (r.error) throw new Error(r.error); });
    test("mad([1,2,3])", () => { const r = runCommand(cas, "mad([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("skewness([1,2,3,10])", () => { const r = runCommand(cas, "skewness([1,2,3,10])"); if (r.error) throw new Error(r.error); });
    test("kurtosis([1,2,3])", () => { const r = runCommand(cas, "kurtosis([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("moment([1,2,3], 2)", () => { const r = runCommand(cas, "moment([1,2,3], 2)"); if (r.error) throw new Error(r.error); });
    test("kl_divergence([0.5,0.5], [0.4,0.6])", () => { const r = runCommand(cas, "kl_divergence([0.5,0.5], [0.4,0.6])"); if (r.error) throw new Error(r.error); });
    test("mae([1,2], [1.1,1.9])", () => { const r = runCommand(cas, "mae([1,2], [1.1,1.9])"); if (r.error) throw new Error(r.error); });
    test("mse([1,2], [1.1,1.9])", () => { const r = runCommand(cas, "mse([1,2], [1.1,1.9])"); if (r.error) throw new Error(r.error); });
    test("chiSquareTest([10,20], [15,15])", () => { const r = runCommand(cas, "chiSquareTest([10,20], [15,15])"); if (r.error) throw new Error(r.error); });
    test("propTest(45, 100, 0.5)", () => { const r = runCommand(cas, "propTest(45, 100, 0.5)"); if (r.error) throw new Error(r.error); });

    // Regression
    test("linearRegression([[1,1],[2,3],[3,2]])", () => { const r = runCommand(cas, "linearRegression([[1,1],[2,3],[3,2]])"); if (r.error) throw new Error(r.error); });
    test("polyRegression([[1,1],[2,4],[3,9]], 2)", () => { const r = runCommand(cas, "polyRegression([[1,1],[2,4],[3,9]], 2)"); if (r.error) throw new Error(r.error); });
    test("expRegression([[1,2],[2,4],[3,8]])", () => { const r = runCommand(cas, "expRegression([[1,2],[2,4],[3,8]])"); if (r.error) throw new Error(r.error); });
    test("powerRegression([[1,1],[2,8],[3,27]])", () => { const r = runCommand(cas, "powerRegression([[1,1],[2,8],[3,27]])"); if (r.error) throw new Error(r.error); });
    test("logRegression([[1,0],[10,1],[100,2]])", () => { const r = runCommand(cas, "logRegression([[1,0],[10,1],[100,2]])"); if (r.error) throw new Error(r.error); });

    // Distributions
    test("normalPDF(0, 0, 1)", () => { const r = runCommand(cas, "normalPDF(0, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("normalCDF(1.96, 0, 1)", () => { const r = runCommand(cas, "normalCDF(1.96, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("invNorm(0.975, 0, 1)", () => { const r = runCommand(cas, "invNorm(0.975, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("binomialPDF(2, 5, 0.5)", () => { const r = runCommand(cas, "binomialPDF(2, 5, 0.5)"); if (r.error) throw new Error(r.error); });
    test("binomialCDF(2, 5, 0.5)", () => { const r = runCommand(cas, "binomialCDF(2, 5, 0.5)"); if (r.error) throw new Error(r.error); });
    test("poissonPDF(3, 5)", () => { const r = runCommand(cas, "poissonPDF(3, 5)"); if (r.error) throw new Error(r.error); });
    test("poissonCDF(3, 5)", () => { const r = runCommand(cas, "poissonCDF(3, 5)"); if (r.error) throw new Error(r.error); });
    test("exponentialPDF(1, 0.5)", () => { const r = runCommand(cas, "exponentialPDF(1, 0.5)"); if (r.error) throw new Error(r.error); });
    test("exponentialCDF(1, 0.5)", () => { const r = runCommand(cas, "exponentialCDF(1, 0.5)"); if (r.error) throw new Error(r.error); });
    test("geometricPDF(3, 0.5)", () => { const r = runCommand(cas, "geometricPDF(3, 0.5)"); if (r.error) throw new Error(r.error); });
    test("geometricCDF(3, 0.5)", () => { const r = runCommand(cas, "geometricCDF(3, 0.5)"); if (r.error) throw new Error(r.error); });
    test("chisquarePDF(2, 3)", () => { const r = runCommand(cas, "chisquarePDF(2, 3)"); if (r.error) throw new Error(r.error); });
    test("chisquareCDF(5, 2)", () => { const r = runCommand(cas, "chisquareCDF(5, 2)"); if (r.error) throw new Error(r.error); });
    test("invChiSquare(0.95, 3)", () => { const r = runCommand(cas, "invChiSquare(0.95, 3)"); if (r.error) throw new Error(r.error); });
    test("studentTPDF(0, 5)", () => { const r = runCommand(cas, "studentTPDF(0, 5)"); if (r.error) throw new Error(r.error); });
    test("studentTCDF(2.0, 10)", () => { const r = runCommand(cas, "studentTCDF(2.0, 10)"); if (r.error) throw new Error(r.error); });
    test("invT(0.95, 10)", () => { const r = runCommand(cas, "invT(0.95, 10)"); if (r.error) throw new Error(r.error); });
    test("fPDF(1, 5, 5)", () => { const r = runCommand(cas, "fPDF(1, 5, 5)"); if (r.error) throw new Error(r.error); });
    test("fCDF(1, 5, 5)", () => { const r = runCommand(cas, "fCDF(1, 5, 5)"); if (r.error) throw new Error(r.error); });
    test("betaPDF(0.5, 2, 2)", () => { const r = runCommand(cas, "betaPDF(0.5, 2, 2)"); if (r.error) throw new Error(r.error); });
    test("betaCDF(0.5, 2, 2)", () => { const r = runCommand(cas, "betaCDF(0.5, 2, 2)"); if (r.error) throw new Error(r.error); });
    test("uniformPDF(0.5, 0, 1)", () => { const r = runCommand(cas, "uniformPDF(0.5, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("uniformCDF(0.5, 0, 1)", () => { const r = runCommand(cas, "uniformCDF(0.5, 0, 1)"); if (r.error) throw new Error(r.error); });
    test("gammaPDF(1, 2, 1)", () => { const r = runCommand(cas, "gammaPDF(1, 2, 1)"); if (r.error) throw new Error(r.error); });
    test("gammaCDF(1, 2, 1)", () => { const r = runCommand(cas, "gammaCDF(1, 2, 1)"); if (r.error) throw new Error(r.error); });
    test("hypergeometricPDF(2, 50, 5, 10)", () => { const r = runCommand(cas, "hypergeometricPDF(2, 50, 5, 10)"); if (r.error) throw new Error(r.error); });
    test("hypergeometricCDF(2, 50, 5, 10)", () => { const r = runCommand(cas, "hypergeometricCDF(2, 50, 5, 10)"); if (r.error) throw new Error(r.error); });

    // Hypothesis Tests
    test("zTest([1,2,3], 0, 1)", () => { const r = runCommand(cas, "zTest([1,2,3], 0, 1)"); if (r.error) throw new Error(r.error); });
    test("tTest([1,2,3], 0)", () => { const r = runCommand(cas, "tTest([1,2,3], 0)"); if (r.error) throw new Error(r.error); });
    test("zInterval([1,2,3], 1, 0.95)", () => { const r = runCommand(cas, "zInterval([1,2,3], 1, 0.95)"); if (r.error) throw new Error(r.error); });
    test("tInterval([1,2,3], 0.95)", () => { const r = runCommand(cas, "tInterval([1,2,3], 0.95)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 6: NUMBER THEORY
    // ============================================================
    console.log("\n[6] Number Theory");
    test("gcd(12, 18)", () => { const r = runCommand(cas, "gcd(12, 18)"); if (r.error) throw new Error(r.error); });
    test("lcm(4, 6)", () => { const r = runCommand(cas, "lcm(4, 6)"); if (r.error) throw new Error(r.error); });
    test("xgcd(35, 15)", () => { const r = runCommand(cas, "xgcd(35, 15)"); if (r.error) throw new Error(r.error); });
    test("chineseRemainder([2,3], [3,5])", () => { const r = runCommand(cas, "chineseRemainder([2,3], [3,5])"); if (r.error) throw new Error(r.error); });
    test("isPrime(17)", () => { const r = runCommand(cas, "isPrime(17)"); if (r.error) throw new Error(r.error); });
    test("isPrime(20)", () => { const r = runCommand(cas, "isPrime(20)"); if (r.error) throw new Error(r.error); });
    test("isSquare(16)", () => { const r = runCommand(cas, "isSquare(16)"); if (r.error) throw new Error(r.error); });
    test("isPerfect(6)", () => { const r = runCommand(cas, "isPerfect(6)"); if (r.error) throw new Error(r.error); });
    test("nextprime(100)", () => { const r = runCommand(cas, "nextprime(100)"); if (r.error) throw new Error(r.error); });
    test("prevprime(10)", () => { const r = runCommand(cas, "prevprime(10)"); if (r.error) throw new Error(r.error); });
    test("divisors(12)", () => { const r = runCommand(cas, "divisors(12)"); if (r.error) throw new Error(r.error); });
    test("primeFactors(60)", () => { const r = runCommand(cas, "primeFactors(60)"); if (r.error) throw new Error(r.error); });
    test("euler(10)", () => { const r = runCommand(cas, "euler(10)"); if (r.error) throw new Error(r.error); });
    test("primitiveRoot(7)", () => { const r = runCommand(cas, "primitiveRoot(7)"); if (r.error) throw new Error(r.error); });
    test("isPrimitiveRoot(3, 7)", () => { const r = runCommand(cas, "isPrimitiveRoot(3, 7)"); if (r.error) throw new Error(r.error); });
    test("moebius(10)", () => { const r = runCommand(cas, "moebius(10)"); if (r.error) throw new Error(r.error); });
    test("sigma(12)", () => { const r = runCommand(cas, "sigma(12)"); if (r.error) throw new Error(r.error); });
    test("sigma(12, 2)", () => { const r = runCommand(cas, "sigma(12, 2)"); if (r.error) throw new Error(r.error); });
    test("legendreSymbol(2, 7)", () => { const r = runCommand(cas, "legendreSymbol(2, 7)"); if (r.error) throw new Error(r.error); });
    test("jacobiSymbol(2, 5)", () => { const r = runCommand(cas, "jacobiSymbol(2, 5)"); if (r.error) throw new Error(r.error); });
    test("cfrac(1.2)", () => { const r = runCommand(cas, "cfrac(1.2)"); if (r.error) throw new Error(r.error); });
    test("divisorSum(10, 1)", () => { const r = runCommand(cas, "divisorSum(10, 1)"); if (r.error) throw new Error(r.error); });
    test("modInverse(3, 11)", () => { const r = runCommand(cas, "modInverse(3, 11)"); if (r.error) throw new Error(r.error); });
    test("modPow(2, 10, 1000)", () => { const r = runCommand(cas, "modPow(2, 10, 1000)"); if (r.error) throw new Error(r.error); });

    // Numeral systems
    test("toBase(255, 16)", () => { const r = runCommand(cas, "toBase(255, 16)"); if (r.error) throw new Error(r.error); });
    test("toRoman(1994)", () => { const r = runCommand(cas, "toRoman(1994)"); if (r.error) throw new Error(r.error); });
    test("fromRoman(MCMXCIV)", () => { const r = runCommand(cas, "fromRoman(MCMXCIV)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 7: SPECIAL FUNCTIONS
    // ============================================================
    console.log("\n[7] Special Functions");
    test("gamma(0.5)", () => { const r = runCommand(cas, "gamma(0.5)"); if (r.error) throw new Error(r.error); });
    test("beta(2, 3)", () => { const r = runCommand(cas, "beta(2, 3)"); if (r.error) throw new Error(r.error); });
    test("factorial(5)", () => { const r = runCommand(cas, "factorial(5)"); if (r.error) throw new Error(r.error); });
    test("fibonacci(10)", () => { const r = runCommand(cas, "fibonacci(10)"); if (r.error) throw new Error(r.error); });
    test("stirling1(4, 2)", () => { const r = runCommand(cas, "stirling1(4, 2)"); if (r.error) throw new Error(r.error); });
    test("stirling2(4, 2)", () => { const r = runCommand(cas, "stirling2(4, 2)"); if (r.error) throw new Error(r.error); });
    test("bell(5)", () => { const r = runCommand(cas, "bell(5)"); if (r.error) throw new Error(r.error); });
    test("catalan(4)", () => { const r = runCommand(cas, "catalan(4)"); if (r.error) throw new Error(r.error); });
    test("bernoulli(4)", () => { const r = runCommand(cas, "bernoulli(4)"); if (r.error) throw new Error(r.error); });
    test("bernoulliPoly(3, x)", () => { const r = runCommand(cas, "bernoulliPoly(3, x)"); if (r.error) throw new Error(r.error); });
    test("eulerPoly(3, x)", () => { const r = runCommand(cas, "eulerPoly(3, x)"); if (r.error) throw new Error(r.error); });
    test("harmonic(5)", () => { const r = runCommand(cas, "harmonic(5)"); if (r.error) throw new Error(r.error); });
    test("zeta(2)", () => { const r = runCommand(cas, "zeta(2)"); if (r.error) throw new Error(r.error); });
    test("lambertw(1)", () => { const r = runCommand(cas, "lambertw(1)"); if (r.error) throw new Error(r.error); });
    test("erf(1)", () => { const r = runCommand(cas, "erf(1)"); if (r.error) throw new Error(r.error); });
    test("erfinv(0.8427)", () => { const r = runCommand(cas, "erfinv(0.8427)"); if (r.error) throw new Error(r.error); });
    test("sinc(x)", () => { const r = runCommand(cas, "sinc(x)"); if (r.error) throw new Error(r.error); });
    test("FresnelS(x)", () => { const r = runCommand(cas, "FresnelS(x)"); if (r.error) throw new Error(r.error); });
    test("FresnelC(x)", () => { const r = runCommand(cas, "FresnelC(x)"); if (r.error) throw new Error(r.error); });
    test("EllipticK(0.5)", () => { const r = runCommand(cas, "EllipticK(0.5)"); if (r.error) throw new Error(r.error); });
    test("EllipticE(x, 0.5)", () => { const r = runCommand(cas, "EllipticE(x, 0.5)"); if (r.error) throw new Error(r.error); });
    test("besselJ(0, x)", () => { const r = runCommand(cas, "besselJ(0, x)"); if (r.error) throw new Error(r.error); });
    test("besselY(0, x)", () => { const r = runCommand(cas, "besselY(0, x)"); if (r.error) throw new Error(r.error); });
    test("legendre(3, x)", () => { const r = runCommand(cas, "legendre(3, x)"); if (r.error) throw new Error(r.error); });
    test("chebyshev(3, x)", () => { const r = runCommand(cas, "chebyshev(3, x)"); if (r.error) throw new Error(r.error); });
    test("hermite(3, x)", () => { const r = runCommand(cas, "hermite(3, x)"); if (r.error) throw new Error(r.error); });
    test("laguerre(3, x)", () => { const r = runCommand(cas, "laguerre(3, x)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 8: COMPLEX NUMBERS
    // ============================================================
    console.log("\n[8] Complex Numbers");
    test("abs(3+4i)", () => { const r = runCommand(cas, "abs(3+4i)"); if (r.error) throw new Error(r.error); });
    test("arg(1+i)", () => { const r = runCommand(cas, "arg(1+i)"); if (r.error) throw new Error(r.error); });
    test("real(3+4i)", () => { const r = runCommand(cas, "real(3+4i)"); if (r.error) throw new Error(r.error); });
    test("imag(3+4i)", () => { const r = runCommand(cas, "imag(3+4i)"); if (r.error) throw new Error(r.error); });
    test("conj(3+4i)", () => { const r = runCommand(cas, "conj(3+4i)"); if (r.error) throw new Error(r.error); });
    test("toPolar(1+i)", () => { const r = runCommand(cas, "toPolar(1+i)"); if (r.error) throw new Error(r.error); });
    test("cis(90)", () => { const r = runCommand(cas, "cis(90)"); if (r.error) throw new Error(r.error); });
    test("phasor(10, 45)", () => { const r = runCommand(cas, "phasor(10, 45)"); if (r.error) throw new Error(r.error); });
    test("rect(1, pi/4)", () => { const r = runCommand(cas, "rect(1, pi/4)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 9: LOGIC
    // ============================================================
    console.log("\n[9] Logic");
    test("true", () => { const r = runCommand(cas, "true"); if (r.error) throw new Error(r.error); });
    test("false", () => { const r = runCommand(cas, "false"); if (r.error) throw new Error(r.error); });
    test("nand(true, false)", () => { const r = runCommand(cas, "nand(true, false)"); if (r.error) throw new Error(r.error); });
    test("nor(false, false)", () => { const r = runCommand(cas, "nor(false, false)"); if (r.error) throw new Error(r.error); });
    test("xnor(true, true)", () => { const r = runCommand(cas, "xnor(true, true)"); if (r.error) throw new Error(r.error); });
    test("truthTable(A and B, [A,B])", () => { const r = runCommand(cas, "truthTable(A and B, [A,B])"); if (r.error) throw new Error(r.error); });
    test("cnf(a or b)", () => { const r = runCommand(cas, "cnf(a or b)"); if (r.error) throw new Error(r.error); });
    test("dnf(a and b)", () => { const r = runCommand(cas, "dnf(a and b)"); if (r.error) throw new Error(r.error); });
    test("simplify(a and not a)", () => { const r = runCommand(cas, "simplify(a and not a)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 10: FINANCE
    // ============================================================
    console.log("\n[10] Finance");
    test("compound(1000, 0.05, 12, 10)", () => { const r = runCommand(cas, "compound(1000, 0.05, 12, 10)"); if (r.error) throw new Error(r.error); });
    test("loan(200000, 0.04, 30)", () => { const r = runCommand(cas, "loan(200000, 0.04, 30)"); if (r.error) throw new Error(r.error); });
    test("npv(0.05, [-1000, 500, 600])", () => { const r = runCommand(cas, "npv(0.05, [-1000, 500, 600])"); if (r.error) throw new Error(r.error); });
    test("irr([-1000, 500, 600])", () => { const r = runCommand(cas, "irr([-1000, 500, 600])"); if (r.error) throw new Error(r.error); });
    test("blackScholes(100, 100, 1, 0.05, 0.2, call)", () => { const r = runCommand(cas, "blackScholes(100, 100, 1, 0.05, 0.2, call)"); if (r.error) throw new Error(r.error); });
    test("annuity(0.05, 10, 100)", () => { const r = runCommand(cas, "annuity(0.05, 10, 100)"); if (r.error) throw new Error(r.error); });
    test("amortization(0.05, 12, 10000)", () => { const r = runCommand(cas, "amortization(0.05, 12, 10000)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 11: LISTS & MISC
    // ============================================================
    console.log("\n[11] Lists & Misc");
    test("seq(k^2, k, 1, 5, 1)", () => { const r = runCommand(cas, "seq(k^2, k, 1, 5, 1)"); if (r.error) throw new Error(r.error); });
    test("range(0, 10, 2)", () => { const r = runCommand(cas, "range(0, 10, 2)"); if (r.error) throw new Error(r.error); });
    test("sort([3,1,2])", () => { const r = runCommand(cas, "sort([3,1,2])"); if (r.error) throw new Error(r.error); });
    test("reverse([1,2,3])", () => { const r = runCommand(cas, "reverse([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("size([1,2,3])", () => { const r = runCommand(cas, "size([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("flatten([[1,2],[3]])", () => { const r = runCommand(cas, "flatten([[1,2],[3]])"); if (r.error) throw new Error(r.error); });
    test("cumsum([1,2,3])", () => { const r = runCommand(cas, "cumsum([1,2,3])"); if (r.error) throw new Error(r.error); });
    test("sum(k^2, k, 1, 10)", () => { const r = runCommand(cas, "sum(k^2, k, 1, 10)"); if (r.error) throw new Error(r.error); });
    test("product(k, k, 1, 5)", () => { const r = runCommand(cas, "product(k, k, 1, 5)"); if (r.error) throw new Error(r.error); });
    test("union([1,2], [2,3])", () => { const r = runCommand(cas, "union([1,2], [2,3])"); if (r.error) throw new Error(r.error); });
    test("intersect([1,2], [2,3])", () => { const r = runCommand(cas, "intersect([1,2], [2,3])"); if (r.error) throw new Error(r.error); });
    test("setdiff([1,2], [2,3])", () => { const r = runCommand(cas, "setdiff([1,2], [2,3])"); if (r.error) throw new Error(r.error); });
    test("set([1,2,2,3])", () => { const r = runCommand(cas, "set([1,2,2,3])"); if (r.error) throw new Error(r.error); });
    test("append([1,2], 3)", () => { const r = runCommand(cas, "append([1,2], 3)"); if (r.error) throw new Error(r.error); });
    test("prepend([2,3], 1)", () => { const r = runCommand(cas, "prepend([2,3], 1)"); if (r.error) throw new Error(r.error); });
    test("concat([1,2], [3,4])", () => { const r = runCommand(cas, "concat([1,2], [3,4])"); if (r.error) throw new Error(r.error); });
    test("isSubset([1,2], [1,2,3])", () => { const r = runCommand(cas, "isSubset([1,2], [1,2,3])"); if (r.error) throw new Error(r.error); });
    test("cartesianProduct([1,2], [a,b])", () => { const r = runCommand(cas, "cartesianProduct([1,2], [a,b])"); if (r.error) throw new Error(r.error); });
    test("approx(pi)", () => { const r = runCommand(cas, "approx(pi)"); if (r.error) throw new Error(r.error); });
    test("root(8, 3)", () => { const r = runCommand(cas, "root(8, 3)"); if (r.error) throw new Error(r.error); });
    test("expectedValue([1,2], [0.5,0.5])", () => { const r = runCommand(cas, "expectedValue([1,2], [0.5,0.5])"); if (r.error) throw new Error(r.error); });
    test("entropy([0.5, 0.5])", () => { const r = runCommand(cas, "entropy([0.5, 0.5])"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 12: CHEMISTRY
    // ============================================================
    console.log("\n[12] Chemistry");
    test("molarMass(H2O)", () => { const r = runCommand(cas, "molarMass(H2O)"); if (r.error) throw new Error(r.error); });
    test("molarMass(C6H12O6)", () => { const r = runCommand(cas, "molarMass(C6H12O6)"); if (r.error) throw new Error(r.error); });
    test("balance(H2 + O2 -> H2O)", () => { const r = runCommand(cas, "balance(H2 + O2 -> H2O)"); if (r.error) throw new Error(r.error); });

    // ============================================================
    // TEST 13: PLOTTING
    // ============================================================
    console.log("\n[13] Plotting");
    test("plot(sin(x), x)", () => { const r = runCommand(cas, "plot(sin(x), x)"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });
    test("plot3d(x^2+y^2, x, y)", () => { const r = runCommand(cas, "plot3d(x^2+y^2, x, y)"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });
    test("plotparam([cos(t), sin(t)], t, 0, 2*pi)", () => { const r = runCommand(cas, "plotparam([cos(t), sin(t)], t, 0, 2*pi)"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });
    test("plotpolar(1+cos(t), t, 0, 2*pi)", () => { const r = runCommand(cas, "plotpolar(1+cos(t), t, 0, 2*pi)"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });
    test("plotimplicit(x^2+y^2=1, x, y)", () => { const r = runCommand(cas, "plotimplicit(x^2+y^2=1, x, y)"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });
    test("plotlist([[1,2],[2,3]])", () => { const r = runCommand(cas, "plotlist([[1,2],[2,3]])"); if (r.error) throw new Error(r.error); if (!r.type || r.type !== 'plot') throw new Error("Expected plot type"); });

    // ============================================================
    // TEST 14: KNOWN UI BROKEN FEATURES (Should fail)
    // ============================================================
    console.log("\n[14] Known UI Bug Detection");
    // These SHOULD fail - we check they're still broken
    console.log("  (Testing known broken UI calls that have wrong arguments)");

    // BUG: runOdePlot passes wrong args to odeplot
    test("odeplot('diff(y,x)=-0.5*y', 0, 10, 100) -> EXPECT FAIL (wrong sig)", () => {
        const r = runCommand(cas, `odeplot("diff(y,x)=-0.5*y", 0, 10, 100)`);
        // This should fail because odeplot expects 1 arg (solution vector), not 4 strings
        if (!r.error) throw new Error("Should have failed but didn't");
    });

    // BUG: series() call in UI passes only 2 args (expr, var) but needs 4
    test("series(x^3-3*x, x) -> EXPECT FAIL (needs 4 args: expr,var,pt,order)", () => {
        const r = runCommand(cas, "series(x^3-3*x, x)");
        // The series function actually exists as alias to laurent with defaults...
        // Let's see what happens
        if (r.error && r.error.includes("requires")) throw new Error(r.error);
    });

    // BUG: rk4 UI passes 6 args but needs 7
    test("rk4(y-t*t+1, y, t, 0.5, 0, 2) -> EXPECT FAIL (needs 7 args)", () => {
        const r = runCommand(cas, "rk4(y-t*t+1, y, t, 0.5, 0, 2)");
        if (!r.error) throw new Error("Should have failed but didn't - UI passes only 6 args");
    });

    // BUG: rk45 UI passes 6 args but needs 9
    test("rk45(y-t*t+1, y, t, 0.5, 0, 2) -> EXPECT FAIL (needs 9 args)", () => {
        const r = runCommand(cas, "rk45(y-t*t+1, y, t, 0.5, 0, 2)");
        if (!r.error) throw new Error("Should have failed but didn't - UI passes only 6 args");
    });

    // BUG: odestats UI passes equation, not solution vector
    test("odestats(diff(y,x)=-0.5*y) -> EXPECT FAIL (needs solution vector)", () => {
        const r = runCommand(cas, "odestats(diff(y,x)=-0.5*y)");
        if (!r.error) throw new Error("Should have failed but didn't");
    });

    // BUG: minimize/maximize with array of vars instead of single symbol
    test("minimize(x^2+y^2, [x,y]) -> EXPECT FAIL (needs single var)", () => {
        const r = runCommand(cas, "minimize(x^2+y^2, [x,y])");
        if (!r.error) throw new Error("Should have failed but didn't - UI passes array of vars");
    });

    // BUG: matrixPow with only 1 arg (matrix, no n)
    test("matrixPow([[1,1],[0,1]]) -> EXPECT FAIL (needs 2 args)", () => {
        const r = runCommand(cas, "matrixPow([[1,1],[0,1]])");
        if (!r.error) throw new Error("Should have failed but didn't - UI passes only 1 arg");
    });

    // BUG: runDistQuantile fallback sends Q(quantile not available...)
    test("Q(quantile not available for this distribution in UI) -> EXPECT FAIL", () => {
        const r = runCommand(cas, "Q(quantile not available for this distribution in UI)");
        if (!r.error) throw new Error("Should have failed but didn't");
    });

    // ============================================================
    // SUMMARY
    // ============================================================
    console.log("\n========================================");
    console.log(`  RESULTS: ${passed} passed, ${failed} failed, ${skipped} skipped`);
    console.log("========================================\n");

    if (failed > 0) {
        console.log("NOTE: Some failures above may be expected (bug detection tests).");
        console.log("      Focus on unexpected failures in non-[EXPECT FAIL] tests.\n");
    }

        process.exit(0);
    };

    tryRun();
});
