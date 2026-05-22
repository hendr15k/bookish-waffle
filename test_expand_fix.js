const fs = require('fs');
const { JSDOM } = require('jsdom');

const path = require('path');
const baseDir = __dirname;

const codeExpression = fs.readFileSync(path.join(baseDir, 'js/expression.js'), 'utf8');
const codeParser = fs.readFileSync(path.join(baseDir, 'js/parser.js'), 'utf8');
const codeChemistry = fs.readFileSync(path.join(baseDir, 'js/chemistry.js'), 'utf8');
const codeHelp = fs.readFileSync(path.join(baseDir, 'js/help.js'), 'utf8');
const codeCas = fs.readFileSync(path.join(baseDir, 'js/cas.js'), 'utf8');

let htmlContent = fs.readFileSync(path.join(baseDir, 'index.html'), 'utf8');

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

function runCommand(cas, Lexer, Parser, cmdStr) {
    try {
        const l = new Lexer(cmdStr);
        const p = new Parser(l);
        const t = p.parse();
        return cas.evaluate(t);
    } catch (e) {
        return { error: e.message };
    }
}

dom.window.addEventListener('load', () => {
    const cas = dom.window.cas;
    const Lexer = dom.window.Lexer;
    const Parser = dom.window.Parser;
    
    console.log("\n========================================");
    console.log("  Testing expand() simplification fix");
    console.log("========================================\n");

    // Test 1: expand((x+1)*(x-1)) should return x^2 - 1
    const test1 = runCommand(cas, Lexer, Parser, 'expand((x+1)*(x-1))');
    const result1 = typeof test1 === 'object' && test1 && test1.toString ? test1.toString() : String(test1);
    console.log('Test 1: expand((x+1)*(x-1))');
    console.log('Expected: (x^2 - 1)');
    console.log('Got:     ', result1);
    console.log('Status:  ', result1 === '(x^2 - 1)' ? 'PASS' : 'FAIL');
    console.log();

    // Test 2: expand((x+2)*(x-3)) should return x^2 - x - 6
    const test2 = runCommand(cas, Lexer, Parser, 'expand((x+2)*(x-3))');
    const result2 = typeof test2 === 'object' && test2 && test2.toString ? test2.toString() : String(test2);
    console.log('Test 2: expand((x+2)*(x-3))');
    console.log('Expected: (x^2 - x - 6)');
    console.log('Got:     ', result2);
    console.log('Status:  ', result2 === '(x^2 - x - 6)' ? 'PASS' : 'FAIL');
    console.log();

    // Test 3: expand((x+1)^2) should return (x^2 + (2 * x) + 1) simplified
    const test3 = runCommand(cas, Lexer, Parser, 'expand((x+1)^2)');
    const result3 = typeof test3 === 'object' && test3 && test3.toString ? test3.toString() : String(test3);
    console.log('Test 3: expand((x+1)^2)');
    console.log('Expected: simplified form (x^2 + 2x + 1)');
    console.log('Got:     ', result3);
    console.log();

    console.log("========================================\n");
    
    process.exit(0);
});
