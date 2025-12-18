
const fs = require('fs');
const jsdom = require('jsdom');
const { JSDOM } = jsdom;

const codeExpression = fs.readFileSync('js/expression.js', 'utf8');
const codeParser = fs.readFileSync('js/parser.js', 'utf8');
const codeChemistry = fs.readFileSync('js/chemistry.js', 'utf8');
const codeHelp = fs.readFileSync('js/help.js', 'utf8');
const codeCas = fs.readFileSync('js/cas.js', 'utf8');

let htmlContent = fs.readFileSync('index.html', 'utf8');

// Replace external script tags with inline code to ensure they load
htmlContent = htmlContent.replace('<script src="js/expression.js"></script>', `<script>${codeExpression}</script>`);
htmlContent = htmlContent.replace('<script src="js/parser.js"></script>', `<script>${codeParser}</script>`);
htmlContent = htmlContent.replace('<script src="js/chemistry.js"></script>', `<script>${codeChemistry}</script>`);
htmlContent = htmlContent.replace('<script src="js/help.js"></script>', `<script>${codeHelp}</script>`);
htmlContent = htmlContent.replace('<script src="js/cas.js"></script>', `<script>${codeCas}</script>`);

// Mock Canvas
const mockCanvasScript = `
<script>
  HTMLCanvasElement.prototype.getContext = function() {
    return {
      clearRect: () => {},
      beginPath: () => {},
      moveTo: () => {},
      lineTo: () => {},
      stroke: () => {},
      fill: () => {},
      arc: () => {},
      strokeRect: () => {},
      fillRect: () => {},
      setLineDash: () => {},
      createImageData: () => ({ data: [] }),
      putImageData: () => {},
      measureText: () => ({ width: 0 }),
      fillText: () => {}
    };
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

// Helper to run tests
function test(desc, fn) {
    try {
        fn();
        console.log(`[PASS] ${desc}`);
    } catch (e) {
        console.error(`[FAIL] ${desc}`);
        console.error(e);
    }
}

// Wait for load
window.addEventListener('load', () => {
    runTests();
});

function runTests() {
    console.log("Running UI Logic Tests...");

    test("CAS Initialization", () => {
        if (!window.cas) throw new Error("CAS object not initialized on window");
        if (!window.cas.evaluate) throw new Error("CAS.evaluate not found");
    });

    test("Labels checked", () => {
        const label = window.document.querySelector('label[for="calc-main-input"]');
        if (!label) throw new Error("Calculator input label missing");
    });

    test("Test runTool: Simplify", () => {
        const input = window.document.getElementById('calc-main-input');
        input.value = "1+1";

        window.runTool('simplify', ['calc-main-input']);

        // Wait for async update (runTool uses setTimeout 50ms)
        setTimeout(() => {
             const entries = window.document.querySelectorAll('.entry');
             const last = entries[entries.length - 1];
             if (!last) console.error("[FAIL] Simplify: No history entry");
             else if (!last.textContent.includes("2")) console.error("[FAIL] Simplify: Result 2 not found");
             else console.log("[PASS] Test runTool: Simplify");
        }, 100);
    });

    test("Test Matrix Builder", () => {
        const rows = window.document.getElementById('mat-rows');
        const cols = window.document.getElementById('mat-cols');
        rows.value = "2";
        cols.value = "2";

        window.generateMatrixGrid();

        const container = window.document.getElementById('matrix-builder-container');
        const inputs = container.querySelectorAll('input');
        if (inputs.length !== 4) throw new Error(`Expected 4 inputs, got ${inputs.length}`);

        inputs[0].value = "1";
        inputs[1].value = "2";
        inputs[2].value = "3";
        inputs[3].value = "4";

        window.insertMatrix();

        const res = window.document.getElementById('lin-expr').value;
        if (res !== "[[1,2],[3,4]]") throw new Error(`Matrix insertion failed. Got: ${res}`);
    });

    test("Test Converter Logic", () => {
        const fromVal = window.document.getElementById('conv-val-from');
        const cat = window.document.getElementById('conv-category');

        cat.value = 'length';
        window.updateConvUnits();

        const fromUnit = window.document.getElementById('conv-unit-from');
        const toUnit = window.document.getElementById('conv-unit-to');

        fromUnit.value = 'm';
        toUnit.value = 'km';

        fromVal.value = "1000";
        window.runConverter('from');

        const toVal = window.document.getElementById('conv-val-to').value;
        if (Math.abs(parseFloat(toVal) - 1) > 1e-5) throw new Error(`Converter failed: 1000m -> ${toVal}km`);
    });
}
