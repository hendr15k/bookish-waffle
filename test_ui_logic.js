
const fs = require('fs');
const jsdom = require('jsdom');
const { JSDOM } = jsdom;

const html = fs.readFileSync('index.html', 'utf8');

// Load libraries code
const code1 = fs.readFileSync('js/expression.js', 'utf8');
const code2 = fs.readFileSync('js/parser.js', 'utf8');
const code3 = fs.readFileSync('js/chemistry.js', 'utf8');
const code4 = fs.readFileSync('js/cas.js', 'utf8');

const dom = new JSDOM(html, {
    runScripts: "dangerously",
    resources: "usable",
    url: "http://localhost/", // To allow localStorage
    beforeParse(window) {
        // Mock MathJax
        window.MathJax = { typesetPromise: () => Promise.resolve() };
    }
});
const window = dom.window;
const document = window.document;

// We need to execute the inline scripts in the HTML, which JSDOM does if runScripts is dangerously.
// However, the external scripts <script src="..."> won't load from disk easily.
// We manually execute the library code first.

const script1 = document.createElement('script'); script1.textContent = code1; document.head.appendChild(script1);
const script2 = document.createElement('script'); script2.textContent = code2; document.head.appendChild(script2);
const script3 = document.createElement('script'); script3.textContent = code3; document.head.appendChild(script3);
const script4 = document.createElement('script'); script4.textContent = code4; document.head.appendChild(script4);

// Wait for inline scripts to execute? They execute synchronously when parsing if they are inline.
// But we need to ensure the main inline script (the last one in index.html) runs.
// JSDOM executes inline scripts during parsing. But since we loaded libraries *after* parsing started (via appendChild),
// the inline script might have failed if it depended on them immediately?
// Actually, `new JSDOM(html)` parses everything. If src scripts fail to load, the global `CAS` class won't be defined.
// Then the inline script `const cas = new CAS();` will throw.

// Strategy: Modify HTML to embed scripts before passing to JSDOM?
// Or just eval the libraries in the window, THEN eval the inline script logic.

// Let's reload the DOM but inject script content into the HTML string first.
let htmlWithScripts = html;
htmlWithScripts = htmlWithScripts.replace('<script src="js/expression.js"></script>', `<script>${code1}</script>`);
htmlWithScripts = htmlWithScripts.replace('<script src="js/parser.js"></script>', `<script>${code2}</script>`);
htmlWithScripts = htmlWithScripts.replace('<script src="js/chemistry.js"></script>', `<script>${code3}</script>`);
htmlWithScripts = htmlWithScripts.replace('<script src="js/cas.js"></script>', `<script>${code4}</script>`);

const dom2 = new JSDOM(htmlWithScripts, {
    runScripts: "dangerously",
    resources: "usable",
    url: "http://localhost/"
});
const window2 = dom2.window;
const document2 = window2.document;

// Wait a tick for async stuff if any?
// Check functions
const check = (name) => {
    if (typeof window2[name] === 'function') {
        console.log(`[PASS] window.${name} is defined.`);
    } else {
        console.error(`[FAIL] window.${name} is not defined.`);
    }
};

check('generateMatrixGrid');
check('insertMatrix');
check('updateDistPlot');
check('runRiemann');
check('runTool');

// Verify labels
const inputs = document2.querySelectorAll('input.form-input');
let labelFailures = 0;
inputs.forEach(input => {
    const id = input.id;
    if (!id) return;
    const formGroup = input.closest('.form-group');
    if (formGroup) {
        const label = formGroup.querySelector('label.form-label');
        if (label) {
            if (label.getAttribute('for') !== id) {
                // console.error(`[FAIL] Label for input #${id} mismatch.`);
                labelFailures++;
            }
        }
    }
});
if (labelFailures === 0) console.log("[PASS] Labels checked.");

// Test logic
// Mock addToHistory to verify
let history = [];
window2.addToHistory = (cmd, res, isErr) => history.push({cmd, res, isErr});
window2.alert = (msg) => console.log("Alert:", msg);

try {
    window2.runTool('plot', ['', 'x', '', '']);
    if (history.length > 0 && history[0].isErr) {
        console.log("[PASS] runTool handles empty args.");
    } else {
        console.error("[FAIL] runTool did not error on empty args.");
    }
} catch(e) {
    console.error("[FAIL] Exception during runTool:", e);
}

// Cleanup
window2.close();
