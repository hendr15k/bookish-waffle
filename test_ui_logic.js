
const fs = require('fs');
const jsdom = require('jsdom');
const { JSDOM } = jsdom;

const html = fs.readFileSync('index.html', 'utf8');

const dom = new JSDOM(html, {
    runScripts: "dangerously",
    resources: "usable",
    url: "http://localhost/" // To allow localStorage
});
const window = dom.window;
const document = window.document;

// Verify window.generateMatrixGrid is defined
if (typeof window.generateMatrixGrid !== 'function') {
    console.error("[FAIL] window.generateMatrixGrid is not defined.");
} else {
    console.log("[PASS] window.generateMatrixGrid is defined.");
}

// Verify window.insertMatrix is defined
if (typeof window.insertMatrix !== 'function') {
    console.error("[FAIL] window.insertMatrix is not defined.");
} else {
    console.log("[PASS] window.insertMatrix is defined.");
}

// Verify labels have 'for' attribute
const inputs = document.querySelectorAll('input.form-input');
let labelFailures = 0;

inputs.forEach(input => {
    const id = input.id;
    if (!id) return;

    // Find label in same form-group usually
    const formGroup = input.closest('.form-group');
    if (formGroup) {
        const label = formGroup.querySelector('label.form-label');
        if (label) {
            if (label.getAttribute('for') !== id) {
                console.error(`[FAIL] Label for input #${id} does not have correct 'for' attribute. Got '${label.getAttribute('for')}'.`);
                labelFailures++;
            }
        } else {
            // Maybe it is a sibling?
            // console.log(`[INFO] No label found in group for #${id}`);
        }
    }
});

if (labelFailures === 0) {
    console.log("[PASS] All labels checked have correct 'for' attributes.");
} else {
    console.error(`[FAIL] ${labelFailures} label issues found.`);
}

// Test runTool error handling (mock alert)
let historyEntries = [];
// Mock addToHistory
window.addToHistory = function(cmd, res, isErr) {
    historyEntries.push({cmd, res, isErr});
};
window.alert = function(msg) {
    console.error("[FAIL] alert() called: " + msg);
};

// Try to run tool with empty args
window.runTool('plot', ['', 'x', '', '']);
if (historyEntries.length > 0 && historyEntries[0].isErr && historyEntries[0].res.includes("fill in")) {
    console.log("[PASS] runTool handles empty args with addToHistory.");
} else {
    console.error("[FAIL] runTool did not log error to history for empty args.");
    console.log(historyEntries);
}

// Cleanup
window.close();
