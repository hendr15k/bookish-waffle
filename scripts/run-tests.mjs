// Test runner: runs every tests/test_*.js in order and fails the suite when a
// file exits non-zero OR reports a failure on its own output.
//
// Many of the older test files print "FAIL: ..." / "[FAIL] ..." without calling
// process.exit(1), so `npm test` used to stay green while a check was failing
// (test_new_features_23.js reported "FAIL: Check output" on every run).
// Lines that legitimately contain "EXPECT FAIL" (tests asserting that a command
// must throw) are ignored.
import { readdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import path from 'node:path';

const dir = 'tests';
const files = readdirSync(dir)
    .filter((f) => /^test_.*\.js$/.test(f))
    .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));

const failurePattern = /^\s*(\[FAIL\]|FAIL:|AssertionError|Uncaught)/;
const failing = [];

for (const file of files) {
    const full = path.join(dir, file);
    console.log('Running', full);
    const result = spawnSync(process.execPath, [full], { encoding: 'utf8' });
    const out = (result.stdout || '') + (result.stderr || '');
    process.stdout.write(out);

    const reported = out.split('\n').filter((line) => failurePattern.test(line) && !line.includes('EXPECT FAIL'));
    if (result.status !== 0) {
        failing.push(`${full} (exit ${result.status})`);
    } else if (reported.length > 0) {
        failing.push(`${full} (${reported.length} failed check${reported.length === 1 ? '' : 's'})`);
    }
}

if (failing.length > 0) {
    console.error(`\n${failing.length} test file(s) failed:`);
    for (const f of failing) console.error('  ' + f);
    process.exit(1);
}
console.log(`\nAll ${files.length} test files passed`);
