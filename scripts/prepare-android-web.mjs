import { cpSync, existsSync, mkdirSync, rmSync, writeFileSync, copyFileSync } from 'node:fs';
import { join } from 'node:path';

const root = process.cwd();
const outDir = join(root, 'dist-android');

rmSync(outDir, { recursive: true, force: true });
mkdirSync(outDir, { recursive: true });
mkdirSync(join(outDir, 'js'), { recursive: true });

const filesToCopy = [
  'index.html',
  'manifest.json',
  'icon-192.png',
  'plot_save.png',
  'statistics_binomial.png'
];

for (const file of filesToCopy) {
  if (existsSync(join(root, file))) {
    copyFileSync(join(root, file), join(outDir, file));
  }
}

for (const file of ['cas.js', 'chemistry.js', 'expression.js', 'help.js', 'parser.js']) {
  copyFileSync(join(root, 'js', file), join(outDir, 'js', file));
}

if (!existsSync(join(outDir, 'icon-512.png')) && existsSync(join(outDir, 'icon-192.png'))) {
  copyFileSync(join(outDir, 'icon-192.png'), join(outDir, 'icon-512.png'));
}

writeFileSync(
  join(outDir, '.nojekyll'),
  ''
);

console.log('Prepared dist-android/ for Capacitor');
