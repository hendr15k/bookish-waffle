const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({
    viewport: { width: 375, height: 812 },
    isMobile: true,
    hasTouch: true
  });
  const page = await context.newPage();

  await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  // Switch to app mode
  await page.click('#mode-toggle');
  await page.waitForTimeout(500);

  // Open calculus
  await page.click('.app-card:has-text("Calculus")');
  await page.waitForTimeout(500);

  // Check all action grids
  const actionGrids = await page.locator('#tool-calculus .action-grid').count();
  console.log('Action grids found:', actionGrids);

  for (let i = 0; i < actionGrids; i++) {
    const grid = page.locator('#tool-calculus .action-grid').nth(i);
    const visible = await grid.isVisible();
    if (visible) {
      const gridCols = await grid.evaluate(el => getComputedStyle(el).gridTemplateColumns);
      const gridDisplay = await grid.evaluate(el => getComputedStyle(el).display);
      console.log(`Action grid ${i}: visible=${visible}, gridTemplateColumns=${gridCols}, display=${gridDisplay}`);
    }
  }

  // Check the exact CSS rule that's supposed to apply
  const cssCheck = await page.evaluate(() => {
    const grids = document.querySelectorAll('#tool-calculus .action-grid');
    const results = [];
    grids.forEach((g, i) => {
      const style = getComputedStyle(g);
      results.push({
        index: i,
        visible: g.offsetParent !== null,
        gridTemplateColumns: style.gridTemplateColumns,
        display: style.display
      });
    });
    return results;
  });
  console.log('CSS Check:', JSON.stringify(cssCheck, null, 2));

  // Also check function grids
  const funcGrids = await page.evaluate(() => {
    const grids = document.querySelectorAll('#tool-calculator .function-grid');
    const results = [];
    grids.forEach((g, i) => {
      const style = getComputedStyle(g);
      results.push({
        index: i,
        gridTemplateColumns: style.gridTemplateColumns
      });
    });
    return results;
  });
  console.log('Function grids:', JSON.stringify(funcGrids, null, 2));

  await browser.close();
})();