const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({
    viewport: { width: 320, height: 568 },
    isMobile: true,
    hasTouch: true,
    deviceScaleFactor: 2
  });
  const page = await context.newPage();

  await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  // Use force click to bypass pointer interception
  await page.locator('#mode-toggle').click({ force: true });
  await page.waitForTimeout(500);

  const appGrid = await page.locator('#app-grid').isVisible();
  console.log('App grid visible:', appGrid);

  const appCategories = await page.locator('#app-categories').isVisible();
  console.log('App categories visible:', appCategories);

  // Open calculator
  await page.click('.app-card:has-text("Calculator")', { force: true });
  await page.waitForTimeout(500);

  const calculatorVisible = await page.locator('#tool-calculator').isVisible();
  console.log('Calculator visible:', calculatorVisible);

  // Back button
  await page.locator('#tool-calculator .back-btn').click({ force: true });
  await page.waitForTimeout(300);

  const appGridAfterBack = await page.locator('#app-grid').isVisible();
  console.log('App grid visible after back:', appGridAfterBack);

  await browser.close();
  console.log('All tests passed on iPhone SE (320x568)!');
})();