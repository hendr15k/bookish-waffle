const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({
    viewport: { width: 320, height: 568 },
    isMobile: true,
    hasTouch: true
  });
  const page = await context.newPage();

  await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  // Check what's intercepting the mode-toggle button
  const modeToggle = page.locator('#mode-toggle');
  
  // Use force click to bypass interception
  await modeToggle.click({ force: true });
  await page.waitForTimeout(500);

  // Check state
  const appGrid = await page.locator('#app-grid').isVisible();
  console.log('App grid visible after force click:', appGrid);

  // Check if app-categories is now visible
  const appCategories = await page.locator('#app-categories').isVisible();
  console.log('App categories visible:', appCategories);

  // Check bottom bar
  const bottomBar = page.locator('.mobile-bottom-bar');
  const bottomBarVisible = await bottomBar.isVisible();
  console.log('Bottom bar visible:', bottomBarVisible);
  if (bottomBarVisible) {
    const box = await bottomBar.boundingBox();
    console.log('Bottom bar position:', box);
  }

  // Check header
  const headerDiv = page.locator('body > div').first();
  const headerBox = await headerDiv.boundingBox();
  console.log('Header position:', headerBox);

  // Check what element is at the location of mode-toggle
  const modeToggleBox = await modeToggle.boundingBox();
  console.log('Mode toggle position:', modeToggleBox);

  await browser.close();
})();