const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  
  // Test on iPhone SE (320x568)
  const context = await browser.newContext({ viewport: { width: 320, height: 568 }, isMobile: true, hasTouch: true });
  const page = await context.newPage();
  await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  const info = await page.evaluate(() => {
    // Get the header container
    const headerContainer = document.querySelector('body > div:first-child');
    if (!headerContainer) return { error: 'no header container' };
    
    const rect = headerContainer.getBoundingClientRect();
    
    // Get all child elements and their positions
    const children = [];
    for (const child of headerContainer.children) {
      const childRect = child.getBoundingClientRect();
      children.push({
        tag: child.tagName,
        id: child.id || '',
        text: child.textContent.substring(0, 30),
        rect: { x: childRect.x, y: childRect.y, width: childRect.width, height: childRect.height }
      });
    }
    
    return {
      containerRect: { x: rect.x, y: rect.y, width: rect.width, height: rect.height },
      children: children,
      viewportHeight: window.innerHeight,
      scrollY: window.scrollY
    };
  });

  console.log('Header container:', JSON.stringify(info, null, 2));

  // Try to scroll mode-toggle into view and click
  await page.locator('#mode-toggle').scrollIntoViewIfNeeded();
  await page.waitForTimeout(500);
  await page.locator('#mode-toggle').click();
  await page.waitForTimeout(500);
  
  const appGrid = await page.locator('#app-grid').isVisible();
  console.log('App grid visible after scroll+click:', appGrid);

  await context.close();
  await browser.close();
})();