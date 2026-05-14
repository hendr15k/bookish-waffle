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

  // Scroll to top and check
  await page.evaluate(() => window.scrollTo(0, 0));
  await page.waitForTimeout(200);

  // Check all visible elements that could overlap the header
  const debugInfo = await page.evaluate(() => {
    const modeToggle = document.getElementById('mode-toggle');
    const rect = modeToggle.getBoundingClientRect();
    
    // Element at center of mode-toggle
    const centerX = rect.left + rect.width / 2;
    const centerY = rect.top + rect.height / 2;
    const elAtPoint = document.elementFromPoint(centerX, centerY);
    
    // Check if any fixed/absolute elements overlap
    const fixedElements = [];
    document.querySelectorAll('*').forEach(el => {
      const style = getComputedStyle(el);
      if ((style.position === 'fixed' || style.position === 'absolute') && 
          style.display !== 'none' && style.visibility !== 'hidden') {
        const r = el.getBoundingClientRect();
        if (r.top < rect.bottom && r.bottom > rect.top && r.left < rect.right && r.right > rect.left) {
          fixedElements.push({
            tag: el.tagName,
            id: el.id,
            classes: el.className.toString().substring(0, 80),
            rect: { top: r.top, left: r.left, width: r.width, height: r.height },
            zIndex: style.zIndex
          });
        }
      }
    });
    
    return {
      modeToggleRect: rect,
      elementAtPoint: elAtPoint ? {
        tag: elAtPoint.tagName,
        id: elAtPoint.id,
        classes: elAtPoint.className.toString().substring(0, 80),
        text: elAtPoint.textContent.substring(0, 30)
      } : null,
      fixedOverlaps: fixedElements
    };
  });

  console.log('Mode toggle rect:', JSON.stringify(debugInfo.modeToggleRect, null, 2));
  console.log('Element at point:', JSON.stringify(debugInfo.elementAtPoint, null, 2));
  console.log('Fixed elements overlapping:', JSON.stringify(debugInfo.fixedOverlaps, null, 2));

  await browser.close();
})();