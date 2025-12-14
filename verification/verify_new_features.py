
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local file
        filepath = os.path.abspath("index.html")
        page.goto(f"file://{filepath}")

        # 1. Switch to App Mode
        page.click("#mode-toggle")

        # 2. Open Statistics Tool
        page.click(".app-card:has-text('Statistics')")

        # 3. Switch to Distributions Tab
        page.click("button:has-text('Distributions')")

        # 4. Select Student-t
        page.select_option("#stat-dist-type", "studentt")

        # 5. Wait for canvas update (Student-t plot)
        page.wait_for_timeout(500)

        # Take screenshot of Statistics tool
        page.screenshot(path="verification/stats_studentt.png", clip={"x": 0, "y": 0, "width": 1000, "height": 800})

        # 6. Close tool - SPECIFIC for Statistics tool
        page.click("#tool-statistics .back-btn")

        # 7. Open Plotting Tool
        page.click(".app-card:has-text('Plotting')")

        # 8. Switch to 3D Plot Tab
        page.click("button:has-text('3D Plot')")

        # 9. Click Plot 3D (default: sin(sqrt(x^2+y^2)))
        page.click("button:has-text('Plot 3D')")

        # 10. Wait for plot
        page.wait_for_timeout(1000)

        # Take screenshot of Plotting tool
        page.screenshot(path="verification/plotting_3d.png", clip={"x": 0, "y": 0, "width": 1000, "height": 800})

        browser.close()

if __name__ == "__main__":
    run()
