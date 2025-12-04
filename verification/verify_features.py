
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load the local index.html
        cwd = os.getcwd()
        page.goto(f"file://{cwd}/index.html")

        # 1. Test 3D Plotting
        # Type "plot3d(x*y, x, y)" into the input and evaluate
        page.fill("#command-input", "plot3d(sin(x)*cos(y), x, y)")
        page.click("#submit-btn")

        # Wait for plot to render
        page.wait_for_timeout(1000)

        # Take screenshot of the 3D plot
        page.screenshot(path="verification/plot3d.png")

        # 2. Test Control Flow (For loop)
        # Type "s:=0; for(i:=1; i<=5; i:=i+1) { s:=s+i }; s"
        page.fill("#command-input", "s:=0; for(i:=1; i<=5; i:=i+1) { s:=s+i }; s")
        page.click("#submit-btn")

        page.wait_for_timeout(500)

        # Take screenshot of the result
        page.screenshot(path="verification/control_flow.png")

        browser.close()

if __name__ == "__main__":
    run()
