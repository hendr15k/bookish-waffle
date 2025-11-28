
from playwright.sync_api import sync_playwright
import os

def run():
    cwd = os.getcwd()
    file_url = f"file://{cwd}/index.html"

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(file_url)

        # Test diff(inv(x), x) which preserves inv as symbolic call
        page.fill("#command-input", "diff(inv(x), x)")
        page.click("#submit-btn")
        page.wait_for_timeout(1000)

        content = page.content()
        # \\left(x\\right)^{-1}
        # MathJax might modify DOM, but we check raw HTML before rendering or if rendering fails

        # Also check error message for inv(x)
        page.fill("#command-input", "inv(x)")
        page.click("#submit-btn")
        page.wait_for_timeout(1000)

        content = page.content()

        # We expect diff(inv(x), x) to produce LaTeX containing inverted format
        if "\\left(x\\right)^{-1}" in content:
             print("Verified: inv latex present in DOM")
        else:
             print("FAILED: inv latex NOT found in DOM")
             # Print nearby content for debugging
             print("Content snippet:")
             print(content[-500:])

        browser.close()

if __name__ == "__main__":
    run()
