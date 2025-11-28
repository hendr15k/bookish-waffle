
from playwright.sync_api import sync_playwright
import os

def run():
    # Get absolute path to index.html
    cwd = os.getcwd()
    file_url = f"file://{cwd}/index.html"

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()

        # Open the local file
        print(f"Opening {file_url}")
        page.goto(file_url)

        # Wait for page to load
        page.wait_for_selector("#command-input")

        # Function to evaluate and wait for result
        def evaluate(command):
            print(f" evaluating: {command}")
            page.fill("#command-input", command)
            page.click("#submit-btn")
            # Wait for result entry to appear (simplistic wait)
            page.wait_for_timeout(500)

        # Test problematic commands
        evaluate("asin(x)")

        # Give some time for MathJax (if it were loading)
        page.wait_for_timeout(2000)

        # Take screenshot
        screenshot_path = "verification/latex_fix_2.png"
        page.screenshot(path=screenshot_path, full_page=True)
        print(f"Screenshot saved to {screenshot_path}")

        # We can also check the innerHTML of the results to verify LaTeX strings are present
        content = page.content()
        if "\\arcsin" in content:
            print("Verified: \\arcsin present in DOM")
        else:
            print("FAILED: \\arcsin NOT found in DOM")

        browser.close()

if __name__ == "__main__":
    run()
