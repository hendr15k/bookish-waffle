from playwright.sync_api import sync_playwright
import os
import time

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        page.on("console", lambda msg: print(f"Console: {msg.text}"))
        page.on("pageerror", lambda exc: print(f"Page Error: {exc}"))

        cwd = os.getcwd()
        file_path = f"file://{cwd}/index.html"
        page.goto(file_path)

        # Check if CAS is defined (main class from js/cas.js)
        try:
            page.evaluate("new CAS()")
            print("CAS initialized")
        except Exception as e:
            print(f"CAS init failed: {e}")

        # Check if openTool is defined
        is_defined = page.evaluate("typeof window.openTool === 'function'")
        print(f"openTool defined: {is_defined}")

        if not is_defined:
            # Maybe scripts are module type? They are not in index.html
            # <script src="js/expression.js"></script> ...
            # Wait a bit?
            time.sleep(1)
            is_defined = page.evaluate("typeof window.openTool === 'function'")
            print(f"openTool defined after wait: {is_defined}")

        if is_defined:
            # Force App Mode
            page.evaluate("document.body.classList.add('mobile-mode')")

            # 2. Verify Plotting
            page.evaluate("window.openTool('plotting')")
            page.wait_for_selector("#plot-mode", state="visible")
            page.select_option("#plot-mode", "parametric")
            page.screenshot(path="verification/plotting_tool.png")
            page.evaluate("window.closeTool()")

            # 3. Verify Geometry
            page.evaluate("window.openTool('geometry')")
            page.click("text=3D Objects")
            page.select_option("#geo-3d-shape", "cone")
            page.screenshot(path="verification/geometry_tool.png")
            page.evaluate("window.closeTool()")

            # 4. Statistics
            page.evaluate("window.openTool('statistics')")
            page.click("text=Distributions")
            page.select_option("#stat-dist-type", "poisson")
            page.screenshot(path="verification/statistics_tool.png")
            page.evaluate("window.closeTool()")

            # 5. Converter
            page.evaluate("window.openTool('converter')")
            page.select_option("#conv-category", "time")
            page.screenshot(path="verification/converter_tool.png")
            page.evaluate("window.closeTool()")

        browser.close()

if __name__ == "__main__":
    run()
