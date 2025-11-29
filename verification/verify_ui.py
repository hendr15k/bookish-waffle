
from playwright.sync_api import sync_playwright
import os

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()

    # Load local file
    file_path = os.path.abspath("index.html")
    page.goto(f"file://{file_path}")

    # 1. Verify Labels have 'for' attributes matching inputs
    print("Verifying Labels...")
    inputs = page.locator("input.form-input").all()
    for input_el in inputs:
        input_id = input_el.get_attribute("id")
        if input_id:
            # Find a label with for=input_id
            label = page.locator(f"label[for='{input_id}']")
            if label.count() > 0:
                print(f"[PASS] Label found for #{input_id}")
            else:
                # Check if it is inside a label (implicit) or if I missed one
                # But my fix was explicit 'for' attributes.
                # Some inputs might not have labels (like matrix cells), but main tool inputs should.
                # Matrix cells don't have IDs usually or are generated.
                # The main inputs like 'calc-expr' should have it.
                if "calc-" in input_id or "alg-" in input_id or "int-" in input_id:
                     print(f"[FAIL] Label missing for #{input_id}")

    # 2. Verify Matrix Builder (generateMatrixGrid)
    print("Verifying Matrix Builder...")
    # Open Linear Algebra tool
    page.evaluate("openTool('linear')")

    # Click Generate Grid
    page.get_by_text("Generate Grid").click()

    # Check if inputs appeared
    matrix_inputs = page.locator("#matrix-builder-container input")
    if matrix_inputs.count() == 4: # 2x2 default
        print("[PASS] Matrix grid generated (2x2)")
    else:
        print(f"[FAIL] Matrix grid generation failed. Found {matrix_inputs.count()} inputs.")

    # 3. Verify Insert Matrix
    print("Verifying Insert Matrix...")
    # Fill first cell
    matrix_inputs.first.fill("5")

    # Click Insert Matrix
    page.get_by_text("Insert Matrix").click()

    # Check input value
    lin_expr = page.locator("#lin-expr")
    val = lin_expr.input_value()
    if val == "[[5,0],[0,0]]":
        print("[PASS] Matrix inserted correctly.")
    else:
        print(f"[FAIL] Matrix insertion incorrect. Got '{val}'")

    # 4. Verify Error Handling (UX)
    print("Verifying Error Handling...")
    # Open Plotting
    page.evaluate("openTool('plotting')")
    # Clear input to force error
    page.fill("#plot-expr", "")
    # Click Plot
    page.locator("#tool-plotting button").filter(has_text="Plot").click()

    # Check history for error
    last_result = page.locator(".result.error").last
    if last_result.is_visible() and "Please fill in" in last_result.inner_text():
        print("[PASS] Error displayed in history instead of alert.")
    else:
        print("[FAIL] Error not displayed in history.")

    # Screenshot
    page.screenshot(path="verification/ui_verification.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
