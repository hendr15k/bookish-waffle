# Benutzeranleitung für Web CAS

Willkommen beim **Web CAS**, einem clientseitigen Computeralgebrasystem, das direkt in Ihrem Browser läuft. Diese Anwendung ermöglicht symbolische und numerische Berechnungen, grafische Darstellungen, Matrixoperationen und vieles mehr, ohne dass eine Internetverbindung oder eine Serverinstallation erforderlich ist.

## Inhaltsverzeichnis

1.  [Erste Schritte](#erste-schritte)
2.  [Benutzeroberfläche](#benutzeroberfläche)
3.  [Syntax und Eingabe](#syntax-und-eingabe)
4.  [Werkzeuge (Tools)](#werkzeuge)
    *   [Rechner (Calculator)](#rechner-calculator)
    *   [Algebra](#algebra)
    *   [Analysis (Calculus)](#analysis-calculus)
    *   [Plotten (Graphen)](#plotten-graphen)
    *   [Lineare Algebra](#lineare-algebra)
    *   [Statistik](#statistik)
    *   [Zahlen (Numbers)](#zahlen-numbers)
    *   [Einheitenumrechner](#einheitenumrechner)
    *   [Geometrie](#geometrie)
    *   [Finanzen](#finanzen)
    *   [Physik](#physik)
5.  [Befehlsreferenz](#befehlsreferenz)

---

## Erste Schritte

Da dieses Programm vollständig im Browser läuft, müssen Sie lediglich die Datei `index.html` in einem modernen Webbrowser (z. B. Chrome, Firefox, Safari, Edge) öffnen.

*   **Keine Installation:** Es ist keine Softwareinstallation notwendig.
*   **Datenschutz:** Alle Berechnungen finden lokal auf Ihrem Gerät statt.

## Benutzeroberfläche

Die Anwendung verfügt über zwei Hauptmodi:

1.  **Desktop-Modus (Standard):**
    *   **Eingabezeile:** Unten finden Sie eine Eingabeaufforderung für direkte Befehle (z. B. `solve(x^2-4, x)`).
    *   **Verlauf (History):** Der mittlere Bereich zeigt Ihre vergangenen Berechnungen und Ergebnisse in schön formatierter mathematischer Notation (LaTeX) an. Der Befehlsverlauf wird lokal im Browser gespeichert und bleibt für die Eingabe per Pfeiltasten nach einem Neuladen verfügbar.
    *   **Seitenleiste:** Rechts sehen Sie definierte Variablen und eine kurze Liste verfügbarer Funktionen.
    *   **App-Modus Button:** Oben rechts können Sie in den "App-Modus" wechseln.

2.  **App-Modus (Mobile/Tablet):**
    *   Klicken Sie auf **"App Mode"**, um eine Rasteransicht mit verschiedenen Werkzeugen (Apps) zu öffnen.
    *   Dieser Modus ist ideal für Touch-Geräte oder wenn Sie eine grafische Benutzeroberfläche für spezifische Aufgaben bevorzugen, anstatt Befehle zu tippen.

## Syntax und Eingabe

Sie können mathematische Ausdrücke wie in den meisten Taschenrechnern eingeben:

*   **Grundrechenarten:** `+`, `-`, `*`, `/`, `^` (Potenz).
*   **Klammern:** `(` und `)` zur Gruppierung.
*   **Funktionen:** `sin(x)`, `cos(x)`, `sqrt(x)` (Wurzel), `ln(x)`, etc.
*   **Variablen:** Sie können Variablen definieren, z. B. `a := 5` oder `f(x) := x^2 + 1`.

**Beispiele:**
*   `2 * (3 + 4)`
*   `sin(pi / 2)`
*   `x^2 + 2*x + 1`

---

## Werkzeuge

Im **App-Modus** stehen Ihnen folgende spezialisierte Werkzeuge zur Verfügung:

### Rechner (Calculator)
Ein wissenschaftlicher Taschenrechner mit erweitertem Funktionsumfang.
*   Tippen Sie Ihre Formel ein.
*   Nutzen Sie die Buttons für komplexe Operationen wie `sin`, `cos`, `Integral` etc.
*   **Calculate:** Berechnet das Ergebnis symbolisch oder numerisch.
*   **Plot:** Zeichnet den Graphen der eingegebenen Funktion.

### Algebra
Werkzeuge zum Umformen von Ausdrücken und Lösen von Gleichungen.
*   **Simplify & Expand:**
    *   *Simplify:* Vereinfacht einen Ausdruck (z. B. `2x + 3x` -> `5x`).
    *   *Expand:* Multipliziert Klammern aus (z. B. `(x+1)^2` -> `x^2+2x+1`).
*   **Solve:**
    *   Löst Gleichungen nach einer Variablen auf (z. B. `x^2 - 4 = 0` nach `x`).

### Analysis (Calculus)
Funktionen für die Differential- und Integralrechnung.
*   **Diff & Limit:**
    *   *Differentiate:* Berechnet die Ableitung einer Funktion.
    *   *Limit:* Berechnet den Grenzwert (Limes) an einer Stelle.
*   **Integration:**
    *   *Indefinite:* Berechnet das unbestimmte Integral (Stammfunktion).
    *   *Definite:* Berechnet das bestimmte Integral mit Unter- und Obergrenze.
*   **Taylor Series:**
    *   Entwickelt eine Funktion in eine Taylorreihe an einem bestimmten Entwicklungspunkt bis zur angegebenen Ordnung.

### Plotten (Graphen)
Ein leistungsstarker Funktionsplotter.
*   Geben Sie die Funktion `y=` ein (z. B. `x^2`).
*   Legen Sie den Bereich für die x-Achse fest (Min X, Max X).
*   Der Graph ist **interaktiv**: Sie können zoomen (Mausrad oder Pinch-Geste) und den Ausschnitt verschieben.

### Lineare Algebra
Arbeiten Sie mit Vektoren und Matrizen.
*   **Matrix Builder:** Erstellen Sie Matrizen visuell über ein Raster.
*   **Operations:** Berechnen Sie Inverse, Determinante, Spur (Trace) oder Transponierte einer Matrix.
*   **Vectors:** Berechnen Sie das Kreuzprodukt zweier Vektoren.

### Statistik
*   **Descriptive:** Berechnen Sie Mittelwert, Varianz, Median und Standardabweichung einer Datenliste (z. B. `[1, 2, 3, 4]`).
*   **Regression:** Führen Sie eine Lineare Regression durch. Geben Sie Punkte als Liste von Listen ein (z. B. `[[1,1], [2,3], [3,5]]`). Das Tool berechnet die Ausgleichsgerade und zeigt den Plot an.
*   **Distributions:** Ein interaktiver Visualisierer für die Normalverteilung, bei dem Sie Mittelwert ($\mu$) und Standardabweichung ($\sigma$) per Schieberegler anpassen können.

### Zahlen (Numbers)
Zahlentheoretische Funktionen.
*   Berechnungen wie Fakultät (`n!`), Primfaktorzerlegung (`Factor`), `GCD` (ggT), `LCM` (kgV) und Binomialkoeffizienten (`nCr`).
*   **Base Converter:** Rechnen Sie Zahlen in Echtzeit zwischen Dezimal-, Binär- und Hexadezimalsystem um.

### Einheitenumrechner
Konvertieren Sie Werte zwischen verschiedenen Einheiten:
*   **Länge:** Meter, Kilometer, Fuß, Meilen, etc.
*   **Masse:** Kilogramm, Gramm, Pfund, Unzen, etc.
*   **Temperatur:** Celsius, Fahrenheit, Kelvin.

### Geometrie
Berechnen Sie Flächen, Umfänge und Volumen.
*   **2D:** Kreis, Dreieck, Rechteck.
*   **3D:** Kugel, Zylinder, Würfel.

### Finanzen
*   **Compound Interest:** Zinseszinsrechnung.
*   **Loan Payment:** Berechnung von Kreditraten (Annuität).

### Physik
Berechnungen für Kinematik, Dynamik und Energie.
*   **Kinematik:** Lösen von SUVAT-Gleichungen (z. B. $v = u + at$, $s = ut + \frac{1}{2}at^2$).
*   **Dynamik:** Newtonsche Gesetze (Kraft $F=ma$, Masse, Beschleunigung).
*   **Energie:** Berechnung von Kinetischer Energie ($E_k$) und Potenzieller Energie ($E_p$).

---

## Befehlsreferenz

Diese Befehle können Sie direkt im "Calculator" oder in der Eingabezeile des Desktop-Modus verwenden.

### Grundfunktionen
*   `simplify(expr)`: Vereinfachen
*   `expand(expr)`: Ausmultiplizieren
*   `factor(n)`: Faktorisieren (Zahlen)
*   `N(expr)` oder `approx(expr)`: Numerische Auswertung

### Analysis
*   `diff(f, x)`: Ableitung von $f$ nach $x$.
*   `diff(f, x, n)`: $n$-te Ableitung.
*   `integrate(f, x)`: Unbestimmtes Integral.
*   `integrate(f, x, a, b)`: Bestimmtes Integral von $a$ bis $b$.
*   `limit(f, x, a)`: Grenzwert von $f$ für $x \to a$.
*   `taylor(f, x, a, n)`: Taylorreihe um $a$ bis Ordnung $n$.
*   `sum(expr, k, start, end)`: Summe ($\sum$).
*   `product(expr, k, start, end)`: Produkt ($\prod$).

### Lösen
*   `solve(eq, x)`: Löst eine Gleichung nach $x$.
*   `desolve(eq, y)`: Löst einfache Differentialgleichungen.
*   `desolve([eq1, eq2], [y1, y2])`: Löst lineare homogene Differentialgleichungssysteme mit konstanten Koeffizienten (z.B. `desolve([y'=z, z'=-y], [y, z])`).

### Lineare Algebra
Matrizen werden als Listen von Listen geschrieben: `[[1,2],[3,4]]`.
*   `det(M)`: Determinante.
*   `inv(M)`: Inverse Matrix.
*   `trans(M)`: Transponierte.
*   `cross(u, v)`: Kreuzprodukt.
*   `dot(u, v)`: Skalarprodukt.
*   `eigenvals(M)`: Eigenwerte.
*   `eigenvects(M)`: Eigenvektoren.
*   `rref(M)`: Zeilenstufenform (Gauß).

### Plotten
*   `plot(f, x, min, max)`: Zeichnet $f(x)$.
*   `plotparam([x(t), y(t)], t, min, max)`: Parametrischer Plot.
*   `plotpolar(r(t), t, min, max)`: Polarplot.

### Statistik & Listen
*   `mean(L)`, `median(L)`, `var(L)`, `std(L)`: Statistische Kennzahlen für Liste $L$.
*   `sort(L)`, `reverse(L)`: Sortieren und Umkehren.
*   `seq(expr, var, start, end, step)`: Erzeugt eine Sequenz.

### Sonstiges
*   `clear()`: Löscht den Verlauf und alle Variablen.
*   `help()`: Zeigt eine Liste von Befehlen an.
