# Lineare Gewöhnliche Differentialgleichungen – MATLAB Simulationen

**ETH Zürich – Elektrotechnik und Informationstechnologie (ITET)**  
Begleitprojekt zur Vorlesung *Analysis 1 & 2* (Ziltener)

---

## Übersicht

Dieses Repository enthält MATLAB-Simulationen zu **linearen gewöhnlichen Differentialgleichungen (GDG)** zweiter Ordnung mit direktem Bezug zur Mechanik und Elektrotechnik. Die Projekte visualisieren den Unterschied zwischen analytischer und numerischer Lösung und machen die verschiedenen Dämpfungsfälle direkt sichtbar.

| Projekt | Typ | DGL |
|---|---|---|
| `pendel_sim.m` | Homogen (freie Schwingung) | Gedämpftes Pendel |
| `schwingkreis_sim.m` | Inhomogen (erzwungene Schwingung) | RLC-Schwingkreis *(coming soon)* |

---

## Theorie: Lineare GDG 2. Ordnung

### Definition

Eine **gewöhnliche Differentialgleichung (GDG)** der Ordnung $n$ für eine Funktion $u: I \to \mathbb{R}$ ist eine Gleichung der Form:

$$\omega\big(t,\, u(t),\, \dot{u}(t),\, \ldots,\, u^{(n)}(t)\big) = 0, \quad t \in I$$

Eine GDG heisst **linear**, wenn sie die Form

$$u^{(n)} + a_{n-1}(t)\,u^{(n-1)} + \cdots + a_1(t)\,\dot{u} + a_0(t)\,u = b(t)$$

hat. Sie heisst **homogen**, wenn $b(t) = 0$, und **inhomogen**, wenn $b(t) \neq 0$.

### Die vollständige Lösung

Die allgemeine Lösung einer inhomogenen linearen GDG setzt sich immer aus zwei Teilen zusammen:

$$u(t) = \underbrace{u_{\text{hom}}(t)}_{\text{klingt ab}} + \underbrace{u_{\text{part}}(t)}_{\text{bleibt bestehen}}$$

- **Homogene Lösung** $u_{\text{hom}}$: Lösung der GDG mit rechter Seite = 0. Sie beschreibt das transiente Verhalten (Einschwingvorgang) und klingt bei Dämpfung ab.
- **Partikuläre Lösung** $u_{\text{part}}$: Eine spezielle Lösung der inhomogenen GDG. Sie beschreibt das stationäre (erzwungene) Verhalten.

---

## Projekt 1: Gedämpftes Pendel (`pendel_sim.m`)

### Physikalisches Modell

Das Pendel der Masse $m$, Länge $L$ in einem Medium mit Dämpfungskonstante $b$ folgt aus dem zweiten Newtonschen Gesetz und dem Stokesschen Reibungsgesetz:

$$m L^2 \ddot{\theta} + b L^2 \dot{\theta} + mgL\sin(\theta) = 0$$

Dividiert durch $mL^2$:

$$\ddot{\theta} + \frac{b}{m}\dot{\theta} + \frac{g}{L}\sin(\theta) = 0$$

Dies ist eine **nichtlineare** homogene GDG wegen des $\sin(\theta)$-Terms.

### Kleinwinkel-Näherung (Linearisierung)

Für kleine Winkel gilt $\sin(\theta) \approx \theta$, womit die GDG **linear** wird:

$$\ddot{\theta} + \frac{b}{m}\dot{\theta} + \frac{g}{L}\theta = 0$$

Rechte Seite = 0 → **homogene** lineare GDG zweiter Ordnung. Dies erlaubt eine geschlossene analytische Lösung über den Ansatz $\theta(t) = e^{\lambda t}$, der auf das **charakteristische Polynom** führt:

$$\lambda^2 + \frac{b}{m}\lambda + \frac{g}{L} = 0$$

mit den Wurzeln:

$$\lambda_{1,2} = -\gamma \pm \sqrt{\gamma^2 - \omega_0^2}$$

wobei $\gamma = \frac{b}{2m}$ die Dämpfungskonstante und $\omega_0 = \sqrt{\frac{g}{L}}$ die Eigenfrequenz ist.

### Die drei Dämpfungsfälle

Das Vorzeichen der Diskriminante $\gamma^2 - \omega_0^2$ bestimmt das qualitative Verhalten:

#### 1. Unterdämpfung: $\gamma < \omega_0$

Die Wurzeln $\lambda_{1,2}$ sind **komplex**, was zu einer gedämpften Schwingung führt. Die gedämpfte Eigenfrequenz ist:

$$\omega_D = \sqrt{\omega_0^2 - \gamma^2}$$

Die analytische Lösung lautet:

$$\theta(t) = e^{-\gamma t}\big(A\cos(\omega_D t) + B\sin(\omega_D t)\big)$$

Das Pendel schwingt mit abnehmender Amplitude aus. $A$ und $B$ werden aus den Anfangsbedingungen bestimmt:

$$A = \theta_0, \qquad B = \frac{\dot{\theta}_0 + \gamma\theta_0}{\omega_D}$$

#### 2. Kritische Dämpfung: $\gamma = \omega_0$

Die Diskriminante ist null, es gibt eine **doppelte reelle Wurzel** $\lambda = -\gamma$. Die Lösungsformel ändert sich grundlegend:

$$\theta(t) = (A + Bt)\,e^{-\gamma t}$$

Das Pendel kehrt **so schnell wie möglich** in die Ruhelage zurück, ohne zu schwingen. Dies ist der Übergangsfall zwischen Schwingung und kriechender Rückkehr.

> **Achtung im Code:** Bei $\omega_D = 0$ gibt es eine Division durch null in der Formel für $B$ der Unterdämpfung. Ein `if`-Sonderfall ist notwendig.

#### 3. Überdämpfung: $\gamma > \omega_0$

Die Wurzeln sind **zwei negative reelle Zahlen** $\lambda_1, \lambda_2 < 0$. Die Lösung ist eine Summe zweier Exponentialfunktionen:

$$\theta(t) = A e^{\lambda_1 t} + B e^{\lambda_2 t}$$

Das Pendel kriecht langsam in die Ruhelage, kein Schwingen. Langsamer als kritische Dämpfung.

### Analytisch vs. Numerisch

| Eigenschaft | Analytische Lösung | Numerische Lösung (ode45) |
|---|---|---|
| Gleichung | Linearisiert ($\sin\theta \approx \theta$) | Exakt ($\sin\theta$) |
| Gültigkeitsbereich | Nur kleine Winkel | Beliebige Winkel |
| Antriebskräfte | Nur bestimmte Formen | Beliebige Funktionen |
| Genauigkeit | Exakt (für das genäherte Modell) | Numerisch (sehr genau) |

Für kleine Startwinkel (z.B. $\theta_0 = 30°$) stimmen beide Lösungen sehr gut überein. Bei grossen Winkeln (z.B. $\theta_0 = 60°+$) weichen sie merklich ab – der Visualisierungseffekt macht die Grenzen der Kleinwinkel-Näherung direkt sichtbar.

### Numerische Methode: Zustandsraumdarstellung

`ode45` löst nur DGLs erster Ordnung der Form $\dot{z} = f(t, z)$. Die DGL zweiter Ordnung wird daher in zwei gekoppelte DGLs erster Ordnung umgeschrieben:

$$z_1 = \theta, \quad z_2 = \dot{\theta}$$

$$\dot{z}_1 = z_2$$
$$\dot{z}_2 = -\frac{b}{m}z_2 - \frac{g}{L}\sin(z_1)$$

Dies ist das Prinzip der **Zustandsraumdarstellung**, das in der Regelungstechnik und Elektrotechnik zentral ist. Allgemein lässt sich jede GDG $n$-ter Ordnung in $n$ gekoppelte GDGs erster Ordnung umschreiben.

### Parameter und ihre physikalische Bedeutung

| Parameter | Bedeutung | Analogie RLC |
|---|---|---|
| `m` – Masse [kg] | Trägheit des Systems | Induktivität $L$ [H] |
| `b` – Dämpfung [kg/s] | Reibung im Medium (Luft < Wasser < Öl) | Widerstand $R$ [Ω] |
| `L` – Länge [m] | Bestimmt $\omega_0$ | Kapazität $C$ [F] |
| `gamma` – Dämpfungskonstante | Wie schnell Einhüllende $e^{-\gamma t}$ abfällt | $\frac{R}{2L}$ |
| `omega0` – Eigenfrequenz | Schwingfrequenz ohne Dämpfung | $\frac{1}{\sqrt{LC}}$ |

---

## Projekt 2: RLC-Schwingkreis *(coming soon)*

Dieses Projekt simuliert einen **Reihenschwingkreis** (Widerstand $R$, Induktivität $L$, Kondensator $C$) mit einer angelegten Wechselspannungsquelle $u(t) = U_0\cos(\Omega t)$.

Die zugehörige GDG für die Stromstärke $I(t)$ ist **inhomogen**:

$$L\ddot{I} + R\dot{I} + \frac{1}{C}I = -U_0\Omega\sin(\Omega t)$$

oder in normierter Form:

$$\ddot{I} + \frac{R}{L}\dot{I} + \frac{1}{LC}I = b(t)$$

Dies ist dieselbe mathematische Struktur wie beim Pendel – mit rechter Seite $\neq 0$. Die vollständige Lösung ist:

$$I(t) = I_{\text{hom}}(t) + I_{\text{part}}(t)$$

Das Projekt zeigt insbesondere das Phänomen der **Resonanz** bei $\Omega = \omega_0 = \frac{1}{\sqrt{LC}}$, wo die Amplitude der partikulären Lösung maximal wird.

---

## Verwendung

```matlab
% Parameter anpassen
m = 1.0;       % Masse [kg]
L = 1.0;       % Länge [m]
b = 0.5;       % Dämpfung — probiere: 0 (ungedämpft), klein (unterdämpft),
               %            2*m*sqrt(g/L) (kritisch), gross (überdämpft)
theta0  = pi/6;   % Startwinkel [rad]
dtheta0 = 0;      % Startgeschwindigkeit [rad/s]

% Script ausführen
run('pendel_sim.m')
```

## Voraussetzungen

- MATLAB (getestet mit R2025a)
