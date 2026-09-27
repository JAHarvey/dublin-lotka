# Dublin–Lotka explorer 

This is an interactive teaching app for exploring how survival and fertility schedules set population growth. I built it for BES 550 Advanced Ecology at the University of Rhode Island, as part of the lecture on population demographics and its drivers.

Open the app: https://JAHarvey.github.io/dublin-lotka/

The first load takes 10–20 seconds because the app runs entirely in your browser.

What it does

You start from a life table: survivorship l(a) and fertility f(a) at each age a. The app calculates:

> R₀ = Σ l(a) f(a): expected lifetime reproduction
> Tc = Σ a l(a) f(a) / R₀: the cohort generation time
> r ≈ ln(R₀) / Tc: the Dublin–Lotka approximation of the intrinsic rate of increase
> Exact r: the solution to the Euler–Lotka equation, Σ e^(−ra) l(a) f(a) = 1
> λ = e^r, and the population's doubling or halving time

The plots show survivorship, the ages that offspring come from, where the Euler–Lotka equation is solved, and a 30-year projection.

How to use it
Pick a starting life table. There are three: the lecture example, a fast life history, and a slow, seabird-like life history.
Change adult survival, fertility, or the age at breeding with the sliders. You can also click a scenario: Cost of immunity, Earlier breeding, or Disease-driven shift.
Compare the baseline and modified values side by side.
To try your own data, enter it in the Edit life table tab.

The About tab has questions to work through.

> Notes and assumptions
> Ages are discrete. f(a) counts female offspring per female, and each table is treated as a single cohort.
> The fast and slow life tables and the scenario values are illustrative. They are not data from real populations.
> The approximation works best when r is near 0 and reproduction is concentrated around Tc. Comparing it to the exact r shows when it breaks down.
