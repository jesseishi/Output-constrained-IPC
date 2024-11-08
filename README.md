# Output-constrained invididual pitch control

- src/Lib/Linearizations from Find Optimal Azimuth Offset repository -> add link.
- Naming conventions
    - U stands for wind speed
    - `xpy` stands for the decimal `x.y` (p stands for point). 
- In all my calcXYSignals, signals must be a string, cannot be a char, otherwise iteration won't work. I use argument validation to ensure this.