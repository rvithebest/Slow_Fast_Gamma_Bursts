# Differential Burst dynamics of Slow and Fast gamma rhythms in Macaque V1
- Computed burst durations of slow and fast gamma (stimulus induced) produced in LFP recorded from V1 of two awake adult feamle monkeys
- Employed MP, OMP-GEAR, and HT to estimate burst durations
- Measured functional connectivity (PLV and WPLI) across all pair of electrodes
- Used WC network model (firing-rate) to simulate bursty slow and fast gamma rhythms
- For further details refer to: https://doi.org/10.1101/2025.10.01.679813
## Results
- Slow gamma rhythm exhibited significantly longer burst durations and longer latencies as compared to fast gamma.
- Slow gamma exhibited higher long-range synchrony compared to fast gamma.
- Results were consistent across different burst estimators.
- Simply changing the firing-rate time-constant of the corresponding inhibitory interneuronal population in the noisy WC model, leads to both slower and longer bursts.
## Requirements
- MATLAB R2024b (MathWorks, RRID: SCR_001622)
- Chronux toolbox (Mitra and Bokil, 2008) (http://chronux.org/, RRID: SCR_005547)
- Fieldtrip Toolbox ((Oostenveld et al., 2011), RRID: SCR_004849)
## Code Structure
- Burst_detection_methods: Codes to estimate burst dynamics of slow and fast gamma, using MP, OMP-GEAR and HT.
- LFP_analysis: Codes for plotting the figures and complete analysis of estimated burst duration and latencies.
- Model_simulation: Used to simulate non-linear self oscillating WC model to generate slow and fast gamma rhythms
- Model_analysis: Analysis of simulated LFP dataset (from model)
## Data availability


