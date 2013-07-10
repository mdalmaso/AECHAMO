variable name in each column:
1. Time (cims time)
2. CIMS(H2SO4) (RC)
3. O3  (RC)
4. Temperature (RC) 
5. relative Humidity (RC) 
6. CPC (RC)
7. MCPC (=PSM) (RC)
8. Condensation SINK (RC) 
9. OH (RC)
10. Monoterpene mixing ratio in Reaction Chamber (RC)
11. Isoprene mixing ratio in Reaction Chamber (RC)
12. Monoterpene mixing ratio in Plant Chamber (PC)
13. Isoprene mixing ratio in plant Chamber (PC)
14. Flow from Plant into Reaction chamber [l/min] (PC->RC)
15. Conditioning flow O3 & H2O [l/min] (PC->RC)
16. Total dilution flow  (Reaction chamber outflow) [l/min] (RC->)

Note: all the data in this matrix were interpolated based on CIMS's time column, except for MCPC data. 
MCPC data have higher time resolution than CIMS. So, MCPC data are the average value between each record CIMS time period.

P.S: All the missing data were filled with NaNs.