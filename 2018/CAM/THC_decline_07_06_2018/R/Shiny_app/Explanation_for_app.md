&nbsp;&nbsp;&nbsp;This app fits a simple exponential decay model to an individual's CN-THCCOOH concentrations measured on different days, controlling for the number of times the person smoked in the last 30 days.
<br/><br/>
&nbsp;&nbsp;&nbsp;To use the app, in the bottom panel input (separated by commas) the CN-THCCOOH concentrations, the associated days since last ingestion, and the number of times the person smoked in the last 30 days. The app then fits a Bayesian regression model with conjugate priors to the data.
<br/><br/>
&nbsp;&nbsp;&nbsp;Once data has been inputted, the app displays the observed data (black points) against model predictions (solid line) with a 95% prediction interval (grey band). For comparison, predictions from the prior information (before fitting the model to data) are also shown (dashed line). Users can select the **Summary** tab to get a numeric summary with range of credible intervals for the estimated concentration of CN-THCCOOH at ingestion, the elimination rate, the influence of the level of use, the half-life (number of days to see a 50% reduction), and the window of detection (the last day before CN-THCCOOH concentrations become undetectable).
<br/><br/>
&nbsp;&nbsp;&nbsp;Advanced users can change the priors used in model estimation via the **Internal settings** tab. Priors are based on the data provided by 70 adolescent regular cannabis users from a study by conducted by the Center for Addiction Medicine at the Massachusetts General Hospital.
<br/><br/>
&nbsp;&nbsp;&nbsp;At least two observations are needed to estimate the slope. We hope this app is informative, but note that the efficacy of the model improves with more data and can still benefit from future data with larger, more diverse samples.
