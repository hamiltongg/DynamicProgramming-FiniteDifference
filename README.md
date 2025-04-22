Heterogeneous-Agent Models in Asset Pricing: The Dynamic Programming Approach and Finite Difference Method (2025)

(1) MainCode: It is a function that solves the PDE of Ak (the equilibrium HJB equation for agent k) in an economy with two agents differing in their RRA. This function implements the Finite Difference method with implicit and upwind schemes.

(2) Risk-Sharing-Rule: This function solves the risk-sharing rule and obtains the optimal consumption for the more risk-averse agent. The MainCode calls this function.  

(3) PolicyFunctions: It uses the function MainCode to obtain Ak, optimal consumption, portfolio, wealth, and equilibrium asset prices.

(4) InspectingRRA-Heterogeneity: It solves a two-agent model with two calibrations: the baseline calibration and an economy with a higher RRA for the more risk-averse agent. This m-file plots the optimal quantities and asset prices for both economies.
