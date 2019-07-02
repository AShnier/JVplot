# JVplot
Analyse current-voltage characteristics of solar cells

# Adam is not responsible for any calculation errors or local minima this script produces.
# In using this program you agree that you are responsible for checking that your data is consistant with any output this program produces.
Data requirements
⦁	IV plot is currently set to only work properly with negative Isc and positive Voc values.
⦁	If the recorded Isc is less than 100 nA the performance parameters are not calculated.
Zoomed graph
⦁	Imin, Imax, Vmin and Vmax are scaling values using Isc and Voc where Imin and Vmin are forced negative and Imax and Vmax forced positive. This is only implemented when bLim is set to True.
Series and shunt resistances
⦁	Makes data sets based on pRs and pRsh from the config file, (eg. 2xpRs + 1)
⦁	Runs a linear regression on voltage vs current (Ohms law: R=V/I)
Current input/output
⦁	There is a multiplier for the user to convert the current from their input files to amps “convert_to_amps” this is found in the config file
⦁	Similarly “r_current_out” is to convert the output from amps to the desired scale (mA default). This is accompanied by the unit name “s_current_out” which is by default the letters “mA”
⦁	“r_power_out” and “s_power_out” follow the previous point.
Max power
⦁	The maximum power is the largest average of two adjacent data points in the range between the Voc and Isc


Run Process flow
Build data arrays in a loop file by files
	Reads a data files to an array
	Fix truncation error on last data point
		Possible issue: Sets this value to 1 volt
	(optional) Divides current by area
	Evaluates data
	(optional) plots graphs
Makes excel file
Updates config file
Opens working directory
