# JVplot <br>
Analyse current-voltage characteristics of solar cells <br>
<br>
Note from author: Hi I've been using this script since about 2017 and have a few people at Wits university (South Africa) who also rely on it to save them time. The input file we have been using is .lvm. I have uploaded this project as is and plan to generalise a lot of the code to allow a greater variety of input files when I have a time. Please edit the code or send a sample file if you'd like me to add compatibility, that will be a good reason for me to get that done sooner. <br>
<br>
Adam is not responsible for any calculation errors or local minima this script produces. <br>
In using this program you agree that you are responsible for checking that your data is consistant with any output this program produces. <br>
Data requirements <br>
⦁	IV plot is currently set to only work properly with negative Isc and positive Voc values. <br>
⦁	If the recorded Isc is less than 100 nA the performance parameters are not calculated. <br>
Zoomed graph <br>
⦁	Imin, Imax, Vmin and Vmax are scaling values using Isc and Voc where Imin and Vmin are forced negative and Imax and Vmax forced positive. This is only implemented when bLim is set to True. <br>
Series and shunt resistances <br>
⦁	Makes data sets based on pRs and pRsh from the config file, (eg. 2xpRs + 1) <br>
⦁	Runs a linear regression on voltage vs current (Ohms law: R=V/I) <br>
Current input/output <br>
⦁	There is a multiplier for the user to convert the current from their input files to amps “convert_to_amps” this is found in the config file <br>
⦁	Similarly “r_current_out” is to convert the output from amps to the desired scale (mA default). This is accompanied by the unit name “s_current_out” which is by default the letters “mA” <br>
⦁	“r_power_out” and “s_power_out” follow the previous point. <br>
Max power <br>
⦁	The maximum power is the largest average of two adjacent data points in the range between the Voc and Isc <br>
<br>

Run Process flow <br>
Build data arrays in a loop file by files <br>
	Reads a data files to an array <br>
	Fix truncation error on last data point <br>
		Possible issue: Sets this value to 1 volt <br>
	(optional) Divides current by area <br>
	Evaluates data <br>
	(optional) plots graphs <br>
Makes excel file <br>
Updates config file <br>
Opens working directory <br>
