...One-way coupled FNPT-RANS model for wave-structure interaction problems...
FNPT ----> IITM-FNPT2D
RANS ----> IITM-RANS3D

Code is currently setup for simulating the interaction between a H/d=0.6 solitary wave (d=1m) and a circular cylinder of diameter 1m in a 20m x 4m x 2m. 
The domain is discretized using a uniform mesh of 500 x 50 x 100 (NX x NY x NZ). 

ALL CODES ARE CURRENTLY SERIAL. The entire simulation involves sequential execution of FIVE different codes.  

Linux is recommended and necessary for fast execution of all FIVE codes. In addition, the executable "fnpt" located in "IITM-FNPT2D" will only run in Linux.  

Linux command for compiling any code in terminal: g++ <CODENAME> -O2
Linux command for running the compiled code in terminal: ./a.out
	Note: Prof. Sriram's FNPT code is already compiled as an executable in the folder "IITM-FNPT2D" so compilation is not required in this case. Execute in terminal using the command: ./fnpt [will only run in Linux] 

Steps to be followed for the SOLITON-CYLINDER INTERACTION simulation:

(1) Go to the folder: "PADDLE-SIGNAL" ; compile and run the code "SOLITARY_WAVE_PADDLE_DISP.CPP". After execution, the code will output "wavePaddle.dat" which is the piston-type wave-paddle displacement history required for soliton generation. Open "wavePaddle.dat" in a text editor, enable display of linenumbers, note the total number of lines in the file and write that number at the top of the file in a newline and save. Alternately "wavePaddle.dat" can also be created using Prof. Sriram's MATLAB code "WaveGen". 
		Note: In "SOLITARY_WAVE_PADDLE_DISP.CPP", the file "wavePaddle.dat" gets written in "a" or "append" mode. So if you make changes in "SOLITARY_WAVE_PADDLE_DISP.CPP" pls. delete the previous "wavePaddle.dat" before re-running the code.  

(2) Copy "wavePaddle.dat" to the folder: "IITM-FNPT2D" and run the compiled executable "fnpt". The simulation is setup by changing the parameters in "FEMINPUT.dat". It's currently setup for H/d=0.6 soliton generation. Once execution of "fnpt" is complete, focus on two files: 
	(a) "Output_bvout.dat" - check surface-elevation time series at different wave-probe locations as specified in "FEMINPUT.dat" ; 
	(b) "Output_PHIT.dat" - contains the Eulerian wave-field (pressure, velocity, acceleration etc.) which is to be input to "IITM-RANS3D" after suitable interpolation. 
		Note: Input time-series in terms of paddle displacement is mandatory for ANY type of wave to be generated using "IITM-FNPT2D". Always naming the time-series file "wavePaddle.dat" is advised. 
Critical: [only applicable to this H/d=0.6 soliton] - Given the strong non-linearity of the soliton, the "fnpt" run (see step (2) above) will PAUSE at: "RUNNING TIMESTEP = 5513". Do not force continue the run beyond this point. Press "CTRL+C" in the terminal to exit. 

(3) Copy "Output_PHIT.dat" to the folder: "NSE-INPUT-EXTRACTION". Compile and run the code "FNPT_to_NS.CPP". After execution is complete, focus on three files: 
	(a) "Interp_pf.dat" - interpolated spatio-temporal series of the pressure and volume-fraction field [COLLOCATED MESH] ; 
	(b) "Interp_Ustag.dat" - interpolated spatio-temporal series of the U-velocity field [STAGGERED] ;
	(c) "Interp_Wstag.dat" - interpolated spatio-temporal series of the W-velocity field [STAGGERED].

(4) Copy "Interp_pf.dat", "Interp_Ustag.dat" and "Interp_Wstag.dat" to the folder: "IITM-RANS3D". 

(5) Go to the folder: "SOLID-PATCH" ; compile and run the code "VOF_INIT.c". Once execution is complete, the code will output two files: 
	(a) "CYLINDER.dat" - contains the entire SOLID_VOF field for patching the cylinder within the RANS domain
	(b) "VOFINIT_tecplot.dat" - contains the same information as "CYLINDER.dat" but this file can be imported in TecPlot/ParaVIEW to check quality/correctness of the initial patch. 
		Note: Inclusion of the -O2 optimization flag during compilation is essential to ensure rapid execution of "VOF_INIT.c". 

(6) Copy "CYLINDER.dat" to the folder: "IITM-RANS3D". Compile and run "3-D_TWO-PHASE-NSE.CPP". Compilation using the -O2 flag is advisable for faster execution. Invoking the -O2 flag during compilation will flash some WARNINGS related to reading from the files "Interp_pf.dat", "Interp_Ustag.dat", "Interp_Wstag.dat" and "CYLINDER.dat" - just ignore these warnings. Import the output files: "000_1.dat", "000_2.dat"....etc. to TecPlot/ParaVIEW for analyzing the 3D-RANS results. 
	Note: (a) "3-D_TWO-PHASE-NSE.CPP" is the main code - all functions and variable declarations are available in "functions.h".
		  (b) The global RANS mesh "NX x NY x NZ" is defined in "functions.h" via the pre-processor directive: #define ; for intance: "#define NY 50". Change the global cell-count in "functions.h" as per requirement. 
		  (c) Everything else in "functions.h" is hard-coded - changing any function/variable-declaration other than the global cell-count (NX x NY x NZ) may result in unpredictable behavior/incorrect solution.
Critical: ENSURE CONSISTENCY IN NS MESH-SIZE ACROSS STEPS (2)-(6). Directly changing NX or NZ in step (6) will lead to a wrong result.  

Contact Dr. Shaswat Saincher for any queries/issues at: shaswat.saincher@gmail.com  