**QUEST** is a tool for automated controller synthesis for incrementally input-to-state stable nonlinear control systems with respect to safety and reachability specifications.
The tool is implemented in C++ and contains two major parts:
- Construction of symbolic abstraction: the tool uses state-space quantization-free approach for construction of symbolic which helps to resolve the issue of so-called curse of dimensionality while modelling systems with high-dimensional state spaces.
- Symbolic controller synthesis: the synthesis of controller is implemented using fixed point computations.

### Requirments:
- A working C/C++ development environment
- A working installation of the CUDD library by Fabio Somenzi, which can be downloaded at http://vlsi.colorado.edu/~fabio/ with
	- the C++ object-oriented wrapper
	- the dddmp library and
	- the shared library

option enabled. The package follows the usual configure, make, and make install installation routine. We use cudd-3.0.0, with the configuration

	`$ ./configure --enable-shared --enable-obj --enable-dddmp--prefix=/opt/local/`
	
On Windows and linux, we experienced that the header files util.h and config.h were
missing in /opt/local and we manually copied them to 
	/opt/local/include.

For further details about windows installation, please refer to the installation_notes_windows.txt

FOR MORE DETAILS KINDLY REFER TO USER MANUAL AND REFERENCES THEREIN!!!
