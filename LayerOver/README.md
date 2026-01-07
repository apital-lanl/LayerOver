# LayerOver
Distributed under BSD-3 license by Los Alamos National Labs (copyright assertion O#: O4929)

Library for handling gcode and producing novel metrics, STL files, and limited visualization. 


NOTES:

Mechanical Testing Data Parsing
	- Best results will be obtained from CSVs; XLSX parsing is not fully implemented
	- If a new mechanical testing machine is used and the column names for stress and strain are unique, variables in DataAnalysis.MechanicalData header should be updated
		-- "example_namerows" and "namerow_example_parsing_dict" have fields to check for appropriate stress/strain column names