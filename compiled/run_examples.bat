set root="."
set directory="def_direct.txt"

REM This will run no intervention
set name="output/run_no_intervention.txt"
.\compiled_model.exe %root% %directory% %name% strategy 0

REM This will run IPTp
set name="output/run_IPTp.txt"
.\compiled_model.exe %root% %directory% %name% strategy 1

REM This will run ISTp 
set name="output/run_ISTp.txt"
.\compiled_model.exe %root% %directory% %name% strategy 2

REM This will run ISTp with EIR=1
set name="output/run_ISTp_EIR1.txt"
.\compiled_model.exe %root% %directory% %name% strategy 2 EIR 1

REM This will run hybrid strategy 3 with a test in the first trimester with non-pregnant adult levels of sensitivity
set name="output/run_strategy_3_first_tri_test.txt"
.\compiled_model.exe %root% %directory% %name% strategy 3 first_tri_visit 1  first_tri_rdt 1

REM This will run hybrid strategy 4 with a test in the first trimester all tests with 90% sensitivity
set name="output/run_strategy_4_first_tri_perfect_test.txt"
.\compiled_model.exe %root% %directory% %name% strategy 4 perfect_test 1 first_tri_visit 1 first_tri_rdt 2

